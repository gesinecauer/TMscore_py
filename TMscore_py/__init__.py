import re
import os
import os.path
import errno
import subprocess
import numpy as np
import pathlib
import tempfile


def round_sigfig(array, sigfig):
    array = np.asarray(array)
    array_positive = np.where(
        np.isfinite(array) & (array != 0), np.abs(array), 10 ** (sigfig - 1))
    mags = 10 ** (sigfig - 1 - np.floor(np.log10(array_positive)))
    return np.round(array * mags) / mags


def convert_input(struct, filename=None, verbose=False, space=0):
    try:
        struct = np.asarray(struct).reshape(-1, 3)
    except ValueError:
        raise ValueError("Must input 3D structure, of shape (n, 3).")
    if (~np.isfinite(struct)).sum() > 0:
        raise ValueError("Structure must not contain nan or inf.")

    struct_round = struct.copy().round(6)
    lt0 = struct_round < 0
    struct_round[lt0] = round_sigfig(struct_round[lt0], sigfig=7 - space)
    struct_round[~lt0] = round_sigfig(struct_round[~lt0], sigfig=8 - space)
    # struct_round = round_sigfig(struct, sigfig=7 - space)

    if struct_round.max() >= 100000000:
        raise ValueError("Maximum value of structure must be below 100000000.")
    elif struct_round.min() <= -10000000:
        raise ValueError("Minimum value of structure must be above -10000000.")

    struct_8char = np.array([x.rjust(8) for x in struct_round.astype(
        f'U{8 - space}').flatten()]).reshape(struct.shape)
    nchar = set([len(x) for x in struct_8char.ravel()])
    if nchar != {8}:
        raise Exception
    struct_to_string = '\n'.join(['C ' + ' '.join(row) for row in struct_8char])
    struct_to_string = f"{struct.shape[0]}\n\n" + struct_to_string

    if filename is not None:
        with open(filename, "w") as f:
            f.write(struct_to_string)
    if verbose:
        print(struct_to_string, flush=True)

    return struct_8char, struct_to_string


def parse_kabsch_output(output):
    try:
        output = output.decode('utf8')
    except (UnicodeDecodeError, AttributeError):
        pass

    re_line = r'[^\n]+\n'
    re_sub = r'^.*-------- rotation matrix to rotate '
    re_sub += f'{re_line}{re_line}({re_line}{re_line}{re_line}).*$'
    kabsch_str = re.sub(re_sub, r'\1', output, count=1, flags=re.DOTALL).strip()
    kabsch = np.array([
        list(map(float, line.split()[1:])) for line in kabsch_str.split('\n')])
    translation = kabsch[:, 0]
    rotation = kabsch[:, 1:].T
    return rotation, translation


def parse_scores_ouput(output):
    try:
        output = output.decode('utf8')
    except (UnicodeDecodeError, AttributeError):
        pass
    res = {k: None for k in [
        'tm_score', 'gdt_ts', 'gdt_ts_c1', 'gdt_ts_c2', 'gdt_ts_c3',
        'gdt_ts_c4', 'gdt_ha', 'gdt_ha_c1', 'gdt_ha_c2', 'gdt_ha_c3',
        'gdt_ha_c4', 'rmsd', 'maxsub']}
    for line in output.split("\n"):
        line = re.sub(r"\s\s+", " ", line)
        line = re.sub(r"\s+=", "=", line)
        x = line.split(' ')
        if x[0] == "TM-score=":
            res['tm_score'] = float(x[1])
        elif x[0] == "GDT-TS-score=":
            res['gdt_ts'] = float(x[1])
            res['gdt_ts_c1'] = float(x[2].split('=')[1])
            res['gdt_ts_c2'] = float(x[3].split('=')[1])
            res['gdt_ts_c3'] = float(x[4].split('=')[1])
            res['gdt_ts_c4'] = float(x[5].split('=')[1])
        elif x[0] == "GDT-HA-score=":
            res['gdt_ha'] = float(x[1])
            res['gdt_ha_c1'] = float(x[2].split('=')[1])
            res['gdt_ha_c2'] = float(x[3].split('=')[1])
            res['gdt_ha_c3'] = float(x[4].split('=')[1])
            res['gdt_ha_c4'] = float(x[5].split('=')[1])
        elif x[0] == "RMSD":
            res['rmsd'] = float(line.split('=')[1].strip())
        elif x[0] == "MaxSub-score=":
            res['maxsub'] = float(x[1])

    res_bad = [k for k, v in res.items() if v is None]
    if len(res_bad) > 0:
        raise Exception(
            f"Error - the following scores were not found:\n{res_bad}")

    res['gdt_ts_info'] = (
        res['gdt_ts_c1'], res['gdt_ts_c2'], res['gdt_ts_c3'], res['gdt_ts_c4'])
    res['gdt_ha_info'] = (
        res['gdt_ha_c1'], res['gdt_ha_c2'], res['gdt_ha_c3'], res['gdt_ha_c4'])
    res = {k: v for k, v in res.items() if not re.match(
        r'^gdt_[a-z][a-z]_c[1234]$', k)}

    return res


def _print_scores(scores, prefix=None):
    to_print = []
    if 'tm_score' in scores:
        to_print.append(f"TM-score:\t{scores['tm_score']:g}")
    if 'gdt_ts' in scores:
        gdt_ts_info = ', '.join([f"{x:g}" for x in scores['gdt_ts_info']])
        to_print.append(f"GDT_TS:  \t{scores['gdt_ts']:g}\t({gdt_ts_info})")
    if 'gdt_ha' in scores:
        gdt_ha_info = ', '.join([f"{x:g}" for x in scores['gdt_ha_info']])
        to_print.append(f"GDT_HA:  \t{scores['gdt_ha']:g}\t({gdt_ha_info})")
    if 'rmsd' in scores:
        to_print.append(f"RMSD:    \t{scores['rmsd']:g}")
    if 'maxsub' in scores:
        to_print.append(f"MaxSub:  \t{scores['maxsub']:g}")
    if prefix is not None:
        to_print = [f"{prefix}{x}" for x in to_print]
    if len(to_print) > 0:
        print('\n'.join(to_print), flush=True)


class TMscore():
    def __init__(self, path=None):
        self._setup_path(path)

        self.rmsd = None
        self.gdt_ts = None
        self.gdt_ts_info = None
        self.gdt_ha = None
        self.gdt_ha_info = None
        self.tm_score = None
        self.maxsub = None
        self.scores = None
        self.rotation = None
        self.translation = None

        self._output = None
        self._scores = None
        self._rotation = None
        self._translation = None
        self._best_tmscore = None

    def _setup_path(self, path):
        if path is None:
            path = os.path.join(os.path.dirname(__file__), "TMscore")

        if os.path.isfile(path):
            self.path = str(pathlib.Path(path).resolve())
        else:
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), path)

    def _run_TMscore(self, structures, check_mirror=True, args=None, i=0):
        if args is None:
            args = []
        if isinstance(args, str):
            args = args.split()

        if i == 2:
            cmd = [self.path] + structures
            for j in range(2):
                if structures[j].lower().endswith('.xyz'):
                    cmd.extend([f'-infmt{j + 1}', '2'])
                elif structures[j].lower().endswith('.spicker'):
                    cmd.extend([f'-infmt{j + 1}', '1'])
            cmd.extend(args)
            output = {'original': subprocess.check_output(cmd).decode('utf8')}
            if check_mirror:
                output['mirror'] = subprocess.check_output(
                    cmd + ['-mirror', '1']).decode('utf8')
            return output

        if not isinstance(structures[i], (str, np.ndarray)):
            raise ValueError(
                f"Argument for struct{i} not recognized (type="
                f"{type(structures[i])}):\n{structures[i]}")
        if isinstance(structures[i], str):
            structures[i] = str(pathlib.Path(structures[i]).resolve())
            if not os.path.isfile(structures[i]):
                raise FileNotFoundError(f"File {structures[i]} does not exist.")
            return self._run_TMscore(
                structures, check_mirror=check_mirror, args=args, i=i + 1)
        elif isinstance(structures[i], (np.ndarray, list)):
            with tempfile.TemporaryDirectory() as dir_tmp:
                struct_file_tmp = os.path.join(dir_tmp, f'struct{i + 1}.xyz')
                convert_input(structures[i], filename=struct_file_tmp)
                structures[i] = struct_file_tmp
                return self._run_TMscore(
                    structures, check_mirror=check_mirror, args=args, i=i + 1)

    def __call__(self, structX, structY, check_mirror=True, args=None,
                 parse_kabsch=True, verbose=False):

        # In TM-score, structY is regarded as native structure, meaning they
        # superpose structX onto structY.

        self._output = self._run_TMscore(
            [structX, structY], check_mirror=check_mirror,
            args=args)
        self._scores = {}
        for x in self._output.keys():
            self._scores[x] = parse_scores_ouput(self._output[x])

        # Find best scores (if check_mirror)
        scores = list(self._scores.values())
        self._best_tmscore = list(self._scores.keys())[
            np.argmax([x['tm_score'] for x in scores])]
        self.tm_score = np.max([x['tm_score'] for x in scores])
        idx_ts = np.argmax([x['gdt_ts'] for x in scores])
        self.gdt_ts = [x['gdt_ts'] for x in scores][idx_ts]
        self.gdt_ts_info = [x['gdt_ts_info'] for x in scores][idx_ts]
        idx_ha = np.argmax([x['gdt_ha'] for x in scores])
        self.gdt_ha = [x['gdt_ha'] for x in scores][idx_ha]
        self.gdt_ha_info = [x['gdt_ha_info'] for x in scores][idx_ha]
        self.rmsd = np.min([x['rmsd'] for x in scores])
        self.maxsub = np.max([x['maxsub'] for x in scores])
        if len(self._output) == 1:
            self.scores = self._scores['original']
        else:
            self.scores = {
                'tm_score': self.tm_score, 'gdt_ts': self.gdt_ts,
                'gdt_ha': self.gdt_ha, 'rmsd': self.rmsd, 'maxsub': self.maxsub,
                'gdt_ts_info': self.gdt_ts_info,
                'gdt_ha_info': self.gdt_ha_info}

        # Optionally print scores
        if verbose:
            _print_scores(self.scores)

        # Optionally parse TM-score rotation matrix & translation vector
        if parse_kabsch:
            if verbose:
                print("\n", flush=True)
            self.parse_kabsch(verbose=verbose)

    def parse_kabsch(self, verbose=False):
        if self._output is None:
            raise ValueError("Scores not yet computed.")

        self._rotation = {}
        self._translation = {}
        for x in self._output.keys():
            self._rotation[x], self._translation[x] = parse_kabsch_output(
                self._output[x])

        self.rotation = self._rotation[self._best_tmscore]
        self.translation = self._translation[self._best_tmscore]

        array2string = lambda x: " " + np.array2string(
            x, precision=3).replace("[", "").replace("]", "")
        if verbose:
            print("Rotation:\n" + array2string(self.rotation), flush=True)
            print("Translation:  " + array2string(self.translation), flush=True)

    def __str__(self):
        if self._output is None:
            raise ValueError("Scores not yet computed.")
        if len(self._output) == 1:
            return self._output['original']
        return '\n'.join([
            f"{'=' * 80}\n" + f" {k} ".center(
                80, '=') + f"\n{'=' * 80}\n{v}" for k, v in self._output.items()])
