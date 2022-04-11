import pytest
from TMscore_py import TMscore
import numpy as np


def test_init():
    tmscore = TMscore()
    assert tmscore is not None


def test_init_wrong_path():
    with pytest.raises(FileNotFoundError):
        TMscore('this path is bananas')


def test_pdb():
    tmscore = TMscore()
    tmscore("./proteins/1crn.pdb", "./proteins/best.pdb", check_mirror=False)

    assert tmscore.gdt_ha == 0.2772
    assert tmscore.gdt_ha_info == (0.1304, 0.1957, 0.2609, 0.5217)
    assert tmscore.gdt_ts == 0.4402
    assert tmscore.gdt_ts_info == (0.1957, 0.2609, 0.5217, 0.7826)
    assert tmscore.maxsub == 0.3286
    assert tmscore.tm_score == 0.2875
    assert tmscore.rmsd == 6.709


def test_pdb_reflected_check_mirror():
    tmscore = TMscore()
    tmscore("./proteins/best.pdb", "./proteins/best.mirror.pdb",
            check_mirror=True)

    assert tmscore.gdt_ha == 1
    assert tmscore.gdt_ha_info == (1, 1, 1, 1)
    assert tmscore.gdt_ts == 1
    assert tmscore.gdt_ts_info == (1, 1, 1, 1)
    assert tmscore.maxsub == 1
    assert tmscore.tm_score == 1
    assert tmscore.rmsd == 0


def test_xyz_reflected_check_mirror():
    tmscore = TMscore()
    tmscore("./proteins/small_prot.xyz", "./proteins/small_prot.mirror.xyz",
            check_mirror=True)

    assert tmscore.gdt_ha == 1
    assert tmscore.gdt_ha_info == (1, 1, 1, 1)
    assert tmscore.gdt_ts == 1
    assert tmscore.gdt_ts_info == (1, 1, 1, 1)
    assert tmscore.maxsub == 1
    assert tmscore.tm_score == 1
    assert tmscore.rmsd == 0


def test_array_reflected_dont_check_mirror():
    np.random.seed(0)
    structX = np.random.normal(0, 1, size=(10, 3))
    structY = structX.copy()
    structY[:, 1] = -structY[:, 1]

    tmscore = TMscore()
    tmscore(structX, structY, check_mirror=False)

    assert tmscore.gdt_ha == 0.8250
    assert tmscore.gdt_ha_info == (0.6, 0.8, 0.9, 1)
    assert tmscore.gdt_ts == 0.9250
    assert tmscore.gdt_ts_info == (0.8, 0.9, 1, 1)
    assert tmscore.maxsub == 0.8703
    assert tmscore.tm_score == 0.6043
    assert tmscore.rmsd == 1.540


def test_array_reflected_check_mirror():
    np.random.seed(0)
    structX = np.random.normal(0, 1, size=(10, 3))
    structY = structX.copy()
    structY[:, 1] = -structY[:, 1]

    tmscore = TMscore()
    tmscore(structX, structY, check_mirror=True)

    assert tmscore.gdt_ha == 1
    assert tmscore.gdt_ha_info == (1, 1, 1, 1)
    assert tmscore.gdt_ts == 1
    assert tmscore.gdt_ts_info == (1, 1, 1, 1)
    assert tmscore.maxsub == 1
    assert tmscore.tm_score == 1
    assert tmscore.rmsd == 0


def test_array_rotated():
    np.random.seed(0)
    structX = np.random.normal(0, 1, size=(10, 3))
    R = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]])
    structY = np.dot(structX, R)

    tmscore = TMscore()
    tmscore(structX, structY, check_mirror=False)

    assert tmscore.gdt_ha == 1
    assert tmscore.gdt_ha_info == (1, 1, 1, 1)
    assert tmscore.gdt_ts == 1
    assert tmscore.gdt_ts_info == (1, 1, 1, 1)
    assert tmscore.maxsub == 1
    assert tmscore.tm_score == 1
    assert tmscore.rmsd == 0
