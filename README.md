# TMscore_py
------------
A simple class to wrap the original TMscore C++ code.
The source code for TMscore (required) can be found here: https://zhanglab.ccmb.med.umich.edu/TM-score/
It is also available in this repository.

## Usage

``` Python
from TMscore_py import TMscore


tmscore = TMscore()
print(tmscore(structX, structY))  # file paths or numpy arrays
```

## License
See LICENSE. For more info on TMscore see https://zhanglab.ccmb.med.umich.edu/TM-score/


## A note about installation

Install via

```
python setup.py develop
```

If you run into the following error upon reinstall, remove the `_skbuild/` directory (see [here](https://github.com/equinor/ert/issues/2408) for details.)

```
CMake Error at CMakeLists.txt:3 (find_package):
  By not providing "FindPythonExtensions.cmake" in CMAKE_MODULE_PATH this
  project has asked CMake to find a package configuration file provided by
  "PythonExtensions", but CMake did not find one.

  Could not find a package configuration file provided by "PythonExtensions"
  with any of the following names:

    PythonExtensionsConfig.cmake
    pythonextensions-config.cmake

  Add the installation prefix of "PythonExtensions" to CMAKE_PREFIX_PATH or
  set "PythonExtensions_DIR" to a directory containing one of the above
  files.  If "PythonExtensions" provides a separate development package or
  SDK, be sure it has been installed.
```


