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
