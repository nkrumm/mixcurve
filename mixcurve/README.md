
## Mixcurve (a python package)

Mixcurve is designed to interpret 1:1 mixing data from the coagulation lab. The data is of the form:


|         | PT | PTT | TT | Fibrinogen |
|---------|----|-----|----|------------|
| Patient | 14 |  40 | 23 |        400 |
| Mix     | 14 |  32 |    |            |
| Pool    | 13 |  29 |    |            |

### Installation

```shell

git clone git@gitlab.labmed.uw.edu:nkrumm/mixcurve.git
pip install .

```

### Usage

```python

import mixcurve
mc = MixingConverter("vitk", "curve_fit")

mc.level_to_pt(30) // converts factor level of 30% to approproate PT
mc.ptt_to_level(44) // converts PTT=44 to factor level

```