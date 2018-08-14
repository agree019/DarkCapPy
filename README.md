# DarkCapPy

Written by: Adam Green (agree019@ucr.edu) and Philip Tanedo (flip.tanedo@ucr.edu)

DarkCapPy is a Python 3/Jupyter package for calculating rates associated with dark matter capture in the Earth, annihilation into light mediators, and the subsequent observable decay of the light mediators near the surface of the Earth. The package includes a calculation of the Sommerfeld enhancement at the center of the Earth and the timescale for capture--annihilation equilibrium. The code is open source and can be modified for other compact astronomical objects and mediator spins. The reference paper can be found at: https://arxiv.org/abs/1509.07525v3

## Example Usage

```python
import numpy as np
import pandas as pd
from scipy.interpolate import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style
mpl.rcParams.update(mpl.rcParamsDefault)
%matplotlib inline
from datetime import datetime

from DarkCapPy import *
import DarkCapPy.DarkPhoton as DP

###################
# Define File Paths
###################
# These paths are specific to the preset folders in Temp;ate_Caluclation
def sommerfeldPath(file):
    path = 'Sommerfeld/' + file
    return path

def branchPath(file):
    path = 'Branching_Ratio/' + file
    return path

def signalPath(file):
    path = 'Signal/' + file
    return path
    
def signalBackupPath(file):
    path = 'Signal/Signal_Backups/' + file
    return path

print ('Complete')
```
