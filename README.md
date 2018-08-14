# DarkCapPy

Written by: Adam Green (agree019@ucr.edu) and Philip Tanedo (flip.tanedo@ucr.edu)

DarkCapPy is a Python 3/Jupyter package for calculating rates associated with dark matter capture in the Earth, annihilation into light mediators, and the subsequent observable decay of the light mediators near the surface of the Earth. The package includes a calculation of the Sommerfeld enhancement at the center of the Earth and the timescale for capture--annihilation equilibrium. The code is open source and can be modified for other compact astronomical objects and mediator spins. The reference paper can be found at: https://arxiv.org/abs/1509.07525v3

## Example Usage

```python
import numpy as np
import pandas as pd
from scipy.interpolate import interpolate

from DarkCapPy import *
import DarkCapPy.DarkPhoton as DP

mx = 1000
ma = 1
epsilon = 1e-8
alpha = 1/137
alphax = 0.035
tauCross = DP.tauCross

cap1 = DP.cCap(mx,ma,epsilon,alpha,alphax)
kappa0 = DP.kappa_0(mx,alpha)
cap2 = DP.cCapQuick(mx,ma,epsilon,alphax,kappa0)

sommerfeld = DP.thermAvgSommerfeld(mx,ma,alphax)
sigma = DP.sigmaVtree(mx,ma,alphax)
ann = DP.cAnn(mx,sigma,sommerfeld)

tau = DP.tau(cap1,ann)

gammaAnn = DP.gammaAnn(cap1,ann)
L = DP.decayLength(mx,ma,epsilon,1)
Edecay = DP.epsilonDecay(L)

signal = DP.iceCubeSignal(gammaAnn,Edecay,DP.yr2s(10))


print ('Capture_1            :', cap1)
print ('Kappa_0              :', kappa0)
print ('Capture_2            :', cap2)
print ('Therm Avg Sommerfeld :', sommerfeld)
print ('Sigma V              :', sigma)
print ('Annihilation         :', ann)
print ('EQ Time              :', tau)
print ('Gamma_ann            :', gammaAnn)
print ('Decay Length         :', L)
print ('Epsilon_decay        :', Edecay)
print ('N_signal             :', signal)
```
