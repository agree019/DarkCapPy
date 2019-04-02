from .Conversions import  yr2s
from .PlanetData  import RCross_Planet
from .PlanetData  import RCrit_Planet

########################
# Capture  
########################
global c       # Speed of Light in cm/s
global G       # Gravitational constant in cm^3/(g s^2) 
global V_dot   # Velocity of Sun relative to galactic center
global V_cross # Velocity of Earth relative to Sun
global V_gal   # Escape velocity of the Galaxy
global u_0     # Characteristic speed of Galactic DM
global k       # Kurtosis of velocity disribution 
global RCross  # Radius of Earth: cm 
global RCrit   # Core-Mantle separator: cm


c = 3.0e10                           # cm/s
G = 6.674e-11 * 100**3 * (1000)**-1  # cm^3/(g s^2)
V_dot = 220.0e5/c                    # Dimensionless
V_cross = 29.8e5/c                   # Dimensionless
V_gal = 550.0e5/c                    # Dimensionless
u_0 = 245.0e5/c                      # Dimensionless
k = 2.5                              # Dimensionless
RCross = RCross_Planet               # cm
RCrit = RCrit_Planet     			 # cm


##########################
# Used in Annihilation
##########################
global Gnat     # Gravitational constant G in natural units GeV^-2
global rhoCross # Density at the center of Earth in GeV^4
global TCross   # Temperature at the center of Earth in GeV
global tauCross # Age of Earth in seconds

Gnat = 6.71e-39        # GeV^-2
rhoCross = 5.67e-17    # GeV^4
TCross = 4.9134e-10    # GeV
tauCross = yr2s(4.5e9) # sec
