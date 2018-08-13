from .Conversions import  yr2s


########################
# Capture  
########################
global c
global G
global V_dot
global V_cross
global V_gal
global u_0
global k
global RCross
global RCrit


c = 3.0e10                          # Speed of Light: cm/s
G = 6.674e-11 * 100**3 * (1000)**-1 # Newtons's Constant: cm^3/(g s^2) 
V_dot = 220.0e5/c                   # Velocity of Sun relative to galactic center
V_cross = 29.8e5/c                  # Velocity of Earth relative to Sun
V_gal = 550.0e5/c                   # Escape velocity of the Galaxy
u_0 = 245.0e5/c                     # 
k = 2.5                             # Dimensionless
RCross = 6.738e8                    # Radius of Earth: cm 
RCrit = 3.48e8                      # Core-Mantle separator: cm


##########################
# Used in Annihilation
##########################
global Gnat     # G-Newton in natural units
global rhoCross # Density at the center of Earth
global TCross   # Temperature at the center of Earth
global tauCross # Age of Earth

Gnat = 6.71e-39        # GeV^-2
rhoCross = 5.67e-17    # GeV^4
TCross = 4.9134e-10    # GeV
tauCross = yr2s(4.5e9) # sec
