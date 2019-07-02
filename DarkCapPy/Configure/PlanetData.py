import numpy as np
import pandas as pd
from scipy import interpolate

from .Constants import *
from .AtomicData import *
from .Conversions import *


##########################
# Taken from: https://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
# This imports the file 'PREM500.csv' within the DarkCapPy package so the user doesn't have to.
import os
this_dir, this_filename = os.path.split(__file__)
# this_dir, this_filename = os.path.split(__file__)



##########################
# Earth radius and mass
##########################
# Planet_Path = os.path.join(this_dir1, "PREM500_Mod.csv")
# VelocityDist_Path = os.path.join(this_dir2, "EarthDMVelDist.csv")
# Planet_Radius = 6.371e8     # cm
# Planet_Mass   = 5.972e27    # grams
# Planet_Life   = yr2s(4.5e9) # 4.5 Gyr -> sec


##########################
# Sun radius and mass
##########################
Planet_Path   = os.path.join(this_dir, "struct_b16_agss09.csv")
Vel_Dist_Path = os.path.join(this_dir, "SunDMVelDist.csv")
Planet_Radius = 69.551e9    # cm 
Planet_Mass   = 1.989e33    # g
Planet_Life   = yr2s(4.5e9) # 4.5 Gyr -> sec

# Variables to be used in DarkPhoton.py
# 1). radius_List
# 2). deltaR_List
# 3). escVel2_List
# 4). element_List



########################################################
#                  Data Input                          #
########################################################


##########################
# Read in Planet Data
##########################
Planet_File = pd.read_csv(Planet_Path,  delim_whitespace=True, header = 8)
Vel_Dist_File = pd.read_csv(Vel_Dist_Path)


radius_List = Planet_File['Radius'] * Planet_Radius
enclosedMass_List = Planet_File['Mass'] * Planet_Mass
element_List = np.asarray(Planet_File.columns[6:-1])


assert len(radius_List) == len(enclosedMass_List), 'Lengths of radius list and enclosed mass list do not match'


##########################
# Shell Thickness
##########################

def deltaR_Func(radiusList):
    # Input is a list if radiii values,
    # output is a list of deltaR values
    # DeltaR = Next radius - current radius
    
    # We define the variable 'radiusListm1' in order to give a deltaRList which has the same length as radiusList
    radiusListm1 = radiusList[0:len(radiusList)-1]
    s = [0] # Temporary variable used to obtain deltaRList. Stores radiusList offset by one index.
    for i in radiusListm1:
        s.append(i)

    deltaRList = radiusList[0:len(radiusList)] - s[0:len(s)]
    return deltaRList

##########################
# Shell Mass
##########################

def shellMass_Func(totalMassList):
    shellMass_List = []
    bigMass = 0
    smallMass = 0
    for i in range(0,len(totalMassList)):
        if i == 0:
            shellMass_List.append(0)
        else:
            bigMass = totalMassList[i]
            smallMass = totalMassList[i-1]
            shellMass_List.append(bigMass-smallMass)
    return shellMass_List



##########################
# Shell Density
##########################

def shellDensity_Func(shellMass, shellRadius, deltaR):
	shellDensity = []
	for i in range(0,len(shellMass)):
		shellDensity.append(shellMass[i]/(4*np.pi*shellRadius[i]**2 * deltaR[i]))
	# Kludge for radius = 0     
	shellDensity[0] = shellDensity[1]
	return shellDensity


##########################
# Number Density of each element
##########################
def numDensity_Func(element):
    numDensityList = []
    for i in range(0,len(shellDensity_List)):
        mf = Planet_File[str(element)][i]
        numDensityList.append(mf * g2GeV(shellDensity_List[i]) / amu2GeV(atomicNumbers[element]))
    return numDensityList




##########################
# Escape Velocity
##########################

def escVel_Func(index, enclosedMassList, radiusList, deltaRList):
    
    G_Newton = 6.674e-11 * 100**3 * (1000)**-1 # in cm^3/(g s^2)
    c = 3e10 # in cm/s
    factor = 2.*G_Newton/c**2 # prefactor
    constant = max(enclosedMassList) / max(radiusList)
    
    assert len(enclosedMassList) == len(radiusList), 'Lengths of Enclosed mass list and radius list do not match' 
    assert len(radiusList) == len(deltaRList), 'Lengths of radius list and delta R list do not match'

    if (index == 0):
        tempSum = 0

    elif (index != 0):
        tempSum = 0    
        for i in range(index, len(radiusList)):
            summand = enclosedMassList[i] * deltaRList[i] / (radiusList[i])**2
            tempSum += summand

    return (factor * (tempSum + constant))



##########################
# Generate all lists
##########################
deltaR_List = deltaR_Func(radius_List)
shellMass_List = shellMass_Func(enclosedMass_List)
shellDensity_List = shellDensity_Func(shellMass_List, radius_List, deltaR_List)


assert len(radius_List) == len(deltaR_List)
assert len(radius_List) == len(shellMass_List)
assert len(radius_List) == len(shellDensity_List)

escVel2_List = []                            #| Construct an array of escape velocities 
for i in range(0,len(radius_List)):          #|
    escVel2_List.append(escVel_Func(i, enclosedMass_List, radius_List, deltaR_List))  #|
    
escVel2_List[0] = escVel2_List[1] # Set the i=0 and i=1 escape velocities equal 


##########################
# DM Velocity Distrubution
##########################
velocity_Range_List = Vel_Dist_File['Velocity_Range']        # A list of velocities between 0 and V_gal
planet_velocity_List   = Vel_Dist_File['VelocityDist_Planet_Frame'] # The DM velocity distrubution in the planet frame


########################
# Interpolate the Velocity Distribution
########################
velRange = velocity_Range_List
fCrossVect = planet_velocity_List
fCrossInterp = interpolate.interp1d(velRange, fCrossVect, kind ='linear')


##########################
# Interpolate
# These are intentionally commented out, we don't actually use them in DarkPhoton.py
# I'm not sure why I made these but they are here if they are usefull
##########################
# Earth_enclosedMassInterp = interpolate.interp1d(radius_List, enclosedMass_List, kind='linear') 
# Earth_escVel2Interp = interpolate.interp1d(radius_List, escVel2_List, kind='linear')           
# Earth_densityInterp = interpolate.interp1d(radius_List,Earth_density_List,kind='linear')


