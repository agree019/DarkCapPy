import numpy as np
import pandas as pd
from scipy import interpolate

# from .Constants import *
from .AtomicData import *
from .Conversions import *


##########################
# Taken from: https://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
# This imports the file 'PREM500.csv' within the DarkCapPy package so the user doesn't have to.
import os
this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, "PREM500.csv")
##########################



##########################
# Read in Earth Data
##########################
data = pd.read_csv(DATA_PATH)

radiusTemp1 = data['Radius[m]']  # Radius in meters
densityTemp1 = data[['Density[kg/m^3]']] # Density in kg/m^3

# # The interpolation function doesn't like these objects, so they need to be massaged into 1-D numpy arrays
radiusListBadUnits = np.asarray(radiusTemp1).squeeze()
densityListBadUnits = np.asarray(densityTemp1).squeeze()


# Convert Units and trim off the zero-index value
# These are the lists we use for the rest of the computation
radiusList = radiusListBadUnits[1:] * 100 # cm
densityList = densityListBadUnits[1:] * (100)**-3 * 1000 #  in g/cm^3



##########################
# Define Planet Constants
# Updated on 4/2/19
# These two vales are imported into the Constants.py file
##########################
RCross_Planet = max(radiusList) # The radius of the capturing body
RCrit_Planet  = radiusList[int(np.floor(len(radiusList)/2))] # The radius which separates the mantle from the core


##########################
# Shell Thickness
##########################
# We define the variable 'radiusListm1' in order to give a deltaRList which has the same length as radiusList
radiusListm1 = radiusList[0:len(radiusList)-1]
s = [0] # Temporary variable used to obtain deltaRList. Stores radiusList offset by one index.
for i in radiusListm1:
    s.append(i)

deltaRList = radiusList[0:len(radiusList)] - s[0:len(s)]



##########################
# Shell Mass
##########################
shellMassList = []
lenRadiusList = range(0,len(radiusList))
for i in lenRadiusList:
    shellMassList.append(4 * np.pi * radiusList[i]**2 * densityList[i] * deltaRList[i])



##########################
# Enclosed Mass
##########################
enclosedMassList = []
tempMassSum = 0 # Temporary mass shell sum
for i in shellMassList:
    tempMassSum = tempMassSum + i # add value of shellMass to previous sum
    enclosedMassList.append(tempMassSum)


##########################
# Escape Velocity
##########################

def escVelFunc(index):
    '''
    escVelFunc(index)

    This is the exact same function as Accumulate in Mathematica
    
    returns the escape velocity at the specified index
    '''
    G_Newton = 6.674e-11 * 100**3 * (1000)**-1 # in cm^3/(g s^2)
    c = 3e10 # in cm/s
    factor = 2.*G_Newton/c**2 # prefactor
    constant = max(enclosedMassList) / max(radiusList)
    
    if (index == 0):
        tempSum = 0
        
    elif (index != 0):
        tempSum = 0    
        for i in range(index, len(radiusList)):
            summand = enclosedMassList[i] * deltaRList[i] / (radiusList[i])**2
            tempSum += summand
    
    return factor * (tempSum + constant)

escVel2List = []                       #| Construct an array of escape velocities 
for i in lenRadiusList:                #|
    escVel2List.append(escVelFunc(i))  #|
    
escVel2List[0] = escVel2List[1] # Set the i=0 and i=1 escape velocities equal 




##########################
# Number Density
##########################
mf = 0 # Mass Fraction

def numDensityList(element):
    '''
    numDensityList(element)

    See ElementList in AtomicData.py for valid elements

    Returns: the number density array of element
    '''
    numDensityList = []
    for i in lenRadiusList:
        
        if radiusList[i] < RCrit_Planet:
            mf = coreMassFrac[element]
            
        elif RCrit <= radiusList[i] <= RCross_Planet:
            mf = mantleMassFrac[element]
            
        elif radiusList[i] > RCross_Planet:
            mf = 0

        ##########################################################
        # [mf] = mass fraction = dimensionless
        # [densityList] = g/cm^3
        # [atomicNumbers] = dimensionless
        # To make life easy, we convert everything into GeV
        #    => densityList -> g2GeV(densitylist)
        #    => atomicNumbers -> amu2GeV(atomicNumbers)
        
        n_i = mf * g2GeV(densityList[i]) /(amu2GeV(atomicNumbers[element]))

        numDensityList.append(n_i)
        
    return numDensityList



##########################
# Interpolate
##########################
enclosedMassInterp = interpolate.interp1d(radiusList,enclosedMassList,kind='linear') # g
escVel2Interp = interpolate.interp1d(radiusList,escVel2List,kind='linear')           # Dimensionless
densityInterp = interpolate.interp1d(radiusList,densityList,kind='linear')           # g/cm^3