########################
#  Unit Conversions
########################

def amu2GeV(par1):
    '''
    1 = 0.938272 GeV / Amu
    '''
    return 0.938272 * par1

def amu2g(par1):
    '''
    1 = 1.66053892e-24 g / Amu
    '''
    return 1.66053892e-24 * par1 

def GeV2s(par1):
    '''
    1 = 6.58e-25 GeV Sec
    '''
    return 1.52e24 * par1  

def s2GeV(par1):
    '''
    1 = 6.58e-25 GeV Sec
    '''
    return (6.58e-25)**-1 * par1

def GeV2cm(par1):
    '''
    1 = 1.97e-14 GeV cm
    '''
    return 5.06e13 * par1

def cm2GeV(par1):
    '''
    1 = 1.97e-14 GeV cm
    '''
    return (0.197e-13)**-1 * par1 


def KeltoGeV(par1):
    '''
    1 = 1.16e13 GeV / K
    '''
    return 8.62e-14 * par1 # GeV

def yr2s(par1):
    '''
    1 = (365)(24)(3600) sec/yr
    '''
    return (3.1536000e7) * par1 # s


def g2GeV(par1):
    '''
    1 = 5.62e23 GeV / g
    '''
    return (5.62e23) * par1 #GeV
