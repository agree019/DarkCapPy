########################
# FilePaths
########################

def photonSommerfeldPath(file):
    '''
    Path to dark photon Sommerfeld Data files
    '''
    path = 'DarkPhotonCapture/SommerfeldData/' + file
    return path

def photonEquilibriumPath(file):
    '''
    Path to dark photon Equilibrium time contour plots
    '''
    path = 'DarkPhotonCapture/EquilibriumPlots/' + file
    return path


def photonCapturePath(file):
    path = 'DarkPhotonCapture/CapturePlots/' + file
    return path

def photonSignalDataPath_Complete(file):
    '''
    Path to complete dark photon Signal Rate Data files
    '''
    path = 'DarkPhotonCapture/SignalRateData/Complete/' + file
    return path

def photonSignalDataPath_Incomplete(file):
    '''
    Path to incomplete dark photon Signal Rate Data files
    '''
    path = 'DarkPhotonCapture/SignalRateData/Incomplete/' + file
    return path

def photonSignalDataPath_Incomplete_Backups(file):
    '''
    Path to derk photon Signal Rate Bakup Files
    '''
    path  = 'DarkPhotonCapture/SignalRateData/Incomplete/Backups/' + file
    return path

def photonSignalPlotPath(file):
    '''
    Path to dark photon Signal Rate Plots
    '''
    path = 'DarkPhotonCapture/SignalRatePlots/' + file
    return path

def photonBranchPath(file):
    '''
    Path to dark photon Branching Ratio Data
    '''
    path = 'DarkPhotonCapture/BranchingRatioData/' + file
    return path