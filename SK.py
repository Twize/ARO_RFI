def SK(data,V4,M=128,upper=1.5,lower=0.6):
    """Returns a binary mask of flagged RFI using the Spectral Kurtosis (SK) algorithm
    
    Keyword arguments:
    data -- Input data set (in terms of power ie. |V|^2)
    V4 -- Data set squared (power squared), if the data-set has been decimated, this must also be decimated with the same factor 
    M -- Number of pixels summed together (default = 128)
    upper -- The upper limit on the spectral kurtosis for thresholding (value around ~1.5 for M=128 seems to work well)
    lower -- The lower limit on the spectral kurtosis for thresholding (value around ~0.6 for M=128 seems to work well)
    
    """
    print 'Beginning SK...'
    
    import numpy as np
    
    # Defining required variables
    dims = np.shape(data) #Dimensions of reduced array for later use
    eps = 1E-5 #small value to add to arrays to negate divide by zero errors
    
    # Creating mask with same dimensions as input data 
    kurtmask = np.zeros_like(data)
    
    #Differencing kurtosis and squared second moment
    Speckurt = ((M+1)/(M-1))*(M*(V4+eps)/((data**2)+eps)-1)
    
    #Loop to do the thresholding
    for y in range(dims[0]):
        for x in range(dims[1]):
            if Speckurt[y][x] > upper or Speckurt[y][x] < lower:
                kurtmask[y][x] = 1
            else:
                continue
                        
    percentflag = np.mean(kurtmask)*100 # Percent of data flagged by SK algorithm 
    print 'Percent flagged by SK algorithm: ', percentflag
    print 'Finished SK.'
    
    return kurtmask