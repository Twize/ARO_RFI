def decimate(data,decfactor=128,polars=1):
    """Decimates data by a given factor and returns an array of power as well as the power squared (for SK algorithm)
    
    Keyword arguments:
    data -- Array, the data set which you want to decimate
    decfactor -- Factor by which to decimate the data, should be a proper fraction of the size of the data (default = 128)
    polars -- The polarization which you want to decimate and return, options are 1, 2 or 'both' (default = 1)
    """
    print 'Beginning Data Decimation...'
    
    # Importing required modules
    import numpy as np
    
    # Defining variables
    dims = np.shape(data)
    nfreqs = dims[1] # Number of frequency channels
    timelen = dims[0]
    
    #Separating the polarizations and getting rid of extra dimension
    if len(np.shape(data)) == 2: #If data only consists of a single polarization, continue 
        data = np.transpose(data)
        
        power = (np.real(data)**2. + np.imag(data)**2.) 
        V4 = power**2 
        
        # Decimating data for greater clarity
        reduced = np.reshape(power, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        V4reduced = np.reshape(V4, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        
        print 'Finished Data Decimation.'
        return reduced,V4reduced
    
    elif polars == 1: 
        polar = np.squeeze(data[:,:,:1])
        polar = np.transpose(polar)
        
        # Getting power and power squared for calculation of kurtosis 
        power = (np.real(polar)**2. + np.imag(polar)**2.) 
        V4 = power**2 
        
        # Decimating data for greater clarity
        reduced = np.reshape(power, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        V4reduced = np.reshape(V4, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        
        print 'Finished Data Decimation.'
        return reduced,V4reduced
        
    elif polars == 2:
        polar = np.squeeze(data[:,:,-1])
        polar = np.transpose(polar)
        
        # Getting power and power squared for calculation of kurtosis 
        power = (np.real(polar)**2. + np.imag(polar)**2.) 
        V4 = power**2 
        
        # Decimating data for greater clarity
        reduced = np.reshape(power, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        V4reduced = np.reshape(V4, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        
        return reduced,V4reduced
        
    elif polars == 'both':
        polar1 = np.squeeze(data[:,:,:1])
        polar1 = np.transpose(polar1)
        polar2 = np.squeeze(data[:,:,-1])
        polar2 = np.transpose(polar2)
        
        # Getting power and power squared for calculation of kurtosis 
        power1 = (np.real(polar1)**2. + np.imag(polar1)**2.) 
        V41 = power1**2 
        power2 = (np.real(polar2)**2. + np.imag(polar2)**2.) 
        V42 = power2**2 
        
        # Decimating data for greater clarity
        reduced1 = np.reshape(power1, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        V4reduced1 = np.reshape(V41, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        reduced2 = np.reshape(power2, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        V4reduced2 = np.reshape(V42, (nfreqs,timelen/decfactor,decfactor)).sum(axis=2)
        
        return reduced1,V4reduced1,reduced2,V4reduced2