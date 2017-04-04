def SumThreshold(data,axis='both'):
    """Returns a binary mask of flagged RFI using A. Offringa's SumThreshold method. 
    
    Uses exponentially increasing subsequence lengths (2**i) to speed up performance. Thresholds are calculated based on a shifted and scaled Chi-squared distribution.
    
    Keyword arguments:
    data -- Input data set (in terms of power ie. |V|^2)
    axis -- Which axis to perform the flagging along, axis = 0 is the frequency-axis, axis = 1 is the time-axis
            and 'both' will run along each axis independently, then combine the masks (default axis = 'both')
    
    """
    print 'Beginning SumThreshold...'
    
    # Importing required modules
    import numpy as np
    from scipy.stats import chi2
    
    # Creating a copy of the input data
    rawdata = np.copy(data) 

    # Creating an initial mask of zeros to contain flagged values 
    Threshmask = np.zeros_like(data) 
    
    if axis == 0: # Running the algorithm along the frequency axis
        flipped = np.transpose(np.copy(rawdata))# Flipping the input data and mask for easier looping
        flippedmask = np.transpose(np.copy(Threshmask))
        
        dims = np.shape(flipped) #Dimensions of input array
        
        # Calulating the number of iterations based on the length of the input data (in exponential factors of 2)
        iterations = 7 #If full number of iterations is desired, use: np.int(np.log(dims[1])/np.log(2))
        
        thresholds = [] # Storing the base threshold for each column, to calculate subsequent thresholds
        for i in range(iterations+1): # Loop to iterate several times and to increase sub-sequence lengths
            sublen = 2**i
            print 'Testing Sub-sequence length: ', sublen
            dof = 2*128*sublen # Degrees of freedom, 256 per pixel (128 pixels summed into one, each with Re and Im components)
            basechi = chi2(dof).rvs(size=1E5) # Calculating a base chi-squared distribution (not scaled or shifted)
            med_chi = np.median(basechi) # Median of chi-sq distribution
            mad_chi = np.median(np.abs(basechi - med_chi)) # MAD of chi-sq distribution
            sumdata = np.reshape(flipped, (dims[0],dims[1]/sublen,sublen)).sum(axis=2)
            for j in range(dims[0]): # Iterating along the columns and flagging in each one
                col = flipped[j] 
                sumcol = sumdata[j]
                med = np.median(sumcol) # Median of channel
                mad = np.median(np.abs(sumcol-med)) # Median absolute deviation (MAD) of channel
                offset = med - med_chi * (mad/mad_chi)
                chidist = chi2(dof,scale = mad/mad_chi, loc=offset)
                thresh = chidist.ppf(1-(1E-15))
                 
                for k in np.arange(0,len(col),sublen): # Loop to slide window along freq channels
                    window = np.copy(col[k:k+sublen])
                    for l in range(len(window)):  # Loop to set values which have been previously flagged as zero to later calculate the mean of the window
                        if flippedmask[j][k+l] == 1:
                            window[l] = 0
                    windowsum = np.sum(window)
                    nonzeros = len(np.nonzero(window)[0]) # Counting the number of non-zero values
                    windowmean = windowsum/nonzeros
                    for m in range(len(window)):
                        if flippedmask[j][k+m] == 1:
                            window[m] = windowmean
                    windowthresh = np.sum(window)
                    flipped[j][k:k+sublen] = window # Setting the values in the data set to the values of the window

                    if windowthresh > thresh: # If the average of the window is greater than the threshold, flag all values inside
                        flippedmask[j][k:k+sublen] = 1
        
        # Flipping the mask back to its original orientation
        reflippedmask = np.transpose(flippedmask)
        
        return reflippedmask
    
    elif axis == 1: # Running the algorithm along the time axis 
       
        #Dimensions of input array
        dims = np.shape(rawdata) 
        
        # Calulating the number of iterations based on the length of the input data (in exponential factors of 2)
        iterations = 7 #If full number of iterations is desired, use: np.int(np.log(dims[1])/np.log(2))
        
        for i in range(iterations+1): # Loop to iterate several times and to increase sub-sequence lengths ('sublen' variable)
            sublen = 2**i
            print 'Testing Sub-sequence length: ', sublen
            dof = 2*128*sublen # Degrees of freedom, 256 per pixel (128 pixels summed into one, each with Re and Im components)
            basechi = chi2(dof).rvs(size=1E5) # Calculating a base chi-squared distribution (not scaled or shifted)
            med_chi = np.median(basechi) # Median of chi-sq distribution
            mad_chi = np.median(np.abs(basechi - med_chi)) # MAD of chi-sq distribution
            sumdata = np.reshape(rawdata, (dims[0],dims[1]/sublen,sublen)).sum(axis=2)
            for j in range(np.shape(rawdata)[0]):
                channel = rawdata[j]
                sumchan = sumdata[j]
                med = np.median(sumchan) # Median of channel
                mad = np.median(np.abs(sumchan-med)) # Median absolute deviation (MAD) of channel 
                offset = med - med_chi * (mad/mad_chi)
                chidist = chi2(dof,scale = mad/mad_chi, loc=offset)
                thresh = 1.2*chidist.ppf(1-(1E-15))
                
                for k in np.arange(0,len(channel),sublen): # Loop to slide window along freq channels
                    window = np.copy(channel[k:k+sublen])
                    for l in range(len(window)): # Loop to set values which have been previously flagged as zero to later calculate the mean of the window
                        if Threshmask[j][k+l] == 1:
                            window[l] = 0
                    windowsum = np.sum(window)
                    nonzeros = len(np.nonzero(window)[0]) # Counting number of non-zeros
                    windowmean = windowsum/nonzeros
                    for m in range(len(window)):
                        if Threshmask[j][k+m] == 1:
                            window[m] = windowmean
                    windowthresh = np.sum(window)
                    rawdata[j][k:k+sublen] = window # Setting the values in the data set to the values of the window
                                                      
                    if windowthresh > thresh: # If the average of the window is greater than the threshold, flag all values inside
                        Threshmask[j][k:k+sublen] = 1 
        
        return Threshmask
    
    
    elif axis == 'both': # Runs the algorithm along both time and frequency axes, then combines the output masks
        
        # Flipping the input data and mask for easier looping in frequency-axis
        flipped = np.transpose(np.copy(rawdata))
        flippedmask = np.transpose(np.zeros_like(Threshmask))
        
        ### First operation, going along time axis ### 
        
        #Dimensions of input array
        dims1 = np.shape(rawdata) 
        
        # Calulating the number of iterations based on the length of the input data (in exponential factors of 2)
        iterations = 6 #If full number of iterations is desired, use: np.int(np.log(dims[1])/np.log(2))
        
        for i in range(iterations+1): # Loop to iterate several times and to increase sub-sequence lengths ('sublen' variable)
            sublen = 2**i
            print 'Testing Sub-sequence length: ', sublen
            dof = 2*128*sublen # Degrees of freedom, 256 per pixel (128 pixels summed into one, each with Re and Im components)
            basechi = chi2(dof).rvs(size=1E5) # Calculating a base chi-squared distribution (not scaled or shifted)
            med_chi = np.median(basechi) # Median of chi-sq distribution
            mad_chi = np.median(np.abs(basechi - med_chi)) # MAD of chi-sq distribution
            sumdata1 = np.reshape(rawdata, (dims1[0],dims1[1]/sublen,sublen)).sum(axis=2)
            for j in range(np.shape(rawdata)[0]):
                channel = rawdata[j]
                sumchan = sumdata1[j]
                med = np.median(sumchan) # Median of channel
                mad = np.median(np.abs(sumchan-med)) # Median absolute deviation (MAD) of channel 
                offset = med - med_chi * (mad/mad_chi)
                chidist = chi2(dof,scale = mad/mad_chi, loc=offset)
                thresh = 1.2*chidist.ppf(1-(1E-15))
                
                for k in np.arange(0,len(channel),sublen): # Loop to slide window along freq channels
                    window = np.copy(channel[k:k+sublen])
                    for l in range(len(window)): # Loop to set values which have been previously flagged as zero to later calculate the mean of the window
                        if Threshmask[j][k+l] == 1:
                            window[l] = 0
                    windowsum = np.sum(window)
                    nonzeros = len(np.nonzero(window)[0])
                    windowmean = windowsum/nonzeros
                    for m in range(len(window)):
                        if Threshmask[j][k+m] == 1:
                            window[m] = windowmean
                    windowthresh = np.sum(window)
                    rawdata[j][k:k+sublen] = window # Setting the values in the data set to the values of the window
                                                      
                    if windowthresh > thresh: # If the sum of the window is greater than the threshold, flag all values inside
                        Threshmask[j][k:k+sublen] = 1
                    
        
        ### Second operation, going along the frequency axis ###
        
        dims2 = np.shape(flipped) #Dimensions of input array
        
        # Calulating the number of iterations based on the length of the input data (in exponential factors of 2)
        iterations = 6 #If full number of iterations is desired, use: np.int(np.log(dims[1])/np.log(2))
        
        for i in range(iterations+1): # Loop to iterate several times and to increase sub-sequence lengths ('sublen' variable)
            sublen = 2**i
            print 'Testing Sub-sequence length: ', sublen
            dof = 2*128*sublen # Degrees of freedom, 256 per pixel (128 pixels summed into one, each with Re and Im components)
            basechi = chi2(dof).rvs(size=1E5) # Calculating a base chi-squared distribution (not scaled or shifted)
            med_chi = np.median(basechi) # Median of chi-sq distribution
            mad_chi = np.median(np.abs(basechi - med_chi)) # MAD of chi-sq distribution
            sumdata2 = np.reshape(flipped, (dims2[0],dims2[1]/sublen,sublen)).sum(axis=2)
            for j in range(dims2[0]):
                col = flipped[j] 
                sumcol = sumdata2[j]
                med = np.median(sumcol) # Median of channel
                mad = np.median(np.abs(sumcol-med)) # Median absolute deviation (MAD) of channel
                offset = med - med_chi * (mad/mad_chi)
                chidist = chi2(dof,scale = mad/mad_chi, loc=offset) #Recalculating a re-scaled and shifted chisq distribution
                thresh = chidist.ppf(1-(1E-15))
                 
                for k in np.arange(0,len(col),sublen): # Loop to slide window along freq channels
                    window = np.copy(col[k:k+sublen])
                    for l in range(len(window)): # Loop to set values which have been previously flagged as zero to later calculate the mean of the window
                        if flippedmask[j][k+l] == 1:
                            window[l] = 0
                    windowsum = np.sum(window)
                    nonzeros = len(np.nonzero(window)[0])
                    windowmean = windowsum/nonzeros
                    for m in range(len(window)):
                        if flippedmask[j][k+m] == 1:
                            window[m] = windowmean
                    windowthresh = np.sum(window)
                    flipped[j][k:k+sublen] = window # Setting the values in the data set to the values of the window

                    if windowthresh > thresh: # If the average of the window is greater than the threshold, flag all values inside
                        flippedmask[j][k:k+sublen] = 1

        # Flipping the mask back to its original orientation
        reflippedmask = np.transpose(flippedmask)
        
        combinedmask = np.zeros_like(data) # Mask to combine both time-axis and freq-axis masks
        print 'Combining the two masks...'
        
        # Loop to iterate through both masks and create a combined mask
        for n in range(np.shape(reflippedmask)[0]): 
            for m in range(np.shape(reflippedmask)[1]):
                if reflippedmask[n][m] == 1 or Threshmask[n][m] == 1:
                    combinedmask[n][m] = 1
        
        print 'Finished SumThreshold.' # For testing purposes
        return combinedmask