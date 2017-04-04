def radiometer(data, chunklen, data2=None, label1=None, label2=None, savefig=False, decfactor=128, bandwidth=0.390625 * (10**6)): 
    '''Function which uses the radiometer equation as a quality metric for comparing RFI excision methods.
    
    The function iterates along frequency channels in chunks, and uses these chunks to compute the value of the 
    radiometer equation. The function produces histograms of the computed values of the radiometer equation and
    counts the number of outliers. 
    
    Keyword arguments:
    data -- Input data set 
    chunklen -- The length of the chunks which the algorithm will use
    data2 -- (OPTIONAL) Second data set, for comparison purposes
    label1 -- (OPTIONAL) Plot labels for first input data set (default = 'Data 1')
    label2 -- (OPTIONAL) Plot labels for second input data set (default = 'Data 2')
    savefig -- (OPTIONAL) Set to 'True' if you want to save the resulting histogram, will prompt for a filename
               and then saves to the working directory (default = 'False')
    decfactor -- (OPTIONAL) Factor by which the data sets have been decimated by (should be same for both sets) (default = 128)
    bandwidth -- (OPTIONAL) Bandwidth of each frequency channel (default for ARO = 0.390625 * (10**6) Hz)
            
    '''
    
    # Importing modules
    import numpy as np
    import matplotlib.pyplot as plt
    from os import getcwd
    
    
    # Defining variables 
    tsamp = 2.560000E-6 # Time per sample 
    t = tsamp * decfactor # Time per pixel
    
    # Conditions checking that chosen chunk length is proper given input data
    if chunklen > np.shape(data)[1] or chunklen == 0:
        raise ValueError('Improper chunk length')
    
    elif (np.shape(data)[1])%chunklen != 0:
        raise ValueError('Chunk length is not a proper fraction of data length')
    
    else:
        # Loop to calculate mean and RMS in chunks in each freq channel for data set 1, then comparing to radiometer equation
        outliercount = 0 # Counter to keep track of outliers
        data1rad = []
        for i in range(np.shape(data)[0]):
            channel = data[i]
            for j in np.arange(0,len(channel),chunklen):
                chunk = channel[j:j+chunklen]
                chunkmean = np.mean(chunk) # Calculating mean and RMS in each time chunk
                chunkRMS = np.std(chunk)
                rad = chunkmean/(chunkRMS*np.sqrt(bandwidth*t)) # radiometer equation
                if np.isnan(rad) == True or np.isinf(rad) == True:
                    data1rad.append(0) # If value of rad is inf or nan, set to 0, then exclude from histogram
                elif rad > 1.125 or rad < 0.825: # Condition to count number of outliers
                    outliercount += 1
                    data1rad.append(rad)
                else:
                    data1rad.append(rad)

        # Converting to array for histogramming
        data1rad = np.array(data1rad)

        # Calculating the percentage of outliers
        data1outlierpercent = (np.float(outliercount)/np.float(len(data1rad)))*100

        # Case for if two data sets are provided
        if data2 != None:
            # Same loop as above to calculate mean and RMS in chunks in each freq channel for second (optional) data set
            outliercount2 = 0 # Counter to keep track of outliers
            data2rad = []
            for i in range(np.shape(data2)[0]):
                channel = data2[i]
                for j in np.arange(0,len(channel),chunklen):
                    chunk = channel[j:j+chunklen]
                    chunkmean = np.mean(chunk[chunk!=np.float('NaN')]) # Calculating mean and RMS in each time chunk
                    chunkRMS = np.std(chunk[chunk!=np.float('NaN')])
                    rad = chunkmean/(chunkRMS*np.sqrt(bandwidth*t)) # radiometer equation
                    if np.isnan(rad) == True or np.isinf(rad) == True:
                        data2rad.append(0) # If value of rad is inf or nan, set to 0, then exclude from histogram
                    elif rad > 1.175 or rad < 0.825: # Condition to count number of outliers
                        outliercount2 += 1
                        data2rad.append(rad)
                    else:
                        data2rad.append(rad)

            # Converting to array for histogramming
            data2rad = np.array(data2rad)

            # Calculating the percentage of outliers
            data2outlierpercent = (np.float(outliercount2)/np.float(len(data2rad)))*100

            # Histogram of Radiometer Equation for data 1 and data 2
            plt.figure(figsize=[10,10])
            plt.title('Radiometer Quality Metric Using Time-Chunks of Length %s'%(chunklen))
            
            if label1 == None: # Setting legend labels based on input, otherwise set to default label
                label1 = 'Data 1'
                
            if label2 == None:
                label2 = 'Data 2'
                
            n, bins, patches = plt.hist(data1rad.flatten()[data1rad!=0],range=(0,4), bins=200, alpha=0.5, log=True, label=label1)
            n, bins, patches = plt.hist(data2rad.flatten()[data2rad!=0], range=(0,4), bins=200, alpha=0.5, log=True, label=label2)
            plt.xlabel('Value of Radiometer Metric')
            plt.ylabel('#')
            plt.legend(loc='upper right')
            
            
            if savefig == True:
                fname = raw_input('Enter a file name...')
                fpath = getcwd() # Extracting working directory using os.getcwd and using this as save location
                plt.savefig('%s/%s'%(fpath,fname))
                print 'File saved to %s/%s.png'%(fpath,fname)
            
            plt.show()

            print 'Outlier percentage of data set 1 using time chunks of length ', chunklen,':', data1outlierpercent
            print 'Outlier percentage of data set 2 using time chunks of length', chunklen,':', data2outlierpercent

            return (data1outlierpercent - data2outlierpercent)
        
        # Case for if only one data set is provided
        else: 
            # Histogram of Radiometer Equation for data 1 and data 2
            plt.figure(figsize=[10,10])
            plt.title('Radiometer Quality Metric Using Time-Chunks of Length %s, Outliers = %s%%'%(chunklen,data1outlierpercent))
            
            if label1 == None: # Setting legend labels based on input, otherwise set to default label
                label1 = 'Data 1'
            
            n, bins, patches = plt.hist(data1rad.flatten()[data1rad!=0],range=(0,4), bins=200, alpha=0.5, log=True, label=label1)
            plt.xlabel('Value of Radiometer Metric')
            plt.ylabel('#')
            plt.legend(loc='upper right')
            
            if savefig == True:
                fname = raw_input('Enter a file name...')
                fpath = getcwd()
                plt.savefig('%s/%s'%(fpath,fname))
                print 'File saved to %s/%s.png'%(fpath,fname)
                
            plt.show()

            print 'Outlier percentage of data set 1 using time chunks of length ', chunklen,':', data1outlierpercent
    
            return