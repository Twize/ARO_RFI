# List of hard-coded frequency channels with narrow-band RFI at ARO for reference
channelarray = [384,541,757,808,809,810,811,816,817,818,867,868,870,889,893,896,909,962]

def upchannelize(Vin, channels, upfactor=128, refchannel=865): 
    '''Function which uses fast fourier transforms to up-channelize frequency channels, remove RFI and then re-transforms back
    
    The function iterates through given frequency channels, up-channelizes, then backfills RFI-affected portions with
    values pulled from a gaussian distribution based on a 'clean' portion of the data, before inverse-FFT'ing the data back 
    and replacing corrupted channels in the input data set with the cleaned channels. Persistent narrow-band RFI specific to ARO
    is harcoded, but if a channel is given that is not hardcoded, it will prompt for user to select 'clean' portion, then the
    portions which should be cut/backfilled. 
    
    Keyword arguments:
    Vin -- Input data, must be voltages (not power) 
    channels -- List of channels which contain RFI and should be up-channelized
    upfactor -- Factor by which to up-channelize by (number of frequencies which each channel will be split into)(default = 128)
    refchannel -- A 'Clean' channel from the original data to use as a reference when re-normalizing outputted cleaned channels (default = 865)
    '''
    print 'Beginning Upchannelizing...'
    
    #Importing required modules
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Reference channel, for re-scaling the cleaned output channels
    reference = Vin[:,refchannel,0]
    refpow = np.real(reference)**2. + np.imag(reference)**2.
    refmean = np.mean(refpow)
    
    # Creating copy of input voltage array, to replace corrupted channels
    cleanedoutput = np.copy(Vin)
    
    for channel in channels:
        print channel
        # Selecting channel to up-channelize
        lineRFI = Vin[:,channel,0] 
        
        # Performing FFT on this frequency channel
        FFT = np.fft.fft(lineRFI.reshape([-1,upfactor]),axis=1)
        
        # Creating copy of FFT output array, to backfill with gaussian data
        cleaned = np.copy(FFT)
        
        # Hard-coding channels with persistent narrow-band RFI
        if channel == 384 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,0:20] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,0:20] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 541 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,45:100])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,45:100]))
            Remean = np.mean(np.real(FFT[:,45:100])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,45:100]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,30:45] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,30:45] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 757 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,20:35] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,20:35] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 808 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,20:60])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,20:60]))
            Remean = np.mean(np.real(FFT[:,20:60])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,20:60]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,75:100] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,75:100] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 809 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,0:40])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,0:40]))
            Remean = np.mean(np.real(FFT[:,0:40])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,0:40]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,40:100] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,40:100] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 810 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,70:110])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,70:110]))
            Remean = np.mean(np.real(FFT[:,70:110])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,70:110]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,40:60] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,40:60] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 811 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,115:128] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,115:128] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 816 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,0:40] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,0:40] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 817 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,0:40] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,0:40] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 818 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,0:40] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,0:40] = ReGauss1 + ImGauss1*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 867 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,0:40])) # STD of both real and imaginary components for back-filling
            Imstd = np.std(np.imag(FFT[:,0:40]))
            Remean = np.mean(np.real(FFT[:,0:40])) # Mean of both real and imaginary components for back-filling
            Immean = np.mean(np.imag(FFT[:,0:40]))
            
            # Hard-coded the backfilling of bad frequencies for corrupted channel 867
            RFI1 = cleaned[:,50:80] 
            dims1 = np.shape(RFI1)
            # Generating a gaussian distribution based on Re and Im components, to back-fill with
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) 
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,50:80] = ReGauss1 + ImGauss1*1j
            
            RFI2 = cleaned[:,78:85] 
            dims2 = np.shape(RFI2)
            ReGauss2 = np.random.normal(loc=Remean,scale=Restd,size=(dims2[0],dims2[1]))
            ImGauss2 = np.random.normal(loc=Immean,scale=Imstd,size=(dims2[0],dims2[1]))
            cleaned[:,78:85] = ReGauss2 + ImGauss2*1j
            
            RFI3 = cleaned[:,87:94] 
            dims3 = np.shape(RFI3)
            ReGauss3 = np.random.normal(loc=Remean,scale=Restd,size=(dims3[0],dims3[1]))
            ImGauss3 = np.random.normal(loc=Immean,scale=Imstd,size=(dims3[0],dims3[1]))
            cleaned[:,87:94] = ReGauss3 + ImGauss3*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2. # Calculating power of inv, for re-normalizing
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
        
        elif channel == 868 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,15:40])) 
            Imstd = np.std(np.imag(FFT[:,15:40]))
            Remean = np.mean(np.real(FFT[:,15:40])) 
            Immean = np.mean(np.imag(FFT[:,15:40]))
            
            # Manually backfilling bad frequencies
            RFI1 = cleaned[:,45:90] 
            dims1 = np.shape(RFI1) # Extracting dimensions of RFI infected chunks
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) # Pulling a Re and Im component from gaussian distribution to backfill
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,45:90] = ReGauss1 + ImGauss1*1j # Adding Re and Im components together
            
            RFI2 = cleaned[:,0:5]
            dims2 = np.shape(RFI2)
            ReGauss2 = np.random.normal(loc=Remean,scale=Restd,size=(dims2[0],dims2[1]))
            ImGauss2 = np.random.normal(loc=Immean,scale=Imstd,size=(dims2[0],dims2[1]))
            cleaned[:,0:5] = ReGauss2 + ImGauss2*1j
            
            RFI3 = cleaned[:,120:128] 
            dims3 = np.shape(RFI3)
            ReGauss3 = np.random.normal(loc=Remean,scale=Restd,size=(dims3[0],dims3[1]))
            ImGauss3 = np.random.normal(loc=Immean,scale=Imstd,size=(dims3[0],dims3[1]))
            cleaned[:,120:128] = ReGauss3 + ImGauss3*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
        elif channel == 870 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,20:60])) 
            Imstd = np.std(np.imag(FFT[:,20:60]))
            Remean = np.mean(np.real(FFT[:,20:60])) 
            Immean = np.mean(np.imag(FFT[:,20:60]))
            
            # Manually backfilling bad frequencies
            RFI1 = cleaned[:,70:80] 
            dims1 = np.shape(RFI1) # Extracting dimensions of RFI infected chunks
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) # Pulling a Re and Im component from gaussian distribution to backfill
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,70:80] = ReGauss1 + ImGauss1*1j # Adding Re and Im components together
            
            RFI2 = cleaned[:,82:95] 
            dims2 = np.shape(RFI2)
            ReGauss2 = np.random.normal(loc=Remean,scale=Restd,size=(dims2[0],dims2[1]))
            ImGauss2 = np.random.normal(loc=Immean,scale=Imstd,size=(dims2[0],dims2[1]))
            cleaned[:,82:95] = ReGauss2 + ImGauss2*1j

            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
        elif channel == 889 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) 
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) 
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Manually backfilling bad frequencies
            RFI1 = cleaned[:,15:30] 
            dims1 = np.shape(RFI1) # Extracting dimensions of RFI infected chunks
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) # Pulling a Re and Im component from gaussian distribution to backfill
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,15:30] = ReGauss1 + ImGauss1*1j # Adding Re and Im components together

            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
        elif channel == 893 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:60])) 
            Imstd = np.std(np.imag(FFT[:,40:60]))
            Remean = np.mean(np.real(FFT[:,40:60])) 
            Immean = np.mean(np.imag(FFT[:,40:60]))
            
            # Manually backfilling bad frequencies
            RFI1 = cleaned[:,70:90] 
            dims1 = np.shape(RFI1) # Extracting dimensions of RFI infected chunks
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) # Pulling a Re and Im component from gaussian distribution to backfill
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,70:90] = ReGauss1 + ImGauss1*1j # Adding Re and Im components together

            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
        elif channel == 896 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,40:80])) 
            Imstd = np.std(np.imag(FFT[:,40:80]))
            Remean = np.mean(np.real(FFT[:,40:80])) 
            Immean = np.mean(np.imag(FFT[:,40:80]))
            
            # Manually backfilling bad frequencies
            RFI1 = cleaned[:,0:10] 
            dims1 = np.shape(RFI1) # Extracting dimensions of RFI infected chunks
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) # Pulling a Re and Im component from gaussian distribution to backfill
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,0:10] = ReGauss1 + ImGauss1*1j # Adding Re and Im components together

            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
        elif channel == 909 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,75:100])) 
            Imstd = np.std(np.imag(FFT[:,75:100]))
            Remean = np.mean(np.real(FFT[:,75:100])) 
            Immean = np.mean(np.imag(FFT[:,75:100]))
            
            # Manually backfilling bad frequencies
            RFI1 = cleaned[:,25:45] 
            dims1 = np.shape(RFI1) # Extracting dimensions of RFI infected chunks
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) # Pulling a Re and Im component from gaussian distribution to backfill
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,25:45] = ReGauss1 + ImGauss1*1j # Adding Re and Im components together

            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
        elif channel == 962 and upfactor == 128:
            # Mean and Std of clean portions of the Up-channelized Sections for Backfilling
            Restd = np.std(np.real(FFT[:,0:40])) 
            Imstd = np.std(np.imag(FFT[:,0:40]))
            Remean = np.mean(np.real(FFT[:,0:40])) 
            Immean = np.mean(np.imag(FFT[:,0:40]))
            
            # Manually backfilling bad frequencies
            RFI1 = cleaned[:,41:60] 
            dims1 = np.shape(RFI1) # Extracting dimensions of RFI infected chunks
            ReGauss1 = np.random.normal(loc=Remean,scale=Restd,size=(dims1[0],dims1[1])) # Pulling a Re and Im component from gaussian distribution to backfill
            ImGauss1 = np.random.normal(loc=Immean,scale=Imstd,size=(dims1[0],dims1[1]))
            cleaned[:,41:60] = ReGauss1 + ImGauss1*1j # Adding Re and Im components together

            RFI2 = cleaned[:,85:128] 
            dims2 = np.shape(RFI2)
            ReGauss2 = np.random.normal(loc=Remean,scale=Restd,size=(dims2[0],dims2[1]))
            ImGauss2 = np.random.normal(loc=Immean,scale=Imstd,size=(dims2[0],dims2[1]))
            cleaned[:,85:128] = ReGauss2 + ImGauss2*1j
            
            # Inverse FFT'ing the corrupted freq channels
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
        else: # If given channel is not hard-coded, ask user for input on clean chunk & portions to cut
            ddd = np.real(FFT)**2 + np.imag(FFT)**2 # Power, for displaying plot before asking for input
            
            # Generating a plot of the power, so user can determine clean chunks and portions to cut
            fig = plt.figure(figsize=[10,10])
            plt.imshow(np.log10(ddd),aspect='auto', interpolation='none', vmin=-3,vmax=6)
            plt.title('Up-Channelized Frequency Channel %s' %channel)
            cbar = plt.colorbar()
            cbar.set_label('Log_10(Power)')
            plt.xlabel('Frequency Channel')
            plt.ylabel('Time Sample')
            plt.show()
            
            while True:
                # Asking user for input on clean portion to use as reference for back-filling
                cleaninput = np.array(list(raw_input('Where is the clean portion for reference (use slicing notation x:y)? ')))
                if np.all(':' in cleaninput) == False: # Raise error if user doesn't use list-slicing notation
                    print 'Improper input, please use list slicing notation --> x:y'
                    continue  # Continue to ask user for proper input until it is received
                else: # When proper input is provided, continue
                    break
            
            mid = np.where(cleaninput == ':') # Finding index where character ':' is located. 
                                              #... this denotes the split between start and end of interval
            
            cleanstart = cleaninput[:mid[0]]
            cleanstart = ''.join(cleanstart)
            cleanstart = np.int(cleanstart) # Converting user input to int's for list-slicing

            cleanend = cleaninput[mid[0]+1:]
            cleanend = ''.join(cleanend)
            cleanend = np.int(cleanend)
            
            # Conditions to ensure input was of correct format for list slicing
            if cleanstart > cleanend:
                print 'Start of interval is greater than the end'
                continue
            elif cleanstart < 0 or cleanend < 0:
                print 'Input value(s) were negative'
                continue
            elif cleanstart == cleanend:
                print 'Start of interval was equal to end of interval'
                continue

            # Finding the STD and mean of the clean portion provided as input
            Restd = np.std(np.real(FFT[:,cleanstart:cleanend])) 
            Imstd = np.std(np.imag(FFT[:,cleanstart:cleanend]))
            Remean = np.mean(np.real(FFT[:,cleanstart:cleanend])) 
            Immean = np.mean(np.imag(FFT[:,cleanstart:cleanend]))
            
            # Asking user for input on portions with RFI to cut and back-fill
            RFIchunks = []
            while True:
                RFIinput = np.array(list(raw_input('Which sections would you like to cut and backfill? (use slicing notation x:y)? ')))
                if np.all(':' in RFIinput) == False:
                    print 'Improper input, please use list slicing notation --> x:y'
                    continue
                else:
                    RFIchunks.append(RFIinput)
                    choice = raw_input('Continue? y/n: ') # Give user option to enter more chunks to cut/backfill
                    if choice == 'y':
                        continue
                    elif choice == 'n':    
                        break
                    else:
                        print'Input not understood, please enter "y" or "n"'
                        continue
                        
            for j in range(len(RFIchunks)):
                mid = np.where(RFIchunks[j] == ':')
                
                RFIstart = RFIchunks[j][:mid[0]]
                RFIstart = ''.join(RFIstart)
                RFIstart = np.int(RFIstart)

                RFIend = RFIchunks[j][mid[0]+1:]
                RFIend = ''.join(RFIend)
                RFIend = np.int(RFIend)

                # Conditions to ensure input was of correct format for list slicing
                if RFIstart > RFIend:
                    print 'Start of interval is greater than the end'
                    continue
                elif RFIstart < 0 or RFIend < 0:
                    print 'Input value(s) were negative'
                    continue
                elif RFIstart == RFIend:
                    print 'Interval of improper size: size 0'
                    continue
                
                RFI = cleaned[:,RFIstart:RFIend]
                RFIdims = np.shape(RFI)
                ReGauss = np.random.normal(loc=Remean,scale=Restd,size=(RFIdims[0],RFIdims[1])) # Pulling a Re and Im component from gaussian distribution to backfill
                ImGauss = np.random.normal(loc=Immean,scale=Imstd,size=(RFIdims[0],RFIdims[1]))
                cleaned[:,RFIstart:RFIend] = ReGauss + ImGauss*1j # Adding Re and Im components together
            
            # Inverse FFT'ing the corrupted freq channel
            inv = np.fft.ifft(np.squeeze(cleaned),axis=1).reshape([65536,-1])
            pow = np.real(inv)**2. + np.imag(inv)**2.
            invmean = np.mean(pow)
            
            # Replacing old channels with new cleaned ones
            cleanedoutput[:,channel,0] = np.squeeze(np.sqrt(refmean/invmean)*inv)
            
    print 'Finished Upchannelizing.' # For testing purposes
    return cleanedoutput