def dataIO(filenum,time=False, filepath='/mnt/scratch-lustre/simard/B0329_080416/20160805T202001Z_aro_raw/*vdif'):
    """Loads in data based on given file number and returns a numpy array of the raw voltages
    
    Keyword arguments:
    filenum -- int , the file-number of the data sets you want to load
    time -- boolean, if True, it will return the time of observation
    filepath -- String, full filepath to the data set (defaults to filepath for b0329+54 data set)
    
    """
    print 'Loading in Data Set...'
    
    # Importing required modules 
    import numpy as np
    from baseband import vdif
    import astropy.units as u
    from astropy.time import Time
    import glob
    
    # Defining variables (default values for ARO data)
    nfreqs = 1024 # Number of frequency channels 
    sample_rate = 800/(2.*nfreqs) * u.MHz
    size = 2**16
    
    # Loading in the data
    filelist = np.sort(glob.glob(filepath))
    fn = filelist[filenum]
    fh = vdif.open(fn, mode='rs', sample_rate=sample_rate)
    d = fh.read(size)
    
    #Rearranging the dimensions of the array
    d = np.einsum('ijk->ikj',d)
    
    #Loop to manually parse the header and extract the start time information
    if time == True:
        head = str(fh)
        timeindex = head.find('time=') + 5
        datestr = []
        for j in np.arange(timeindex,timeindex+29,1):
            datestr.append(head[j])
   
        timestr = ''.join(datestr)
        tstart = Time(timestr, precision=9, format = 'isot')
        print 'Finished loading in Data Set.'
        return d, timestr
    
    else:
        print 'Finished loading in Data Set.'
        return d