
path = '/Users/joe/GMOS_redux/'
cd, path

rawpath = path + 'Raw/'
irafpath = path + 'IRAF/'

biasfiles = ['S20130904S0108.fits' $
             , 'S20130904S0109.fits' $
             , 'S20130904S0110.fits' $
             , 'S20130904S0111.fits' $
             , 'S20130904S0112.fits']

darkfiles =  ['S20130905S0038.fits' $
              , 'S20130905S0039.fits' $
              , 'S20130905S0040.fits' $
              , 'S20130905S0070.fits' $
              , 'S20130905S0071.fits' $
              , 'S20130905S0072.fits' $
              , 'S20130905S0073.fits' $
              , 'S20130905S0074.fits'] 

flatfiles = ['S20130903S0032.fits' $
             , 'S20130903S0033.fits' $
             , 'S20130903S0034.fits' $
             , 'S20130903S0035.fits' $
             , 'S20130903S0036.fits' $
             , 'S20130903S0037.fits' $
             , 'S20130903S0038.fits']


datafiles = ['S20130905S0044.fits' $
             , 'S20130905S0046.fits' $
             , 'S20130905S0048.fits' $
             , 'S20130905S0051.fits' $
             , 'S20130905S0053.fits' $
             , 'S20130905S0055.fits']

gmos_imag_reduce, biasfiles, darkfiles, flatfiles, datafiles $
                  , rawpath = rawpath, irafpath = irafpath

END






