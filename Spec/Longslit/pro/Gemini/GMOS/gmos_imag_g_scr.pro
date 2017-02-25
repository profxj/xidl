
path = '/Users/joe/GMOS_redux/'
cd, path

rawpath = path + 'Raw/'
irafpath = path + 'IRAF/'

biasfiles = ['S20130904S0108.fits' $
             , 'S20130904S0109.fits' $
             , 'S20130904S0110.fits' $
             , 'S20130904S0111.fits' $
             , 'S20130904S0112.fits']

darkfiles =  ['S20130906S0066.fits' $
              , 'S20130906S0067.fits' $
              , 'S20130906S0068.fits' $
              , 'S20130906S0069.fits' $
              , 'S20130906S0070.fits']

;; These are archival flats from 11-2012. Replace with updated ones
flatfiles = ['S20121120S0053.fits' $
             , 'S20121120S0054.fits' $
             , 'S20121120S0055.fits' $
             , 'S20121120S0056.fits' $
             , 'S20121120S0057.fits']    


datafiles = ['S20130905S0043.fits' $
             , 'S20130905S0045.fits' $
             , 'S20130905S0047.fits' $
             , 'S20130905S0049.fits' $
             , 'S20130905S0050.fits' $
             , 'S20130905S0052.fits' $
             , 'S20130905S0054.fits' $
             , 'S20130905S0056.fits']

gmos_imag_reduce, biasfiles, darkfiles, flatfiles, datafiles $
                  , rawpath = rawpath, irafpath = irafpath

END






