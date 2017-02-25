PRO nirspec_superflat, infiles, darkfiles $
                       ,superdarkfile, pixflatfile, illumflatfile $
                       , slitfile = slitfile $
                       , objfile = objfile $
                       , verbose = verbose, indir = indir $
                       , tempdir = tempdir $
                       , use_illum = use_illum $
                       , use_pixel = use_pixel $
                       , filter = filter, MXSHFT=mxshft $
                       , npoly = npoly, CHK = CHK $
                       , _EXTRA = extra 

tset_slits=mrdfits(slitfile,1)
;; Create the superdark
darkimg = niri_superdark(darkfiles,/NIRSPEC)
mwrfits,darkimg,superdarkfile,/create
;; Create the wavelength image
nirspec_proc,objfile,image,hdr=hdr
filetemp=fileandpath(objfile,path=path)
wavefile = 'wave-' + filetemp
waveimg=nirspec_waveimg(image,tset_slits,hdr,outfile=wavefile,filter=filter,mxshft=mxshft)

long_superflat,infiles,pixflatfile,illumflatfile $
               ,slitfile=slitfile,wavefile=wavefile $
               ,use_illum=use_illum,use_pixel=use_pixel $
               ,tempdir=tempdir,slitsamp = 5.0 $
               ,biasfile=superdarkfile,/NIRSPEC
spawn, 'rm -f ' + superdarkfile
spawn, 'rm -f ' + wavefile

END
