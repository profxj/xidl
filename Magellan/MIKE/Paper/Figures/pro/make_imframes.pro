
modelfile = 'Extract/Obj_mr0006-Model.fits.gz'
imagefile = ''Final/f_mr0006.fits.gz'
im = mrdfits(imagefile)
objmodel = mrdfits(modelfile,3)
skymodel = mrdfits(modelfile,4)



x1 = 100
x2 = 259
y2 = 940
y1 = 341

; first one just to show the SLLS in PKS2000
outfile = 'data/PKS2000_im1.fits'
mwrfits, im[x1:x2, y1:y2], outfile, /create
mwrfits, objmodel[x1:x2, y1:y2], outfile
mwrfits, skymodel[x1:x2, y1:y2], outfile
residual = im[x1:x2, y1:y2] - objmodel[x1:x2, y1:y2] - skymodel[x1:x2, y1:y2]
mwrfits, residual, outfile


; second one just to show the sky subtraction in PKS2000

x1 = 814
x2 = x1 + 159
y1 = 1000
y2 = y1 + 599

outfile = 'data/PKS2000_im2.fits'
mwrfits, im[x1:x2, y1:y2], outfile, /create
mwrfits, objmodel[x1:x2, y1:y2], outfile
mwrfits, skymodel[x1:x2, y1:y2], outfile
residual = im[x1:x2, y1:y2] - objmodel[x1:x2, y1:y2] - skymodel[x1:x2, y1:y2]
mwrfits, residual, outfile

end


