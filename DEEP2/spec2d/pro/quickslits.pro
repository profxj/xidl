pro quickslits, filename, slitnumbers
;+
; NAME:
;   QUICKSLITS
;
;
; PURPOSE:
;   Parse a file from the quicklook widget into a list of slits to reduce
;
;
; CATEGORY:
;   Quicklook
;
;
; CALLING SEQUENCE:
;   quickslits,filename,slitnumbers
;
; INPUTS:
;   filename - the name/path of a file output by the quicklook widget
;
; OUTPUTS:
;   slitnumbers - an array of the numbers of the slitlets to reduce
; MODIFICATION HISTORY:
;   JAN 8/7/02
;-

  if file_test(filename) eq 0 then begin 
     message, ' Cannot find the file - do you have the path right?'
  endif


  read_planfile, filename, maskno, rawdatadir, outdatadir, $
    flatnames, arcnames, sciencenames, listofslitnums, $
    nbright, ndim, nsimp, nlong

  bfile = '*bintab*.fits'
  if file_test(bfile) ne 0 then bfile = bfile[0] $
    else bfile = arcnames[0]

  header = headfits(bfile)
  
   deimos_grating,header,g_rule,grangle,lambda_c

   slider=sxpar(header,'GRATEPOS')

   if slider eq 3 then tltval=sxpar(header,'G3TLTVAL') $
      else tltval=sxpar(header,'G4TLTVAL') 

   grname = floor(g_rule+1E-5)

; ; determine parameters given grating, slider, and tilt
   omodel_params, grname, slider, tltval, roll, o3, mu


; ; set up system variables
   deimos_map=sysinit(mu,gmm=g_rule,o3=o3,roll3=roll)


simple_tables, bfile, slitnames=slitnames, mags=mags, slitwid=slitwid, $
  slitlen=slitlen, xmm=xmm, ymm=ymm, objnames=objnames
  

 npoints = 200
 nslits = n_elements(xmm)
 lambda = findgen(npoints)*30.  +4000.



; it looks like that error no longer is occurring.

; pick the appropriate amap & bmap for this data
choose_amapbmap,header,amapfile,bmapfile

; now call optical model for an array of points

;xi1/yi1 = mosaic coordinate system
;xp1/yp1 = CCD coordinate system

qmodel, deimos_map, replicate(1, npoints)#xmm, $
    replicate(1, npoints)#ymm, lambda#replicate(1, nslits),  $
     xi1, yi1,  ccdnum1, xp1, yp1, amapfile, $
      bmapfile, /cubic, goodloc=goodloc

xp1 = reform(xp1, npoints, nslits, /over)
yp1 = reform(yp1, npoints, nslits, /over)
ccdnum1 = reform(ccdnum1, npoints, nslits, /over)
goodloc = reform(goodloc, npoints, nslits, /over)

medx = fltarr(nslits)

for i=0, nslits-1 do begin
   whgood = where(goodloc[*, i], goodct)
   if goodct gt 0 then medx[i] = median(xp1[whgood, i])
endfor


  lengths = slitlen


  slitnums = slitnames

  nobj = n_elements(objectcat)
  nslit = n_elements(slitnames)
 
  isastar = slitwid gt 2*median(slitwid)


  isdeep = total(stregex(objnames,'^[0-9]{7}$') LT 0) EQ 0

  issky = slitwid*0

  if isdeep then begin

     issky = (long(objnames) MOD 1E6) ge 5E5
  endif

  
  slitnumbers = slitnames[where(isastar)]

; don't want first or last slit!
  sortbright = sort(mags+1E10*isastar+1E10*issky $
          +1E10*(findgen(nslit) eq 0) +1E10*(findgen(nslit) eq (nslit-1)))

  sortdim = sort(-mags+1E10*isastar+1E10*issky $
         +1E10*(findgen(nslit) eq 0) +1E10*(findgen(nslit) eq (nslit-1)))

  sortlong = sort(-lengths +1E10*isastar+1E10*issky +1E10*(lengths gt 12) $
           +1E10*(findgen(nslit) eq 0) +1E10*(findgen(nslit) eq (nslit-1)))

; a basic score definition, may want to refine this.  lower='simpler'
; to reduce.  this does seem to work OK, even though it looks
; absolutely horrible.

; penalties:
  score = 1E3*isastar+1E3*issky ; don't want star/sky slits
;  score = score+ 1E3*(nobjonslit ne 1) ; don't want multiobject slits
  score = score+ 1E3*(lengths lt 6) ; don't want short slits
  score = score+ 2*(lengths gt 11)  ; super-long slits aren't ideal
  score = score+abs(mags-23.25)*(1.+(mags gt 23.25))
 score = score +5*(long(slitnames) le 10)+5*(long(slitnames) ge (nslit-11)) ;don't want first/last slits
; don't want slits at edges of chips
 score = score+5*(medx lt 75)+5*(medx gt 1972)


  sortscore = slitnames[sort(score)]
  magscore = mags[sort(score)]

  if ndim gt 0 then slitnumbers = [slitnumbers, sortdim[0:ndim-1]]
  if nbright gt 0 then slitnumbers = [slitnumbers, sortbright[0:nbright-1]]
  if nlong gt 0 then slitnumbers = [slitnumbers, sortlong[0:nlong-1]]
  if nsimp gt 0 then slitnumbers = [slitnumbers, sortscore[0:nsimp-1]]

  slitnumbers = long(slitnumbers)

; get only the unique elements
  if total(listofslitnums) ge 1 then slitnumbers = [slitnumbers, listofslitnums]

  slitnumbers = slitnumbers[sort(slitnumbers)]
  slitnumbers = slitnumbers[uniq(slitnumbers)]
  
  return
end







