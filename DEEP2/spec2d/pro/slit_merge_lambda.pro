;+
; NAME:
;    slit_merge_lambda
;
; PURPOSE:
;    combines slit*.fits.gz files into one large fits file, 
;    after rectifying and aligning them by wavelength, for viewing
;    all of a DEIMOS mask at once. 
;    last 10 columns of returned array specify the slit number
;    top 10 rows of returned array specify the wavelength
;
; CALLING SEQUENCE:
;    image = slit_merge_lambda(filelist, hdr, [/rebin, $
;                              blue_end=blue_end, red_end=red_end])
;
; INPUTS:
;    filelist -- list of slit.*.fits files
;
; OPTIONAL INPUTS:
;    /rebin -- If specified, the data will be smoothed by 3 pixels and
;              rebinned to 
;		half the size in each dimension.
;    blueend,redend -- The blue and red end cuts for wavelength range. 
;		This range should be wider than the wavelength range
;		of real data. 
;	  	If not specified, blueend=5500(in angstrom),redend=10000, 
;		for DEEP 1200 grating use only.
;
; KEYWORDS:
;   rebin -- set if you want to shrink array by 2x in each dimension
;
; OUTPUTS:
;   image-- a full frame of the DEIMOS data
;   hdr -- the header prepared for the FITS file, wavelength info
;            is included
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;   IDL>list=findfile('slit*')
;   IDl>image=slit_merge_lambda(list,hdr)
;   IDL>writefits,'Allslits.lambda.masknumber.grating.fits',image,hdr
; COMMENTS:
;   this is intended to be used to merge into one fits file all the 2d
;   sky subtracted slitlets for a given DEIMOS mask
;
; REVISION HISTORY:
;   Jun.20, 2002 by Renbin Yan 
;   Based on slit_merge.pro by Marc Davis
;   Modified on Jun.22,2002 by Renbin Yan
;   Revised 06oct2002 by mcc
;
;----------------------------------------------------------------------


function  slit_merge_lambda, filelist, hdr, rebin=rebin, $
                             blue_end=blue_end, red_end=red_end, $
                             nonlocal=nonlocal

; check that filelist was passed.
  if n_elements(filelist) eq 0 then begin
      print, '(slit_merge_lambda.pro) No file list passed!'
      filelist = findfile('slit*')
  endif
  nfiles = n_elements(filelist)
  if nfiles eq 0 then $
    message, 'No slit files found in working directory!'

; check if the red_end and blue_end optional parameters were
; passed. if not, use default values.
  if keyword_set(blue_end) then blue_end = blue_end[0] $
  else blue_end = 5500.
  if keyword_set(red_end) then red_end = red_end[0] $
  else red_end = 10000.

; define a small region in which to mark the slit number at the edge
; of the image.
  tail_l = 15.
; sort the file list! this routine assumes that the slitlist is
; sorted properly!!! 
  filelist = filelist[sort(filelist)]
; initialize the variable nrows. this array will be filled with the
; number of rows in each slit.
  nrows = intarr(n_elements(filelist))

; extract the design table from the bin table file.
  binfile = findfile('*bintab*', count=nbin)
  if nbin eq 0 then $
    message, 'No bin table file found in directory!' $
  else binfile = binfile[0]
  len = strlen(binfile)
  if strmid(binfile, len-3, 3) eq '.gz' then begin
      spawn, 'gunzip -v ' + binfile
      binfile = strmid(binfile, 0, len-3)
  endif
  fits_open, binfile, fcb
  dex = where(fcb.extname eq 'DesiSlits', tabcnt)
  if tabcnt eq 0 then begin
      fits_close, fcb
      message, 'No DesiSlits bin table found!'
  endif
  table = mrdfits(binfile, dex[0], /silent)
; establish a logical (called isanobj) which is 1 for program objects
; and 0 for alignment stars and sky slits.
  isanobj = table.slittyp eq 'P'
; within the design table, the objects are denoted by 'A' for
; alignment stars and 'P' for program/science objects. Note that
; sky-only slits are also calissified as 'P'. here we want to exclude
; the star slits since in the past their wavelength solution has not
; always been reliable. we don't want the goofy lambda range for those
; slits to throw off the routine. 
  objdex = where(table.slittyp eq 'P', objcnt)
  if objcnt eq 0 then $
    message, 'No objects found in design file!'
  objslits = string(table.slitname, format='(I3.3)')
; extract the slit number corresponding to each slit file from the
; name of the file.
  BRpos = strpos(filelist[0],'B.fits') > strpos(filelist[0],'R.fits')
  slitnames = strmid(filelist,BRpos-3,3)
  nslit = n_elements(objslits)
  for kk=0,nslit-1 do begin
      mdex = where(slitnames eq objslits[kk], mcnt)
      if mcnt eq 0 then isanobj[kk] = 0
  endfor
  whobj = where(isanobj eq 1, whcnt)
  if whcnt eq 0 then $
    message, 'No objects found!' $
  else whobj = where(objslits[whobj[0]] eq slitnames, whcnt)
  if whcnt eq 0 then $
    message, 'No objects found!!' $
  else whobj = whobj[0]

; determine the wavelength interval (angstroms per pixel)
  test_data = mrdfits(filelist[whobj], 1, /silent)
  lambda2d = lambda_eval(test_data)
  lamb = lambda2d[*, 0]
  lamb_intvl = (max(lamb)-min(lamb))/n_elements(lamb)
; determine the number of pixels needed to cover the given wavelength
; range. recall that the variables blue_end and red_end denote the
; lower and upper limits of our lambda range, respectively.
  num = ceil((red_end-blue_end)/lamb_intvl)
; make this an even number of pixels.
  num = num/2*2 ; Making it an even number.
; construct the wavelength table for the image.
  lambda = findgen(num)*lamb_intvl + blue_end
; construct a separating bar to clearly mark the edges of slits on the
; image.
  bar = fltarr(num,3) -1000.
; initialize the image. 
  image = bar
; define the variable previous_slitn. it will redefined in the for/do
; loop below.
  previous_slitn = '-1'

; iterate thru the list of slit files and build the image.
  for i=0,nfiles-1 do begin
; read-in the ith slit file.
      if keyword_set(nonlocal) then begin
          fits_info, filelist[i], /silent, n_ext=n_ext
          if n_ext gt 2 then begin
            ss = mrdfits(filelist[i], 3, /silent) 
          endif else begin 
              ss = mrdfits(filelist[i], 1, /silent)
              ss.flux = (ss.flux * 0.)-1
          endelse
      endif else begin 
          ss = mrdfits(filelist[i], 1, /silent)
      endelse
; determine the number of rows in the ith file.
      nrows[i] = (size(ss.flux, /dimen))[1]
; print info to display.
;      print, 'reading file, rows: ', filelist[i], nrows[i]
; get slitnumber from file name, again. this was already done
; before...why do it again???
      filename = strsplit(filelist[i], '.', /extract) 

      slitn = strmid(filelist[i],BRpos-3,3)

      slitnum = round(float(slitn))
; check for NaN's in the slit data. set NaN's to be -1000.
      ind_NaN = where(finite(ss.flux) eq 0, cnt) 
      if cnt ne 0 then ss.flux(ind_NaN) = -1000. 

      if i ne 0 then begin
; if this file is does not match the previous file, then we must put
; the other file into the image.
          if slitn ne previous_slitn then begin
; place the slit number in the "tail" region.
              slit_aligned[(num-tail_l):(num-1), *] = previous_slitnum 
; also, place the object number in the "tail" region so as to denote
; the object position on the image.
              objpos = round(find_objpos(previous_slitnum, nrows[i-1], $
                                         objnum=objnum))
              if objpos[0] gt 0 then begin
                  dex = where(objpos lt nrows[i-1] and objpos ge 0, gdcnt) 
                  if (gdcnt gt 0) and strnumber(objnum[0]) then begin 
                      for qq=0,gdcnt-1 do begin
                          objstr = strcompress(string(objnum[dex[qq]]), /rem)
                          if strlen(objstr) gt 3 then $
                            objnum[dex[qq]] = $
                            long(strmid(objstr, strlen(objstr)-3)) + 1000.
                          slit_aligned[(num-tail_l):(num-1),$
                                       objpos[dex[qq]]] = objnum[dex[qq]]
                      endfor
                  endif
              endif
; add the data to the image and put a bar on after it.
              image = [[image],[slit_aligned],[bar-previous_slitnum]]
; initiate the next slit array.
              slit_aligned = fltarr(num,nrows[i])
          endif else begin
; if we have a pair of slit files (red and blue pair) then place them
; together but don't add them to the image.
; this only accounts for cases where R file (second file) has more
; rows then the B file (first file)!!! this ain't robust!
              if nrows[i] gt nrows[i-1] then $ ;R image has more rows than B's
                for j =0,nrows[i]-nrows[i-1]-1 do $
                slit_aligned=[[slit_aligned],[bar[*,0]-slitnum]] ; Increasing the height
          endelse
      endif else slit_aligned = fltarr(num,nrows[i]) ;i=0,the first slitfile

; get the lambda solution for this slit file.
      lambda2d = lambda_eval(ss)

; test whether wavelength solution has zeros. why???
      test =  where(lambda2d eq 0, test_count)
      if test_count ne 0 then begin
          print, 'Error encounted in slit file: ', filelist[i]
          print, 'Reason: There are zeroes in wavelength solution.'
          print, 'This file was skipped.'
      endif
; determine where to place the data in the image. we want to do this
; so that all the files are properly aligned in lambda space.
      indices = round((lambda2d-blue_end)/lamb_intvl)

; fill an array with the slit data.
      for j=0, nrows[i]-1 do begin
;        slit_aligned[indices[*,j],j] = ss.flux[*,j]
          lamrange=minmax(lambda2d[*,j])
          wh=where(lambda gt lamrange[0] and lambda lt lamrange[1],whct)
          if whct gt 0 then $
            slit_aligned[wh,j]=interpol(ss.flux[*,j],lambda2d[*,j],lambda[wh])
      endfor
; store the slitn and slitnum values before iterating.
      previous_slitn = slitn
      previous_slitnum = slitnum
  endfor

; now, for the last slit file, complete the image array. 
  slit_aligned[(num-tail_l):(num-1), *] = slitnum
  objpos = round(find_objpos(slitnum, nrows[nfiles-1], objnum=objnum))
  if objpos[0] gt 0 then begin
      dex = where(objpos lt nrows[nfiles-1] and objpos ge 0, gdcnt) 
      if gdcnt gt 0 and size(objnum[0], /tname) ne 'STRING' then begin
          for qq=0,gdcnt-1 do begin
              objstr = strcompress(string(objnum[dex[qq]]), /rem)
              if strlen(objstr) gt 3 then $
                objnum[dex[qq]] = long(strmid(objstr, strlen(objstr)-3)) + 1000.
              slit_aligned[(num-tail_l):(num-1),objpos[dex[qq]]] = objnum[dex[qq]]
          endfor
      endif
  endif
  image = [[image],[slit_aligned],[bar-slitnum]] ; add the last slit


; finally, trim the array of blank space from the image.
  print, 'trimming the ends...'
  step = 50
  excluderow = isanobj[image[num-1, *]] eq 0 or $
    (abs(image[num-1, *]) gt 990 AND abs(image[num-1, *]) lt 1300)
  excludearray = (intarr(step)+1) # excluderow


; save the tail information in which we specified the slitno.
  tail = image[num-tail_l:num-1, *]
; trim the left side of the image.
  left = -step
  repeat begin
      left = left+step
      temp = where(image[left:(left+step-1)<(num-1), *] gt -1000. and $
                   image[left:(left+step-1)<(num-1), *] ne 0. and $
                   (excludearray eq 0), cnt)
  endrep until cnt ne 0

; trim the right side of the image.
  right = num-1-tail_l+step
  repeat begin
      right = right-step
      temp = where(image[(right-step+1)>0:right, *] gt -1000. and $
                   image[(right-step+1)>0:right, *] ne 0. and $
                   (excludearray eq 0), cnt)
  endrep until cnt ne 0

; with the variables left and right defined, extract the trimmed image
; from the large image and add the tail info back onto the image.
  if left-tail_l lt 0 then dex1 = 0 $
  else dex1 = left-tail_l
  npix = n_elements(image[*,0])
  if right+tail_l gt npix-1 then dex2 = npix-1 $
  else dex2 = right+tail_l
  image = [tail, image[dex1:dex2, *], tail]

; get the lambda values at the left and right edges.
  lambda_left = left*lamb_intvl+blue_end
  lambda_right = lambda_left+(num-1)*lamb_intvl
; trim the lambda array as we trimmed the image array and pad it with
; zeros so that it has the same length as the image with its tails.
  lamlo = findgen(tail_l)*lamb_intvl + lambda[dex1] - lamb_intvl*tail_l
  lamhi = findgen(tail_l)*lamb_intvl + lambda[dex2] + lamb_intvl*tail_l
  lambda = [lamlo, lambda[dex1:dex2], lamhi]
;Add wavelength indicator on the top
;  lamb_table = findgen(num)*lamb_intvl + lambda_left
;  lamb_table = lamb_table#(fltarr(10)+1)
  lamb_table = lambda # (fltarr(10)+1)

; float compress the image.
  print,'floatcompressing...'
  print, 'It takes a while, please wait.'
  ind_NaN = where(finite(image) eq 0, count)  ;Looking for NaN data
  if count ne 0 then image[ind_NaN] = -1000. 
  image = floatcompress(image, ndig=11)


  image = [[lamb_table], [image]]
; redefine the num variables which gives the number of pixels in the
; wavelength direction.
  num = (size(image, /dimen))[0]
; why make num even?
  if fix(num/2*2) ne num then begin ; if num is an odd number
      image = image[0:num-2, *]
      num = num-1               ; make it an even number
  endif

; smooth and rebin the image.
  if keyword_set(rebin) then begin
      print,'smoothing...'
      image = smooth(image,3,/NAN)
      print,'rebinning...'
      dimen = round(size(image,/dimen)/2.)*2
      bar = fltarr(num) -1000.  ; the length (num) has changed.
      if dimen[1] gt (size(image,/dimen))[1] then image=[[image],[bar]]
      image = rebin(image,dimen[0]/2,dimen[1]/2)
                                ;Side effect on num and lamb_intvl
      num = n_elements(image[*, 0])
      lamb_intvl = (lambda_right-lambda_left)/(num-1.)
  endif

;Write a header for the FITS file
  print, 'Making the header ...'
  mkhdr, hdr, image, /extend
  delt=[lamb_intvl,0]
  CRPIX=[0,0]
  CRVAL=[image[0,0],0]
  CTYPE=['LAMBDA','LAMBDA']
  make_astr,astr,delt=delt,crpix=CRPIX,crval=CRVAL,ctype=CTYPE
  putast,hdr,astr

  return,  image

end


