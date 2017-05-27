;+
; NAME:
;   deimos_2dreduce
;
; PURPOSE:
;   top level routine to test 2d analysis code, reduces one DEIMOS
;   mask to separate FITS files for each slitlet.
;
; CALLING SEQUENCE:
;   deimos_2dreduce, planfile, chiplist=chiplist, outdir=outdir
; 
; INPUTS:
;   planfile -- an ascii file detailing what flats, arcs, science
;               frames to use for this mask
; 
; KEYWORDS:
;   chiplist -- if set, what subset of chips (1-8) to use, otherwise all
;            8 are reduced.
;   outdir -- directory for output of data (current if not specified)
;
; OUTPUTS:
;   separate FITS files for each slitlet, one HDU per science frame,
;   separate files for R and B sides of spectrum
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;  Doug Finkbeiner, Marc Davis  Aug01-May02
;
;----------------------------------------------------------------------

pro deimos_2dreduce,  planfile,  chiplist=chiplist, outdir=outdir, $
          file=file, quick=quick
; script to call deimos_spslit




; call astrolib if needed
;  defsysv, '!TEXTOUT', exist=exist
;  if NOT exist then astrolib

  if n_elements(quick) eq 0 then quick = 0

  if NOT keyword_set(outdir) then outdir = ''

  device, pseudo=8
  verbset, 4-3*(quick gt 0)
  ybin = 8
  flat = 1
  anamorph = 1.6 ; anamorphic factor, which is a function of grating and tilt
  spline =  0 ; do not do spline fit

; -------- Get environment variables, set paths
  deimos_data = getenv('DEIMOS_DATA')+'/'
  calib_data = getenv('CALIB_DATA')+'/'

  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'
  if calib_data eq '/' then message, 'You need to set $CALIB_DATA!'

  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'You need to set $DEEP_DIR!'

; input the plan file for this mask

  if n_elements(planfile) eq 0 then planfile = '*.plan'
  
  vprint, 2, '------------------------------------------------------------'
  vprint, 2, systime()
  vprint, 2, 'Plan file: ', planfile
  read_planfile, planfile, maskname, rawdatadir, outdatadir, $
             flatnames, arcnames, sciencenames, chips=chiptext



   chipnums=fix(strsplit(chiptext,',',/extract))
   if NOT keyword_set(chiplist) then chiplist=chipnums


  if strpos(maskname, '.E') ge 0 OR strpos(maskname,  '.W') ge 0 then reads, maskname, maskno, format='(I4.4)' else maskno = maskname;get mask_number
  if strpos(rawdatadir, '/') ne 0 then maskdir = deimos_data+ rawdatadir $
    else maskdir = rawdatadir


  if quick gt 0 and n_elements(file) gt 0 then begin
     sciencenames = concat_dir(maskdir, file)
     vprint, 3, 'Data directory: ', maskdir
     vprint, 3, 'Output directory: ', outdatadir
     vprint, 3, 'Science files: '
     vprint, 3, sciencenames
  endif

  calib_flist = findfile(outdatadir+'calib*.fits', count=filecount)
  if filecount eq 0 then $
     calib_flist = findfile(outdatadir+'calib*.fits.gz', count=filecount)

  if filecount eq 0 then $
     message, 'generate the calib files first!!'



; perhaps things like this should be versioned???

  badmaskname = calib_data+'deimos_badmask.fits.Z'
  pixflatname = calib_data+'processed_pix_mult_flat.fits.gz'


  if file_test(badmaskname) eq 0 then begin 
     ; generate file
     vprint, 2, ' Cannot find badmask file - Generating now...'
     spawn, 'mkdir -p '+calib_data
     deimos_badchip, calib_data
     print, ' Done.'
  endif


  chip_calib = intarr(filecount)
  for i=0, filecount -1 do $
    chip_calib[i] = sxpar(headfits(calib_flist[i],  ext=1), 'CHIPNO')  

  for i=0, n_elements(chiplist) -1  do begin 
     chipno = chiplist[i]
  

     jj = where(chip_calib eq chipno, nff)
     if nff gt 0 then begin 
       print
       vprint,1, 'Chipno: ', chipno
       vprint, 1, 'Number of slitlets: ', nff
       fflist = calib_flist[jj]

; -------- read subimage
      if quick le 0 then print, 'Reading bad mask: ', badmaskname
      deimos_badchip = mrdfits(badmaskname, chipno,  /silent) 

      if file_test(pixflatname) eq 0 OR quick gt 0 then begin 
                  vprint, 2, ' Cannot find pixflat file - Setting to 1...'
                  deimos_pixflat_mul = 1.
      endif else begin
	  vprint, 2, 'Using pixflat!'
          deimos_pixflat_mul = mrdfits(pixflatname, chipno, /silent)+1.
;          whzero = where(deimos_pixflat_mul eq 0, zeroct)
;          if zeroct gt 0 then deimos_pixflat_mul[whzero] = 1.
      endelse
 
       deimos_spslit, chipno, maskno, fflist, $
          deimos_badchip, sciencenames, $
         outdatadir=outdatadir, pixmap=deimos_pixflat_mul, /flat, quick=quick
    endif  
  endfor

; verbosity levels

; 0 nothing
; 1  per mask
; 2  per chip
; 3  per slit
; 4  everything!
  return
end









