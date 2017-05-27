;+
;
; NAME
;     make_bintab_file.pro
;
; PURPOSE
;     Reads the bin tables from a raw DEIMOS frame and writes the
;     tables to a new .fits file. It also reads-in the pcat .fits file
;     corresponding to the given mask and stores the pcat structure as
;     the last extension in the new .fits file.
;
;
; SYNTAX
;     make_bintab_file, planfile
;
; INPUTS
;     planfile = the name of a file which is assumed to be in current
;                directory unless the full path is specified.
;
; KEYWORDS
;     None.
;
; OUTPUTS
;     A file named "xxxx.bintabs.fits" where xxxx denotes the mask
;     number (e.g. 2101 or 1145). This file contains (as .fits
;     extensions) the bin tables that were read from the raw DEIMOS
;     frame. It also contains the pcat data structure. The file is
;     written to the current directory (for the DEEP spec2d reductions
;     this directory should be the same directory as that which will
;     contain the reduced 2d slit files. 
;
; PROCEDURES CALLED 
;     read_planfile.pro
;     getphot.pro
;     fits_open
;     mwrfits
;     headfits
;     sxaddpar 
;
; EXAMPLES
;
;
; COMMENTS
;     This routine is part of the spec2d package. It is intended to be
;     called by the routine domask.pro. Note that this routine assumes
;     that the raw DEIMOS frame pointed to by the planfile is UNcompressed. 
;
; HISTORY
;     Created July 18, 2002 by mcc.
;     Revised July 25, 2002 by mcc - added to routine so that the
;        corresponding pcat structure is also read-in and stored in
;        the bintabs file. 
;     Revised August 12, 2002 by mcc - revised so that entire pcat
;        file is not saved in the xxxx.bintabs.fits file. Instead only
;        the needed subset of the photometry catalog is included.
;     Revised August 26, 2002 by mcc - revised routine so that a
;        portion of the 1HS mask file (1HS.xxxx.fits) is saved to the
;        bin table file (xxxx.bintabs.fits). In particular, the
;        ObjectList bin table from the mask file is saved so as to
;        provide information regarding the position of multiple
;        objects in a single slit --- This information is employed
;        in the do_extract routine (in the 1-d extraction portion of
;        the pipeline.  
;     Revised October 21, 2002 by mcc - correct for bug in mask design
;        (switching of topdist and botdist tags in slitobjmap bin
;        table).
;     Revised December 10, 2002 by mcc - removed keyword
;     /nophot. Routine is now smart enough to check if the mask is a
;     DEEP2 1HS mask before attempting to locate pcat information.
;-

pro make_bintab_file, planfile,quick=quick

  if n_elements(quick) eq 0 then quick = 0

; if the plan file isn't passed then search for it.
  if n_params() eq 0 then begin
    planfile = findfile('*.plan', count=nplans)
    if nplans gt 0 then planfile = planfile[0] $
      else message, '(make_bintab_file.pro) ERROR: no planfile found!'
  endif
; read all of the smack in the .plan file.
  read_planfile, planfile, maskname, rawdatadir, outdatadir, $
    flatnames, arcnames, sciencenames

; select a raw frame from which to read the bin tables.
  file = sciencenames[0]

 if strlen(file) eq 0 then file = arcnames[0]

; construct the name of the output .fits bin tab file.
  bin_name = STRCOMPRESS(maskname, /rem)
  deimos_isdeep, isdeep, bin_name, file=file
  bin_file = bin_name + '.bintabs.fits'
; open the science frame .fits file.  
  fileheader=headfits(file)
  fits_open, file, fcb
; determine which extensions are bin tables.
  dex = where(fcb.xtension eq 'BINTABLE' or $
              fcb.xtension eq 'TABLE', tabcnt)
; if there are no bin table extensions appended to the science frame
; then check on the arc. this is useful for cases where a slitmask
; wasn't used in observing the science target but used in doing the
; arcs and flats.
  if tabcnt eq 0 then begin 
      fits_close, fcb
      file = arcnames[0]
      fits_open, file, fcb
      dex = where(fcb.xtension eq 'BINTABLE' or $
                  fcb.xtension eq 'TABLE', tabcnt)
  endif
; catch possible error in mask design!!! we must determine the
; creation date of the mask so as to adjust for a bug in Drew's mask
; design code.
  maskdex = where(fcb.extname eq 'MaskDesign', maskcnt)
  if maskcnt eq 0 then begin
    fits_close, fcb
    message, 'MaskDesign bin table not found!' 
  endif else maskdex = maskdex[0]
  masktable = mrdfits(file, maskdex, /silent)
  date = (strsplit(masktable.DESDATE, 'T', /extract))[0] ;date has form 'yyyy-mm-dd'
  subdate = strsplit(date, '-', /extract)
  print, '(make_bintab_file.pro) Printing mask design date: ' + $
    strcompress(date, /rem)
  if n_elements(subdate) ne 3 then begin
    print, '(make_bintab_file.pro) Error in extraction of mask design date!!' 
    ddumb = 0
  endif else begin
; all masks made after the date July 25, 2002 have their botdist and
; topdist mixed up! 
      if ( long(subdate[0]) eq 2002 and long(subdate[1]) ge 7 ) or $
        ( long(subdate[0]) ge 2003 ) then begin
          if long(subdate[0]) eq 2002 and $
            long(subdate[1]) eq 7 and $
            long(subdate[2]) lt 25 then ddumb = 0 $
          else ddumb = 1
      endif else ddumb = 0
  endelse

; write science frame header to bintab file  
  mwrfits,junk,bin_file,fileheader,/silent,/create

; now read the bin tables from the raw frame and write them to another
; .fits file which will be stored in the same directory as the reduced
; data.
  for i=0,tabcnt-1 do begin 
; read the ith bin table and the corresponding header.
    bintab = mrdfits(file, dex[i], hdr, /silent) 
; parse the header for the extension name.
    extname = strcompress(sxpar(hdr, 'EXTNAME'), /rem)
; FIX ERROR WITH TOPDIST and BOTDIST!
    if extname eq 'SlitObjMap' and ddumb then begin
      print, '(make_bintab_file.pro) Correcting error in mask design...'
      print, '(make_bintab_file.pro) Switching topdist and botdist.'
      tmpdist = bintab.topdist
      bintab.topdist = bintab.botdist
      bintab.botdist = tmpdist
    endif
; write the bin table and header to the output file
; (xxxx.bintabs.fits).

      mwrfits, bintab, bin_file, hdr, /silent
  endfor

; now get the 'ObjectCat' bin table.
  objdex = where(fcb.extname eq 'ObjectCat', cnt)
  if cnt eq 0 then begin
    fits_close, fcb
    message, 'ObjectCat bin table not found!' 
  endif else objdex = objdex[0]
  objcat = mrdfits(file, objdex, /silent)
  
; check if this is a DEEP2 1HS mask. if it is not, then do not include
; the pcat files in the bintab file and do not include the pcat files
; in the bintab file and do not include the 1HS.xxxx.mask.fits file
; either.
  if isdeep AND quick le 0 then begin
; read the first object number taking only the first 2 digits in the
; string - this denotes the photometry field.
      obj0 = string(objcat[0].object, format='(I8.8)')
      field = strmid(obj0, 0, 2)
; check if the mask is from the Groth strip. if so, get pcat 10.
      if strmid(field, 0, 1) eq '1' then field = '10'
; check if the mask is from the Keck Telescope Redshift Survey
; (KTRS). if so, get pcat 00.
      if strmid(field, 0, 1) eq '0' then field = '00'
; now read the corresponding pcat file.
      pcat = getphot(field, fieldlabel, hdr)
; modify the header to give it a proper extension name.
      sxaddpar, hdr, 'EXTNAME', 'PhotCat', AFTER='XTENSION'
; only take the portion of the pcat needed for this mask. this will
; help keep the xxxx.bintab.fits file small in size.
      object_nums = long(objcat.object)
      n =n_elements(object_nums)
      for i=0,n-1 do begin
          dex = where( pcat.objno eq object_nums[i], cnt )
          if cnt gt 0 then begin
              if i eq 0 then pdex = dex $
              else pdex = [pdex,dex]
          endif
      endfor
      pcat = pcat[pdex]
; now, write the pcat structure as the next extension in the bin table
; .fits file. 
      mwrfits, pcat, bin_file, hdr, /silent
      
; lastly, let's read the 1HS mask file and store a portion of it in
; the bin table file.
      maskfile = findfile('1HS.*.fits', count=nfiles)
; if the maskfile is found in the current directory, then read the
; object list bin table and save it to the bin tab file.
      if nfiles gt 0 then begin
          maskfile = maskfile[0]
          objlist = mrdfits(maskfile, 2, hdr, /silent)
          mwrfits, maskfile, bin_file, hdr, /silent
      endif
  endif
  
; close the fcb.
  fits_close, fcb


end





