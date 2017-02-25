;+ 
; NAME:
; hires_proc   
;     Version 1.0
;
; PURPOSE:
;    Process a raw data frame.  Bias subtract and flatten.  Output is
;    written in Final/ and is a mutli-extension fits file.
;    Extension-0 is the processed image, Extension-1 is the inverse
;    variance.
;
; CALLING SEQUENCE:
;   
;  There are several ways to call hires_proc
; 
;  hires_proc, hires, indx    [processes the image with index=indx]
;  hires_proc, hires, SETUP=1, OBJ=4  [processes images with setup=1 and
;     obj_id=4]
;  rslt = hires_proc_sngl('Rawfil')    [process a single image]
;
; INPUTS:
;   hires     -  HIRES structure
;   [indx]   -  Index (or indices) of the image(s) in the HIRES
;              structure to process
;
; RETURNS:
;
; OUTPUTS:
;  Fully processed image
;
; OPTIONAL KEYWORDS:
;  /CLOBBER - Clobber existing image
;  SETUP    - Setup number
;  EXP      - Exposure numbers
;  OBJ      - Obj_id value
;  /ARC     - Identifies the image as an Arc explicitly.  This is
;             important for bias subtraction on the red chip.
;  /NOFLAT  - Do not flatten
;  FLATFIL  - Name of flat file to use.
;  /SVOV    - Save the OV subtracted image
;  /REDOOV  - Overwrite any existing OV images
;  IOUTFIL  - Name for output file (array for mult
;  /STD     - Processes standard stars corresponding to an input SETUP
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The inverse variance is given this set of values in an attempt to
;   deal with the Poisson limit
;
;   ivar = 1.0/(abs(img * gain - sqrt(2.0)*readno) + readno^2 ) 
;
; EXAMPLES:
;   hires_proc, hires, [20L,21L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_proc_work
;  hires_subbias_sngl
;  hires_getfil
;  hires_delov
;
; REVISION HISTORY:
;   30-Apr-2002 Written by JXP
;   19-Aug-2003 Added work routine
;   20-Feb-2005 Taken from mike_proc
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hires_proc_sngl, rawfil, chip, exten, ARC=arc, $
  REDOOV=redoov, CLOBBER=clobber, SETUP=setup, $
  SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
  IOUTFIL=ioutfil, GAIN=gain, READNO=readno, $
  IMG_FLAT=img_flat, FLAT_IVAR=img_flat_ivar, $
  FRAME=frame, RBIN=rbin, CBIN=cbin, $
  TRIM_CLM=trim_clm

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = hires_proc_sngl(rawfil, chip, exten, /SVOV, /REDOOV, /CLOBBER, ' + $
        ' SETUP=, GAIN=, READNO=, SETUP=, IMG_FLAT= /NOFLAT, [v1.1]'
      return, -1
  endif 

; Optional Keywords
  if not keyword_set( GAIN ) then gain = 1.0
  if not keyword_set( READNO ) then readno = 3.0
  if not keyword_set( SETUP ) and NOT keyword_set( FLATFIL ) and $
    not keyword_set( NOFLAT ) and NOT keyword_set( IMG_FLAT) then begin
      print, 'hires_proc_sngl: Must give info for flat file!!'
      stop
      return, -1
  endif
  if not keyword_set( FRAME ) then begin
      pos = strpos(rawfil, '.fits')
      frame = long(strmid(rawfil,pos-4,4))
  endif
  if not keyword_set(CBIN) or not keyword_set(RBIN) then begin
      head = xheadfits(rawfil, exten=exten)
      cbin = round(2048./sxpar(head,'NAXIS1'))
      if CHIP NE -1 then rbin = round(4096./sxpar(head,'NAXIS2')) $
      else rbin = round(2048./sxpar(head,'NAXIS2'))
  endif
      
  ;; Check for existing image
  flg_proc = 0
  if not keyword_set( CLOBBER ) then begin
      ;; Chk for outfil
      if keyword_set(IOUTFIL) then outfil = ioutfil else stop
;        outfil = hires_getfil('fin_fil', SUBFIL=rawfil, /name)
      if x_chkfil(outfil, /silent) NE 0 then begin
          print, 'hires_proc_sngl: File ', outfil, ' exists, returning...'
          return, -1
      endif
  endif 

  ;; Grab the info
  if chip EQ -1 then CCD='SINGLE'
  hires_strct, hires, FILE_LIST=[rawfil], /nomkdir, /NOFILE, CCD=ccd
      
; BIAS sub if necessary
  print, 'hires_proc_sngl: Bias subtracting'
  rslt = hires_subbias_sngl( rawfil, chip, exten, CLOBBER = REDOOV, $
                             USEBIAS=usebias, BADROW=badrow, $
                             CBIN=cbin, RBIN=rbin, FRAME=frame)

; Flatten
     
  if NOT keyword_set ( NOFLAT )  then begin
    if NOT keyword_set( IMG_FLAT ) then begin
        if not keyword_set( FLATFIL ) then begin
            tmp = {hiresstrct}
            tmp.chip = chip
            tmp.colbin = cbin
            flg_f = hires_getpixflat(hires[(exten-1)>0], FIL=flatfil)
            if flg_f EQ -1 then $
              img_flat = hires_getfil('nqtz_fil', setup, $
                                      CHIP=chip, FIL_NM=flatfil, /name) 
        endif
        print, 'hires_proc:  Using Flat fil ', flatfil
       img_flat = xmrdfits(flatfil, 0, fhead, /silent)
       img_flat_ivar = xmrdfits(flatfil, 1, /silent)
       
   endif

    gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, $
                  complement=badpix, ncomplement=nbad)
  endif
      


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP
  
  ;; Outfil
  if keyword_set(IOUTFIL) then outfil = ioutfil else stop
;    outfil = hires_getfil('fin_fil', SUBFIL=rawfil, /name)
  ;; Read IMG
  img = hires_getfil('ov_fil', CHIP=chip, FRAME=frame, head=head)
  sz_img = size(img, /dimensions)
          
  ;; Create VAR image
;  ivar = 1.0/(abs(img * gain - sqrt(2.0)*readno) + readno^2 ) 
  ivar = 1.0/(abs(img*gain  - 7.) + 7. + readno^2 )   ;; JXP 2/10/06
;  ivar = 1.0/((img*gain  - 7.)>0. + 7. + readno^2 )   ;; JXP January 2010
  ;; Should consider a better method and possibly couple with
                                ; extraction
;  ivar = 1.0/((abs(img * gain)>1.) + readno^2 + 20.) 
          
  ;; Flatten
  img_nrm = img

  if NOT keyword_set(NOFLAT) then begin
      print, 'hires_proc_sngl: Flattening!'
      img_nrm[gdflt] = img[gdflt]  / img_flat[gdflt]
  
      ivar_correction = 0 
      if keyword_set(img_flat_ivar) then begin
          ivar_correction = img_nrm^2 * ivar * $
            ((1./(img_flat_ivar + (img_flat_ivar EQ 0))) > 1.0e-6) $
            * (img_flat_ivar GT 0)
      endif
      ivar[gdflt] = ivar[gdflt] * img_flat[gdflt] * img_flat[gdflt]
      print, 'Mean ivar_correction is ', mean(ivar_correction)
      ivar = ivar / (1 + ivar_correction)
      ;; Need to leave the interorder gaps in!!
;      if nbad NE 0 then ivar[badpix] = 0.
      
  endif else print, 'hires_proc_sngl: Not flattening!'
  delvarx, img

  ;; Deal with badrows
  if keyword_set( BADROW ) then begin
      for ii=0L,n_elements(badrow)-1 do ivar[*,badrow[ii]] = 0.
  endif

  ;; 'Trim' blue (RHS)
  if chip EQ 1 then begin
      if not keyword_set(TRIM_CLM) then trim_clm = 1980L / cbin
      ivar[trim_clm:*,*] = 0.
  endif

  ;; Bad pixel mask for Old chip
  if chip EQ -1 then begin
     hires_badpix_single, badpix
     ;; Rebin
     msk = rebin(badpix, sz_img[0], sz_img[1])
     bad = where(msk GT 0.)
     ivar[bad] = 0.
  endif
     

  ;; Add in gain
  fin_img = temporary(img_nrm)*gain
          
  ;; OUTPUT
  print, 'hires_proc_sngl: Writing... ', outfil
  mwrfits, fin_img, outfil, head, /create, /silent
  mwrfits, ivar, outfil, /silent
  print, 'hires_proc_sngl: Compressing... '
  spawn, 'gzip -f '+outfil
      
  ;; DEL OV
  if not keyword_set( SVOV ) then begin
      ovfil = hires_getfil('ov_fil', CHIP=chip, FRAME=frame, /name)
      spawn, '\rm '+ovfil
  endif

  delvarx, fin_img, var_nrm, img, var, img_nrm

  return, 1
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_proc, hires, indx, SETUP=setup, ARC=arc, CHIP=chip, $
               REDOOV=redoov, CLOBBER=clobber, OBJ=obj, $
               SVOV=svov, FLATFIL=flatfil, NOFLAT=noflat, $
               IOUTFIL=ioutfil, STD=std, EXP=exp, XOFF=xoff, $
                TRIM_CLM=trim_clm

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_proc, hires, [indx], /SVOV, /REDOOV, /CLOBBER, ' + $
        ' CHIP=, OBJ=, SETUP=, /STD, /ARC [v1.1]'
      return
  endif 
  
;  Optional Keywords
   if NOT keyword_set(redoov) then redoov = 0L

   sand = keyword_set(chip) ? long(total(chip)) : 0L
   mfa = hires.flg_anly
   cchip = mfa AND ((keyword_set(chip) ? hires.chip EQ sand : 1) GT 0)

   if keyword_set(obj) then begin 
     objmask = mfa * 0
     for i=0, n_elements(obj)-1 do begin
       any = where(obj[i] EQ hires.obj_id)
       if any[0] NE -1 then objmask[any] = 1
     endfor
     objmask = objmask AND mfa
   endif else objmask = mfa AND $
              ((hires.obj_id GT 0) OR (strtrim(hires.type,2) EQ 'STD'))

   setupmask = mfa AND ((keyword_set(setup) ? $
                         hires.setup EQ setup[0] : 1) GT 0)
   stdmask = mfa AND ((keyword_set(std) ? $
                       strtrim(hires.type,2) EQ 'STD':1) GT 0)


   if keyword_set(obj) OR keyword_set(std) OR keyword_set(indx) EQ 0 then begin
     indx = where( setupmask AND stdmask AND objmask AND cchip, nindx)
     print, 'hires_proc: Using obj# and overriding indx if input'
   endif 

   ;; Exposure
   if keyword_set(EXP) then indx = indx[exp]
   ;; Nindx
   nindx = n_elements(indx)

   if nindx EQ 0 then begin
     print, 'No objects found to process '+ string(fix([total(setupmask), $
                        total(stdmask), total(objmask), total(cchip)])) 
     return
   endif

   fullindx = indx
   allchip = hires[fullindx].chip + 10*hires[fullindx].setup
   allchip = allchip[sort(allchip)]
   uchip = allchip[uniq(allchip)]
   flg_indx = 1
      
   ;; CHIP
  for kk=0L,n_elements(uchip)-1 do begin
      cur_chip = uchip[kk] mod 10
      cur_setup = uchip[kk] / 10
      ;; Single chip
      if cur_chip EQ 9 then begin 
         cur_chip = -1 
         cur_setup = cur_setup + 1
      endif

      indx = where(hires[fullindx].chip EQ cur_chip AND $
                   hires[fullindx].setup EQ cur_setup, nindx)
      if nindx GT 0 then indx = fullindx[indx]
      if nindx EQ 0 then continue

      ;; NIMG
      nimg = n_elements(indx)
      if nimg EQ 0 then begin
        print, 'hires_proc: No files to process. Returning!'
        continue
      endif


;;;;;;;;;;;;;;;;
; Flat
      flg_fset = 0
      if NOT keyword_set( NOFLAT) then begin
          if not keyword_set( FLATFIL ) then begin
              flg_f = hires_getpixflat(hires[indx[0]], FIL=flatfil)
              if flg_f EQ -1 then $
                img_flat = hires_getfil('nqtz_fil', cur_setup, $
                                        CHIP=cur_chip, FIL_NM=flatfil, /name) 
              flg_fset = 1
          endif
          
          print, 'hires_proc:  Using Flat fil ', flatfil
          img_flat = xmrdfits(flatfil, 0, fhead, /silent)
          img_flat_ivar = xmrdfits(flatfil, 1, /silent)
          gdflt = where(img_flat NE 0. AND finite(img_flat) EQ 1, $
                        complement=badpix, ncomplement=nbad)
          hires[indx].flat_fil = flatfil ; Kludge
      endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP

      for ii=0L,nimg-1 do begin
           ;; Chk for outfil
           if keyword_set(IOUTFIL) then outfil = ioutfil[ii] else $
             outfil =hires_getfil('fin_fil', /name, CHIP=hires[indx[ii]].chip,$
                                  CHKFIL=chkf, $
                                  FRAME=hires[indx[ii]].frame)

             if keyword_set(chkf) AND NOT keyword_set(clobber) then begin
                print, 'hires_proc: File ', outfil, ' exists, continuing...'
                hires[indx[ii]].img_final = outfil
                continue
             endif

          ;; Call proc_singl
          rslt = hires_proc_sngl(hires[indx[ii]].rootpth+hires[indx[ii]].img_root,$
                                 cur_chip, hires[indx[ii]].exten, $
                                 IMG_FLAT=img_flat, REDOOV=redoov, $
                                 IOUTFIL=outfil, GAIN=hires[indx[ii]].gain, $
                                 READNO=hires[indx[ii]].readno, /SVOV, $
                                 /clobber, FLAT_IVAR = img_flat_ivar, $
                                 NOFLAT=keyword_set(noflat), $
                                 CBIN=hires[indx[ii]].colbin, $
                                 RBIN=hires[indx[ii]].rowbin, $
                                 FRAME=hires[indx[ii]].frame, $
                                 TRIM_CLM=TRIM_CLM $
                                )
          
          ;; Set xyoff
          if keyword_set(XOFF) then $
            hires[indx[ii]].arc_xyoff[0] = xoff
;          arcs = where(hires.setup EQ cur_setup AND $
;                       hires.chip EQ cur_chip AND $
;                       strtrim(hires.type,2) EQ 'ARC' AND $
;                       hires.flg_anly NE 0, narc)
;          if narc NE 0 then begin
;              diff = abs( (hires[indx[ii]].date + $
;                           hires[indx[ii]].exp/2./3600./24.) - $
;                          (hires[arcs].date + hires[arcs].exp/2./3600./24.) )
;              mn = min(diff, gdarc)
;          endif
;;          diff = abs( (hires[indx[ii]].ut + $
;;                       hires[indx[ii]].exp/2.) - $
;;                      (hires[arcs].ut + hires[arcs].exp/2.) )
;;          b = where( diff GT 21.*3600, nb)
;;          if nb NE 0 then diff[b] = abs(diff[b] - 24.*3600.)
;          mn = min(diff, gdarc)
;          hires[indx[ii]].arc_xyoff = hires[arcs[gdarc]].arc_xyoff
          
          ;; OUTPUT
          hires[indx[ii]].img_final = outfil
      endfor  ;; Image loop
  
      if keyword_set( FLATFIL ) and FLG_FSET then delvarx, flatfil
      ;; DEL OV
      if not keyword_set( SVOV ) then hires_delov, hires, indx
  endfor  ;; CHIP loop

  print, 'hires_proc: All done!'

  return
end


