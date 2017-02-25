;+ 
; NAME:
; imacsls_extract   
;   Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;  imacsls_extract, imacsls, setup, obj_id, side, exp, $
;                 APER=, SKYNORD=, /STD=, /SILENT, $
;                 YMODEL=, /BOXCAR, /CHK
;
; INPUTS:
;   imacsls  -  IMACS structure
;   setup    -  Setup ID value
;   side     -  Specific CCD (1='blue', 2='red')
;   obj_id   -  Object ID value
;   [exp]    -  Expsoure index values
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  APER=   -- Aperture size (boxcar)
;  /BOXCAR -- Do boxcar extraction
;  /CHK    -- 
;  YMODEL= -- 2D solution from extract_image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro imacsls_extract, imacsls, setup, obj_id, side, exp, $
                  APER=aper, SKYNORD=skynord, STD=std, SILENT=silent, $
                  YMODEL=ymodel, BOXCAR=boxcar, CHK=chk


;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'imacsls_extract, imacsls, setup, side, obj_id, [exp], APER=  [v1.1]'
      print, '    SKYNOR=, /STD, /SILENT, YMODEL=, /BOXCAR, /CHK'
    return
  endif 

;  Optional Keywords
  if not keyword_set( SKYNORD ) then skynord = 2
  if not keyword_set( SETUP ) then setup = 0
  if not keyword_set( SIDE ) then side = 1
  if not keyword_set( OBJ_ID ) then obj_id = 0
  if not keyword_set( AOFF ) then aoff = 5.

; Find objects
  if not keyword_set( STD ) then begin
      indx = where(imacsls.flg_anly NE 0 AND imacsls.mode EQ 1 AND $
                   imacsls.side EQ side AND imacsls.setup EQ setup AND $
                   imacsls.obj_id EQ obj_id AND imacsls.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'imacsls_extract: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = where(imacsls.flg_anly NE 0 AND imacsls.mode EQ 1 AND $
                   imacsls.side EQ side AND imacsls.setup EQ setup AND $
                   imacsls.type EQ 'STD', nindx)
      if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
          print,'Syntax - ' + $
            'imacsls_extract, imacsls, setup, side, EXP, /STD   [v1.0]'
          return
      endif 
;      indx = indx[exp]
;      nindx = 1L
      ;; Set boxcar
      BOXCAR=1
      aper = [15.,15.]
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;  Loop
  for q=0L,n_elements(exp)-1 do begin

      ;; Open objfil
      objfil = imacsls[indx[exp[q]]].obj_fil 
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'imacsls_extract: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
      nobj = n_elements(objstr)

      ;; Read IMG+VAR
      imgfil = 'Final/f_'+imacsls[indx[exp[q]]].img_root
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'imacsls_extract: No Final file ', imgfil
          stop
      endif
      print, 'imacsls_extract: Reading in the images..'
      fimg = xmrdfits(imgfil, 2, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz = size(fimg, /dimensions)

      ;; Set Sigma
      if not keyword_set( APER ) then begin
          sigma = fltarr(nobj)
          if not keyword_set( SILENT) then $
            print, 'imacsls_extract: Define the Aperture'
          for kk=0L,nobj-1 do begin
              cline = objstr.xcen
              spec = djs_median(img[cline-15L:cline+15L, *],1)
              aper = x_setaper(spec, objstr[kk].ycen, 0.05, RADIUS=20L)
              ;; STD
              if keyword_set(STD) then aper = aper + 3.
              objstr[kk].aper = aper
              sigma[kk] = total(aper)/4.
          endfor
      endif else objstr.aper = aper

      if not keyword_set( BOXCAR ) then begin ;; OPTIMAL
          ;; TRIM
          img = img[*,20:160]
          var = var[*,20:160]

          ;; IVAR
          ivar = 1. / transpose(var)

          ;; Setup the trace
          xcen = objstr.trace[0:sz[0]-1] - 20.  ; Trim offset

          ;; Kludging!
          if nobj EQ 1 then xcen = [[xcen], [xcen*0.+9999.]]
          
          ;; Optimal
          extract_image, transpose(img), ivar, xcen, $
            sigma, flux, finv, ymodel=ymodel, nPoly=skynord, wfixed=[1,1,1]

          if keyword_set( CHK ) then xatv, ymodel, /block

          ;; Fill up objstr
          for kk=0L,nobj-1 do begin
              objstr[kk].fx[0:sz[0]-1] = flux[*,kk]
              objstr[kk].var[0:sz[0]-1] = 1./finv[*,kk]
              objstr[kk].trace[0:sz[0]-1] = xcen[*,kk]
              objstr[kk].npix = sz[0]-1
          endfor

      endif else begin ;; BOXCAR
          flux = fltarr(sz[1])
          fvar = fltarr(sz[1])
          print, 'imacsls_extract: Boxcar extracting'
          for jj=0L,sz[1]-1 do begin
              flux[jj] = total(fimg[objstr.trace[jj]-objstr.aper[0]:$
                                   objstr.trace[jj]+objstr.aper[1],jj])
              fvar[jj] = total(1./ivar[objstr.trace[jj]-objstr.aper[0]:$
                                   objstr.trace[jj]+objstr.aper[1],jj])
          endfor

          ;; Fill it up
          objstr.fx[0:sz[1]-1] = flux
          objstr.var[0:sz[1]-1] = fvar
          objstr.npix = sz[1]-1
      endelse

      if keyword_set( CHK ) then x_splot, flux, ytwo=sqrt(fvar), /block
          
      ;; Write
      print, 'imacsls_extract: Updating ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil
                     
  endfor


  print, 'imacsls_extract: All Done!'
  return
end
