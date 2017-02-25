;+ 
; NAME:
; imacsls_flux   
;   Version 1.1
;
; PURPOSE:
;    Flux an IMACS 1D spectrum
;
; CALLING SEQUENCE:
;  imacsls_flux, imacsls, setup, obj_id, side, [exp], /STD, SENSFIL=
;
; INPUTS:
;  imacsls -- IMACS structure
;  setup   -- Setup ID value
;  obj_id  -- Object ID value
;  side    -- Refers to CCD (1=blue, 2=red)
;  [exp]   -- Exposure index number
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SENSFIL  -- Sensitivity file
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_flux, imacsls, 0, 1, 0
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_flux, imacsls, setup, obj_id, side, exp, STD=std, SENSFIL=sensfil

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'imacsls_flux, imacsls, setup, obj_id, side, [exp]  [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( NFIN ) then nfin = 5L
  if not keyword_set( SETUP ) then setup = 0
  if not keyword_set( SIDE ) then side = 1
  if not keyword_set( OBJ_ID ) then obj_id = 0

; Find objects
  if not keyword_set( STD ) then begin
      indx = where(imacsls.flg_anly NE 0 AND imacsls.mode EQ 1 AND $
                   imacsls.side EQ side AND imacsls.setup EQ setup AND $
                   imacsls.obj_id EQ obj_id AND imacsls.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'imacsls_flux: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = where(imacsls.flg_anly NE 0 AND imacsls.mode EQ 1 AND $
                   imacsls.side EQ side AND imacsls.setup EQ setup AND $
                   imacsls.type EQ 'STD', nindx)
      if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
          print,'Syntax - ' + $
            'imacsls_flux, imacsls, setup, side, EXP, /STD   [v1.0]'
          return
      endif 
      indx = indx[obj_id]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

; Grab sensititivty file
  if not keyword_set( SENSFIL ) then begin
      if side EQ 1 then begin
          case strtrim(imacsls[indx[0]].grising,2) of 
              '300l': begin
                  sensfil = getenv('XIDL_DIR')+'/Magellan/IMACS/LSlit/Calibs/sensB_ltt9491.fits'
              end
              else: stop
          endcase
      endif else begin
          case strtrim(imacsls[indx[0]].grising,2) of 
              '300l': begin
                  sensfil = getenv('XIDL_DIR')+'/Magellan/IMACS/LSlit/Calibs/sensR_ltt9491.fits'
              end
              else: stop
          endcase
      endelse
  endif
  ;; Sensitivity function
  bset = xmrdfits(sensfil, 1, /silent)

;  Loop
  for q=0L,n_elements(exp)-1 do begin

      ;; Open objfil
      objfil = imacsls[indx[exp[q]]].obj_fil 
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'imacsls_extract: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      nobj = n_elements(objstr)
      
      for kk=0L,nobj-1 do begin
          npix = objstr[kk].npix
          wave = objstr[kk].wave[0:npix-1]
          ;; Sens
          sens = bspline_valu(wave,bset)
          
          ;; Dwv
          if wave[1] GT wave[0] then begin
              dwv = (shift(wave,-1)-shift(wave, 1))/2.
              dwv[0] = wave[1]-wave[0]
              dwv[npix-1] = wave[npix-1]-wave[npix-2]
          endif else begin
              dwv = (shift(wave,1)-shift(wave, -1))/2.
              dwv[0] = wave[0]-wave[1]
              dwv[npix-1] = wave[npix-2]-wave[npix-1]
          endelse

          ;; Flux
          objstr[kk].flux[0:npix-1] = objstr[kk].fx[0:npix-1] * sens / $
            imacsls[indx[exp[q]]].exp / dwv
          objstr[kk].sig[0:npix-1] = sqrt(objstr[kk].var[0:npix-1]) * sens / $
            imacsls[indx[exp[q]]].exp / dwv
          objstr[kk].flg_flux = 2
      endfor

      ;; Write
      print, 'imacsls_wavesol: Updating ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil

  endfor


  print, 'imacsls_extract: All Done!'
  return
end
