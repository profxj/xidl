;+ 
; NAME:
; kast_flux   
;   Version 1.0
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   spec = x_apall(ydat, [head])
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   wave       - wavelength array
;   DISPLAY    - Display the sky subtracted image with xatv
;   OVR        - String array for ov region:  '[2050:2100, *]'
;   ERROR      - Variance array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_flux, kast, 0, 1, 0
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_flux, kast, setup, side, obj_id, exp, STD=std, SENSFIL=sensfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_flux, kast, setup, side, obj_id, [exp]  [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( NFIN ) then nfin = 5L
  if not keyword_set( SETUP ) then setup = 0
  if not keyword_set( SIDE ) then side = 1
  if not keyword_set( OBJ_ID ) then obj_id = 0

; Find objects
  if not keyword_set( STD ) then begin
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.obj_id EQ obj_id AND kast.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'kast_wavesol: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.type EQ 'STD', nindx)
      if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
          print,'Syntax - ' + $
            'kast_extract, kast, setup, side, EXP, /STD   [v1.0]'
          return
      endif 
      indx = indx[obj_id]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

; Grab sensititivty file
  if not keyword_set( SENSFIL ) then begin
      stop
      if side EQ 1 then begin
          case strtrim(kast[indx[0]].grising,2) of 
              '452/3306': begin
                  templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_g1.idl'
              end
              else: stop
          endcase
      endif else begin
          case strtrim(kast[indx[0]].grising,2) of 
              '300/7500': begin
                  templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_3007500.idl'
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
      objfil = kast[indx[exp[q]]].obj_fil 
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'kast_extract: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      nobj = n_elements(objstr)
      
      for kk=0L,nobj-1 do begin
          npix = objstr[kk].npix
          ;; Sens
          sens = bspline_valu(objstr[kk].wave[0:npix-1],bset)
          ;; Flux
          objstr[kk].flux[0:npix-1] = objstr[kk].fx[0:npix-1] * sens / $
            kast[indx[exp[q]]].exp
          objstr[kk].sig[0:npix-1] = sqrt(objstr[kk].var[0:npix-1]) * sens / $
            kast[indx[exp[q]]].exp
          objstr[kk].flg_flux = 2
      endfor

      ;; Write
      print, 'kast_wavesol: Updating ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil

  endfor


  print, 'kast_extract: All Done!'
  return
end
