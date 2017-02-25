;+ 
; NAME:
; kast_flux   
;   Version 1.1
;
; PURPOSE:
;    Flux a set of Kast spectra given a sensitivity file.
;
; CALLING SEQUENCE:
;  kast_flux, kast, setup, side, obj_id, [exp], /STD, SENSFIL=
;
; INPUTS:
;   kast -- Kast IDL structure
;  setup -- Setup value 
;   side -- Side of Kast (blue=1, red=2)
; obj_id -- Object ID value
;  [exp] -- Exposure indices
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /STD  -- Flux a standard
; SENSFIL -- Name of sensitivity file
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
;   12-Nov-2005 added blue g2 (600/4310) default (KLC)
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_flux, kast, setup, side, obj_id, exp, STD=std, SENSFIL=sensfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_flux, kast, setup, side, obj_id, [exp], /STD, SENSFIL=  [v1.1]'
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
          print, 'kast_flux: No images to find obj for!', obj_id
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
      if side EQ 1 then begin
          case strtrim(kast[indx[0]].grising,2) of 
              '452/3306': begin
                  sensfil = getenv('XIDL_DIR')+ $
                    '/Lick/Kast/Calibs/kastsens_g191g1d46.fits'
              end
              '600/4310': begin
                  sensfil=getenv('XIDL_DIR') + $
                          '/Lick/Kast/Calibs/kastsens_g2.fits'
              end 
              else: stop
          endcase
      endif else stop
;          case strtrim(kast[indx[0]].grising,2) of 
;              '300/7500': begin
;                  templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_3007500.idl'
;              end
;              else: stop
;          endcase
;      endelse
  endif
  ;; Sensitivity function
  bset = xmrdfits(sensfil, 1, /silent)

;  Loop
  for q=0L,n_elements(exp)-1 do begin

      ;; Open objfil
      objfil = kast[indx[exp[q]]].obj_fil 
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'kast_flux: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      nobj = n_elements(objstr)
      
      for kk=0L,nobj-1 do begin
          npix = objstr[kk].npix
          if npix eq 0 then continue
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
      print, 'kast_flux: Updating ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil

  endfor


  print, 'kast_flux: All Done!'
  return
end
