;+ 
; NAME:
; wfccd_editspec
;    Versioedit.0
;
; PURPOSE:
;   Plots a series of spectra to allow a quick check
;
; CALLING SEQUENCE:
;   
;   wfccd_editspec, wfccd, maskid, expsr, XSIZE=, YSIZE=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_editspec, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_editspec, wfccd, mask_id, expsr, XSIZE=xsize, YSIZE=ysize, $
                    FSPEC=fspec


;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_editspec, wfccd, maskid, [expsr], XSIZE=, YSIZE=,OUTFIL='
    print, '        /FSPEC (v1.0)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then xsize = 1200
  if not keyword_set( YSIZE ) then ysize = 600

  ; WFCCD OBJ
  if not keyword_set(FSPEC) then begin
      
      allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
                     wfccd.mask_id EQ mask_id, nexp)
      if keyword_set(expsr) then exp = allexp[expsr] else exp=allexp[0]

      ;  Object Structure
      wfobj = xmrdfits(wfccd[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)

      ; Tmp arrays
      gdobj = where(wfobj.flg_anly NE 0)
  
      npix = wfobj[gdobj[0]].npix
      wv = wfobj[gdobj].wave[0:npix-1]
      fx = wfobj[gdobj].fx[0:npix-1]
      var = wfobj[gdobj].var[0:npix-1]
      title = strtrim(wfobj[gdobj].slit_id,2)+strtrim(wfobj[gdobj].obj_id,2)

      ; Pass it on
      x_editspec, wv, fx, var, title, NEWVAR=newvar, /block, $
        FLG=wfobj[gdobj].flg_anly, NEWFLG=newflg

      ; SAVE
      wfobj[gdobj].var[0:npix-1] = newvar
      wfobj[gdobj].flg_anly = newflg
      mwrfits, wfobj, wfccd[exp].obj_fil, /create
  endif else begin
;;;;;;;; FSPEC ;;;;;;;;;;      

      wfccd_wrfspec, wffspec, wfccd, /read
      obj = x_getobjnm(wffspec,mask_id)
      
      npix = wffspec[obj].npix
      x_editspec, wffspec[obj].wave[0:npix-1], $
        wffspec[obj].fx[0:npix-1], $
        wffspec[obj].var[0:npix-1], mask_id, NEWVAR=newvar, /block
      ; SAVE
      wffspec[obj].var[0:npix-1] = newvar
      wfccd_wrfspec, wffspec, wfccd
  endelse

  return
end
