;+ 
; NAME:
; wfccd_objsumm
;    Version 1.0
;
; PURPOSE:
;    Write a wfccd_struct to an ASCII file
;
; CALLING SEQUENCE:
;   
;  write_wfccdstr, struct, ANONLY=, OUTFIL=, FITS=
;
; INPUTS:
;   struct        - A wfccd_struct
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   outfil = Output file (default is image.list)
;
; COMMENTS:
;
; EXAMPLES:
;   write_wfccdstr, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_objsumm, wfccd, mask_id, exp_id, OUTFIL=outfil, ANONLY=anonly

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        print, 'wfccd_objsumm, wfccd, mask_id, exp_id [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

;;;;;;;;;;
; Open Obj file
  wfobj = xmrdfits(wfccd[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)
;;;;;;;;;;
; Open Slit file
  wfslit = xmrdfits(wfobj[0].slit_fil, 1, STRUCTYP='mslitstrct', /silent)

  if not keyword_set( OUTFIL ) then begin
      ipos = strpos(wfccd[exp].obj_fil, '.')
      outfil = strmid(wfccd[exp].obj_fil, 0, ipos)+'.txt'
  endif

  nobj = n_elements(wfobj)
  
  close, /all
  openw, 1, outfil

  for q=0,nobj-1 do begin
      
      ; Find the slit
      slit = where(wfslit.id EQ wfobj[q].slit_id)
      
      a = where(wfobj[q].var GT 0., na)
      if na NE 0 then mnwv = min(wfobj[q].wave[a], max=mxwv) else begin
          mnwv = 0.
          mxwv = 0.
      endelse
      
      if (not keyword_set( ANONLY ) OR wfobj[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i3,1x,a12,1x,i1,1x,i6,1x,a2,2i5,i2,2f5.2,2f9.2,i2,1x,a14)',$
        q, $
        wfobj[q].field, $
        wfobj[q].flg_anly, $
        wfobj[q].slit_id, $
        wfobj[q].obj_id, $
        round(wfslit[slit].yedg_flt), $
        wfobj[q].flg_aper, $
        wfobj[q].aper, $
        mnwv, mxwv, $
        wfobj[q].flg_flux, $
        wfobj[q].slit_fil
  endfor
  close, 1

  return

end
