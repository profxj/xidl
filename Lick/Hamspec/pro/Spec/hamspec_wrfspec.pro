;+ 
; NAME:
; x_wrechfspec
;    Version 1.0
;
; PURPOSE:
;    Reads and writes the Hamspec echfspec structure
;
; CALLING SEQUENCE:
;   hamspec_wrfspec, echfspec, outfil, /READ, /UVES
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Sep-2005 Written by JXP
;   24-May-2013 Slightly modified from hires_wrfspec by ENK
;-
;------------------------------------------------------------------------------

pro hamspec_wrfspec, echfspec, outfil, READ=read, UVES=uves

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'hamspec_wrfspec, echfspec, outfil, /READ, /UVES [v1.1]'
    return
  endif 

;  Optional Keywords

  if keyword_set(UVES) then begin
      uves_wrfspec, echfspec, outfil, READ=read
      return
  endif

  if n_elements(ecfspec) GT 1 then stop

  if not keyword_set( READ ) then begin
;;;;;;; WRITE ;;;;;;;;;
;  Create anonymous structure
      anon = { $
               field: ' ', $
               obj_id: ' ',        $ ; ID value (a=primary, b-z=serendip, x=NG)
               flg_anly: 0,      $ ;  0=No analy, 1=Extracted, 2=Fluxed 
               obj_type: 0, $
               xyimg: fltarr(2), $ ; xy pix of original image
               mag: 0., $       ; Usually R mag
               phot_fil: ' ', $
               img_fil: ' ', $  ; Img file
               nexp: 0L, $
               wvmnx: dblarr(100,2), $ ; Wave min/max for each exposure
               texp: dblarr(100), $ ; t_exp for each exposure
               obj_fil: strarr(100), $
               flg_fin: 0, $    ; 
               flg_flux: 0, $   ; 1=fnu, 2=flambda
               npix: lonarr(160), $ ; npix
               phys_ordr: lonarr(160), $ ; ORDERS
               wave: dblarr(6000,160), $
               fx: fltarr(6000,160), $
               var: dblarr(6000,160), $ ; <=0 :: rejected pix
               novar: dblarr(6000,160), $ ; <=0 :: rejected pix
               class:  '', $     ;;;;;; ALL ZANS BELOW HERE  ;;;;;;;;;
               subclass: '', $
               z: 0.0, $
               z_err: 0.0, $
               rchi2: 0.0, $
               dof:  0L, $
               rchi2diff: 0.0, $
               tfile: '', $
               tcolumn: lonarr(20) - 1L, $
               npoly:       0L, $
               theta:      fltarr(20), $
               vdisp:     0.0, $
               vdisp_err:  0.0  $
            }

      copy_struct, echfspec, anon, EXCEPT_TAGS=['zans']
      copy_struct, echfspec.zans, anon

      ; Catch zero strings
      a = where(strlen(anon.field) EQ 0, na)
      if na NE 0 then anon[a].field = ' '
      a = where(strlen(anon.obj_id) EQ 0, na)
      if na NE 0 then anon[a].obj_id = ' '
      a = where(strlen(anon.class) EQ 0, na)
      if na NE 0 then anon[a].class = ' '
      a = where(strlen(anon.subclass) EQ 0, na)
      if na NE 0 then anon[a].subclass = ' '
      a = where(strlen(anon.tfile) EQ 0, na)
      if na NE 0 then anon[a].tfile = ' '
      a = where(strlen(anon.obj_fil[0]) EQ 0, na)
      if na NE 0 then anon[a].obj_fil[0] = ' '
      a = where(strlen(anon.phot_fil) EQ 0, na)
      if na NE 0 then anon[a].phot_fil = ' '
      a = where(strlen(anon.img_fil) EQ 0, na)
      if na NE 0 then anon[a].img_fil = ' '

      ; Write
      mwrfits, anon, outfil, /create
      spawn, 'gzip -f '+outfil
      delvarx, anon
  endif else begin
;;;;; READ ;;;;;;;;
      tmp = {hamspecfspecstrct} 
      anon = xmrdfits(outfil, 1, /silent)
      nanon = n_elements(anon) 
      echfspec = replicate(tmp, nanon)

      ; COPY
      copy_struct, anon, echfspec, /recur_to, except=["WAVE", "FX", "VAR", "NOVAR"]
      if tag_exist(anon[0], 'NOVAR') then flg_nov = 1 else flg_nov = 0
      
      ;; Get the arrays right!
      sz = size(anon.wave, /dimen)
      for ii=0L,nanon-1 do begin
         for jj=0L,sz[1]-1 do begin
            echfspec[ii].wave[0:sz[0]-1,jj] = anon[ii].wave[*,jj]
            echfspec[ii].fx[0:sz[0]-1,jj] = anon[ii].fx[*,jj]
            echfspec[ii].var[0:sz[0]-1,jj] = anon[ii].var[*,jj]
            if flg_nov then $
               echfspec[ii].novar[0:sz[0]-1,jj] = anon[ii].novar[*,jj]
         endfor
      endfor

      delvarx, anon, tmp
  endelse
  return
end
