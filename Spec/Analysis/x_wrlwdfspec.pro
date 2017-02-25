;+ 
; NAME:
; x_wrlwdfspec
;    Version 1.0
;
; PURPOSE:
;    Reads and writes the lwdfspec structure
;
; CALLING SEQUENCE:
;   x_wrlwdfspec, lwdfspec, outfil, /READ
;
; INPUTS:
;   lwdfspec -- IDL structure
;   outfil -- Name of FITS file
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /READ -- Read instead of write the structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_wrlwdfspec, lwdfspec, 'Extract/2015+657_lwd.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_wrlwdfspec, lwdfspec, outfil, READ=read

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'x_wrlwdfspec, lwdfspec, outfil, /READ [v1.0]'
    return
  endif 

;  Optional Keywords

;;;;;;; WRITE ;;;;;;;;;

  if not keyword_set( READ ) then begin
;  Create anonymous structure
      nlwdfspec = n_elements(lwdfspec)
      tmp = { $
              field: ' ', $
              slit_id: 0L, $
              obj_id: ' ',        $ ; ID value (a=primary, b-z=serendip, x=NG)
              flg_anly: 0,      $ ;  0=No analy, 1=Extracted, 2=Fluxed 
              obj_type: 0, $
              xyimg: fltarr(2), $ ; xy pix of original image
              mag: 0., $        ; Usually R mag
              phot_fil: ' ', $
              img_fil: ' ', $
              nexp: 0L, $
              wvmnx: fltarr(100,2), $ ; Wave min/max for each exposure
              texp: dblarr(100), $ ; t_exp for each exposure
              obj_fil: strarr(100), $
              npix: 0L, $
              flg_flux: 0, $    ; 1=fnu
              wave: fltarr(5000), $
              fx: fltarr(5000), $
              var: dblarr(5000), $ ; <=0 :: rejected pix
              class:  '', $     ; ALL ZANS BELOW HERE
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

      anon = replicate(tmp, nlwdfspec)
      copy_struct, lwdfspec, anon, EXCEPT_TAGS=['zans']
      copy_struct, lwdfspec.zans, anon

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
      delvarx, anon, tmp
  endif else begin
;;;;; READ ;;;;;;;;
      tmp = {lwdfspecstrct} 
      anon = mrdfits(outfil, 1, /silent)
      lwdfspec = replicate(tmp, n_elements(anon))

      ; COPY
      copy_struct, anon, lwdfspec, /recur_to
      delvarx, anon, tmp
  endelse
  return
end
