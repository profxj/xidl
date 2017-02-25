;+ 
; NAME:
; prs_cldygrid
;  V1.1
;
; PURPOSE:
;    Parses a standard binary fits file into the CLOUDY struct
;    for Photoionized models
;
; CALLING SEQUENCE:
;   
;   prs_cldygrid, stucture, filename
;
; INPUTS:
;   filename   - Cloudy File
;   ROOT=  - Path to Cloudy files [default: /u/xavier/Cloudy/Grid/Output]
;
; RETURNS:
;   structure  - IDL structure strctcldy containing the Cloudy info
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
;   prs_cldygrid, struct, 'HM01A.fits'
;
;
; PROCEDURES CALLED:
;  xmrdfits
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro prs_cldygrid, supstrc, infil, ROOT=root, STRUCT=struct

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'prs_cldygrid, struct, infil, ROOT=, /STRUCT (v1.0)'
    return
  endif 

;
  if not keyword_set( ROOT ) then root = '/u/xavier/Cloudy/Grid/Output/'

  file = root+infil

; Read in the Binary Fits file

  tmp1 = xmrdfits(file,1,/silent)
  if keyword_set(STRUCT) then begin
      supstrc = tmp1
      return
  endif

; Create the Structure

  tmp2 = { strctcldy }
  supstrc = replicate(tmp2,n_elements(tmp1))

; Fill it up

  supstrc.z = tmp1.z
  supstrc.NHI = tmp1.NHI
  supstrc.FeH = tmp1.FeH
  supstrc.U = tmp1.U
  supstrc.nH = tmp1.nH
  supstrc.Jnu = tmp1.Jnu
  supstrc.Spec = tmp1.Spec
  if tag_exist(tmp1, 'FLG_CLDY') then supstrc.flg = tmp1.flg_cldy $
  else supstrc.flg = tmp1.flg

; Ions

  cnt = 0
  for i=1,30 do begin
      for j=1,8 do begin
          supstrc.X[i,j] = tmp1.X[cnt]
          cnt = cnt+1
      endfor
  endfor

  delvarx, tmp1

return
end
