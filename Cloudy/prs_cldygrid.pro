;+ 
; NAME:
; prs_cldygrid
;  V1.0
;
; PURPOSE:
;    Parses a standard binary fits file into the CLOUDY struct
; CALLING SEQUENCE:
;   
;   prs_cldygrid, stucture, filename
;
; INPUTS:
;   filename       - File
;
; RETURNS:
;   structure      - IDL structure strctcldy
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ION - Input ionic column densities
;  NOELM - Supress inputting Elemental values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_dlalst, struct, '/home/xavier/DLA/Lists/tot_dla.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro prs_cldygrid, supstrc, infil, ROOT=root

; prs_cldygrid -- Reads in CLOUDY grid

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'prs_cldygrid, struct, infil, ROOT= (v1.0)'
    return
  endif 

;
  if not keyword_set( ROOT ) then root = '/u/xavier/Cloudy/Grid/Output/'

  file = root+infil

; Read in the Binary Fits file

  tmp1 = mrdfits(file,1)

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
  supstrc.flg = tmp1.flg_cldy

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
