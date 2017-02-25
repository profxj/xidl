;+ 
; NAME:
; prs_cldycoll
;  V1.1
;
; PURPOSE:
;    Parses a standard binary fits file into the CLOUDY struct
;   for collisional ionization.  The default file was created by
;   J.C. Howk using Cloudy
;
; CALLING SEQUENCE:
;   
;   prs_cldycoll, stucture, [filename]
;
; INPUTS:
;  [filename]  - Input fits file  
;        [default: $XIDL_DIR/Cloudy/cloud_collisons.fits]
;
; RETURNS:
;   structure      - IDL structure strctcldy
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   03-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro prs_cldycoll, supstrc, infil

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'prs_cldycoll, struct, [infil] (v1.1)'
    return
  endif 

  if keyword_set( DS96 ) then $
    infil = getenv('XIDL_DIR')+'/Cloudy/cloudy_collisions.fits' $
  else infil = getenv('XIDL_DIR')+'/Cloudy/gnat_CIE.fits'

; Read in the Binary Fits file

  tmp1 = xmrdfits(infil,1,/silent)
  ncoll = n_elements(tmp1)
  
; Create the Structure
  
  tmp2 = { strctcldy }
  supstrc = replicate(tmp2, ncoll)


  if keyword_set(DS96) then begin
; Fill it up
      
      for qq=0L,ncoll-1 do begin
          supstrc[qq].NHI = tmp1[qq].H1COLUMN
;  supstrc.FeH = tmp1.FeH
          supstrc[qq].U = tmp1[qq].log_U
          supstrc[qq].nH = tmp1[qq].DENSITY
          supstrc[qq].Jnu = tmp1[qq].log_Jnu
;  supstrc.Spec = tmp1.Spec
;  supstrc.flg = tmp1.flg_cldy
          supstrc[qq].Tval = tmp1[qq].temperature
          
          ;; Ions
          supstrc[qq].X[1,1:2] = tmp1[qq].H[1:2]
          supstrc[qq].X[2,1:3] = tmp1[qq].He[1:3]
          supstrc[qq].X[6,1:7] = tmp1[qq].C[1:7]
          supstrc[qq].X[7,1:8] = tmp1[qq].N[1:*]
          supstrc[qq].X[8,1:9] = tmp1[qq].O[1:*]
          supstrc[qq].X[10,1:11] = tmp1[qq].Neon[1:*]
          supstrc[qq].X[12,1:13] = tmp1[qq].Mg[1:*]
          supstrc[qq].X[14,1:15] = tmp1[qq].Si[1:*]
          supstrc[qq].X[15,1:16] = tmp1[qq].P[1:*]
          supstrc[qq].X[16,1:17] = tmp1[qq].S[1:*]
          supstrc[qq].X[18,1:18] = tmp1[qq].Ar[1:*]
          supstrc[qq].X[20,1:20] = tmp1[qq].Ca[1:*]
          supstrc[qq].X[22,1:22] = tmp1[qq].Ti[1:*]
          supstrc[qq].X[24,1:24] = tmp1[qq].Cr[1:*]
          supstrc[qq].X[25,1:25] = tmp1[qq].Mn[1:*]
          supstrc[qq].X[26,1:26] = tmp1[qq].Fe[1:*]
          supstrc[qq].X[28,1:28] = tmp1[qq].Ni[1:*]
          supstrc[qq].X[30,1:30] = tmp1[qq].Zn[1:*]
      endfor
      
      ;; Sort on Temperature
      srt = sort(supstrc.Tval)
      supstrc = supstrc[srt]
  endif else begin  ;; Gnat & Sternberg 2007
      ;; Temperature
      supstrc.Tval = tmp1.T
          
      ;; Ions
      for qq=0L,ncoll-1 do begin
          supstrc[qq].X[1,1:2] = alog10(tmp1[qq].H)
          supstrc[qq].X[2,1:3] = alog10(tmp1[qq].He)
          supstrc[qq].X[6,1:7] = alog10(tmp1[qq].C)
          supstrc[qq].X[7,1:8] = alog10(tmp1[qq].N)
          supstrc[qq].X[8,1:9] = alog10(tmp1[qq].O)
          supstrc[qq].X[10,1:11] = alog10(tmp1[qq].Neon)
          supstrc[qq].X[12,1:13] = alog10(tmp1[qq].Mg)
          supstrc[qq].X[14,1:15] = alog10(tmp1[qq].Si)
          supstrc[qq].X[16,1:17] = alog10(tmp1[qq].S)
          supstrc[qq].X[26,1:27] = alog10(tmp1[qq].Fe)
      endfor
  endelse
      
return
end
