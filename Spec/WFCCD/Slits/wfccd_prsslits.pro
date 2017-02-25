;+ 
; NAME:
; wfccd_prsslits
;    Version 1.0
;
; PURPOSE:
;    Parses the mask file for slit information
;
; CALLING SEQUENCE:
;   
;   slitstr = wfccd_prsslits( mskfil, [flg] )
;
; INPUTS:
;   mskfil     - Mask file (ASCII) to parse
;
; RETURNS:
;
; OUTPUTS:
;   slittr      -  Creates a structure of slits
;
; OPTIONAL KEYWORDS:
;   flg      - Designates the mask file :: Default = 0 -> Andy's
;              output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   slitstr = wfccd_prsslits( mskfil, [flg])
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function wfccd_prsslits, mskfil, flg

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'slitstr=wfccd_prsslits( mskfil, [flg]) (v1.0)'
    return, -1
  endif 

;  Optional Keywords
   if not keyword_set( flg ) then flg = 0

; Check for mask file
   a = findfile(mskfil, count=nfil)
   if nfil EQ 0 or strlen(mskfil) EQ 0 then begin
       print, 'wfccd_prsslits: Mask file doesnt exist: ', mskfil
       stop
   endif
; Dummy variables

   dums = ''
   duml = 0L
   dumf = fltarr(8)

   arcpix =  0.774

; Allow for various input files
   close, /all
   case flg of 
       0: begin  ; ANDYs OUTPUT file  (circa 2003 or earlier)
           openr, 1, strtrim(mskfil,2)
           readf, 1, dums       ; Dummy line
           readf, 1, duml
           nslits = duml
           ; CREATE SLIT STRUCTURE
           tmp = {mslitstrct}
           slitstr = replicate(tmp, nslits)
;          these are to check vignetting below (MRB/MG 02/04)
;          (note reversal of x/y definition)
           xmask=fltarr(nslits)
           ymask=fltarr(nslits)
           
           for i=0,3 do readf, 1, dums       ; Dummy lines
           for i=0,nslits-1 do begin
               readf, 1, duml, dumf
               slitstr[i].id = duml
               slitstr[i].nobj = 1
               slitstr[i].xpos = dumf[1]/arcpix + 1024
               ypospix = -dumf[0]/arcpix + 1024
               slitstr[i].ypos[0] = dumf[4]  ; Primary science object
               slitstr[i].width = dumf[2]
               slitstr[i].length = dumf[3]
               slitstr[i].priority = dumf[7]
               ; Bottom edge
               slitstr[i].yedg[0] = ypospix - dumf[6]/arcpix 
               ; Top edge
               slitstr[i].yedg[1] = ypospix + dumf[5]/arcpix 
               xmask[i]=dumf[0]
               ymask[i]=dumf[1]
           endfor
           close, 1
           ; Set yedg_flt to yedg
           slitstr.yedg_flt = slitstr.yedg

           ; CONVERT to PIXELS
       end
       1: begin  ; ANDYs OUTPUT file + MRB/MG  (circa 2004 and later)
           openr, 1, strtrim(mskfil,2)
           readf, 1, dums       ; Dummy line (extra dummy line MRB/MG)
           readf, 1, dums       ; Dummy line (extra dummy line MRB/MG)
           readf, 1, dums       ; Dummy line
           readf, 1, duml
           nslits = duml
           ; CREATE SLIT STRUCTURE
           tmp = {mslitstrct}
           slitstr = replicate(tmp, nslits)
;          these are to check vignetting below (MRB/MG 02/04)
;          (note reversal of x/y definition)
           xmask=fltarr(nslits)
           ymask=fltarr(nslits)
           
           for i=0,3 do readf, 1, dums       ; Dummy lines
           for i=0,nslits-1 do begin
               readf, 1, duml, dumf
               slitstr[i].id = duml
               slitstr[i].nobj = 1
               slitstr[i].xpos = dumf[1]/arcpix + 1024
               ypospix = -dumf[0]/arcpix + 1024
               slitstr[i].ypos[0] = dumf[4]  ; Primary science object
               slitstr[i].width = dumf[2]
               slitstr[i].length = dumf[3]
               slitstr[i].priority = dumf[7]
               ; Bottom edge
               slitstr[i].yedg[0] = ypospix - dumf[6]/arcpix 
               ; Top edge
               slitstr[i].yedg[1] = ypospix + dumf[5]/arcpix 
               xmask[i]=dumf[0]
               ymask[i]=dumf[1]
           endfor
           close, 1
           ; Set yedg_flt to yedg
           slitstr.yedg_flt = slitstr.yedg

           ; CONVERT to PIXELS
       end
       else: begin
           print, 'wfccd_prsslits: Not prepared for flg = ', flg
           stop
       end
   endcase

   ; Arcpix
   slitstr.arcpix = arcpix

   ; Field
   slitstr.field = ' '

   ; PA?

   ; Turn off stars!
   obj = where(slitstr.priority NE -1)
   slitstr[obj].flg_anly = 1

; ??? MG determined vignetting limits;
; use these until we have a better solution to the vignetting issue
   ivignette=where((ymask gt  450. and xmask lt -650.) or $
                   (ymask lt -450. and xmask lt -650.) or $
                   (ymask gt  450. and xmask gt  650.) or $
                   (ymask lt -450. and xmask gt  650.) or $
                   (xmask lt -720.) or $
                   (xmask gt 720.), nvignette)
   if(nvignette gt 0) then begin
       splog,'excluding '+string(nvignette)+' slits due to possible vignetting'
       slitstr[ivignette].flg_anly=0
   endif
       
  return, slitstr
end
