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
   dumf = fltarr(7)

   arcpix =  0.774

; Allow for various input files
   close, /all
   case flg of 
       0: begin  ; ANDYs OUTPUT file
           openr, 1, strtrim(mskfil,2)
           readf, 1, dums       ; Dummy line
           readf, 1, duml
           nslits = duml
           ; CREATE SLIT STRUCTURE
           tmp = {mslitstrct}
           slitstr = replicate(tmp, nslits)
           
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
               ; Bottom edge
               slitstr[i].yedg[0] = ypospix - dumf[6]/arcpix 
               ; Top edge
               slitstr[i].yedg[1] = ypospix + dumf[5]/arcpix 
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
   obj = where(slitstr.id LT 9000)
   slitstr[obj].flg_anly = 1
       
  return, slitstr
end
