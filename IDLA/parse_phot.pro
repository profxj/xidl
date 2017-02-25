;+ 
; NAME:
;  parse_phot
;   Version 1.1
;
; PURPOSE:
;    Reads in the BRI ASCII file from DLA photometery to an IDL
;    structure.  Not too useful for the common folk.
;
; CALLING SEQUENCE:
;   
;   parse_phot, file, nobj, supstrc, FLAG=
;
; INPUTS:
;
; RETURNS:
;
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
; PROCEDURES/FUNCTIONS CALLED:

; REVISION HISTORY:
;   
;
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro parse_phot, file, nobj, supstrc, FLAG= flag, ZFILE=zfile

; parse_phot -- Reads in BRI data to a structure

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'parse_phot, filename, nobj, struct, FLAG= (1=BRI), ZFILE='
    return
  endif 

  if not keyword_set( FLAG ) then    flag    = 1
  if keyword_set( ZFILE ) then    flgzph = 1 else flgzph = 0

; Define the structure

  if(flag EQ 1) then begin

      strctnm = {photBRI, $
                 id: 0,$
                 xpix: 0.0,$
                 ypix: 0.0,$
                 B: 0.0,$
                 Bs: 0.0,$
                 R: 0.0,$
                 Rs: 0.0,$
                 I: 0.0,$
                 Is: 0.0,$
                 zphot: 0.0,$
                 zprob: 0.0,$
                 ztemp: 0 }

      nobj = numlines(file)
      supstrc = replicate(strctnm,nobj)

      id = 0
      xpix = 0.0
      ypix = 0.0
      B = 0.0
      Bs = 0.0
      R = 0.0
      Rs = 0.0
      I = 0.0
      Is = 0.0

      close, 1
      openr, 1, file
      for k=0,nobj-1 do begin
          readf, 1, '(i,8f)', id, xpix, ypix, B, Bs, R, Rs, I, Is
          supstrc[k].id = id
          supstrc[k].xpix = xpix
          supstrc[k].ypix = ypix
          supstrc[k].B = B
          supstrc[k].Bs = Bs
          supstrc[k].R = R
          supstrc[k].Rs = Rs
          supstrc[k].I = I
          supstrc[k].Is = Is
      endfor
      close, 1

;  Zphot
      if(flgzph EQ 1) then begin
          temp = 0
          lik = 0.0
          zphot = 0.0
          zprob = 0.0
          
          openr, 1, zfile
          for k=0,nobj-1 do begin
              readf, 1, id, xpix, ypix, zphot, lik, temp, zprob
              if(abs(supstrc[k].xpix - xpix) GT 1) then break
              supstrc[k].zphot = zphot
              supstrc[k].zprob = zprob
              supstrc[k].ztemp = temp
          endfor
          close, 1
      endif

  endif

return
end

