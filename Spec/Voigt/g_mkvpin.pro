;+ 
; NAME:
; vpparse
;   Version 1.1
;
; PURPOSE:
;   given a vpstrct, creates a fort.13 input file for vpfit
; 
;
; CALLING SEQUENCE:
;   mkvpin, vpstr, FIL=
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
;   Written by GEP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro g_mkvpin, vpstr, FIL=fil

    star = '*'
    sp = ' '
    bsp = '   '
    spcol = replicate(sp, 2+vpstr.nreg+vpstr.nion)
    line = strarr(vpstr.nreg+2+vpstr.nion)
   
    line[0] = star+bsp
    line[vpstr.nreg+1] = star+bsp
    for i=0,vpstr.nreg-1 do $
        line[i+1] = sp + vpstr.fluxfil[i] + bsp + '1'+ bsp+$
                    strtrim(string(vpstr.reg_beg[i]),2)+$
                    bsp+strtrim(string(vpstr.reg_end[i]),2)

    for i=0,vpstr.nion-1 do $
        line[i+vpstr.nreg+2] = sp+sp+vpstr.ion[i]+bsp+strtrim(string(vpstr.n[i]),2)+bsp+$
                    strtrim(string(vpstr.z[i]),2)+bsp+strtrim(string(vpstr.b[i]),2)+bsp+$
                    '0.0   1.00E+00   0   !  ' + strtrim(string(i+1),2)

    writecol, fil, line, spcol

end
