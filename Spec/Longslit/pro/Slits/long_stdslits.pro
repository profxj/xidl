;+
; NAME:
;   long_stdslits
;
; PURPOSE:
;   Find slits for a standard star?
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------
PRO long_stdslits, stdfile, slitfile_old, slitfile_new, STDID = STDID

t0 = systime(1)
IF NOT KEYWORD_SET(STDID1) THEN STDID = 0 $
ELSE STDID = STDID1-1L
slitmask_old = mrdfits(slitfile_old, 0)
tset_old = mrdfits(slitfile_old, 1)
traceset2xy, tset_old[0], yy1, leftedge
traceset2xy, tset_old[1], yy2, rightedge
if (size(leftedge, /n_dimen) EQ 1) then nslit = 1 $
else nslit = (size(leftedge, /dimens))[1]
obj = mrdfits(stdfile, 5)
std = obj[STDID]
mask = (std.IVAR_OPT GT 0.0)
func = 'legendre'
ncoeff = 3
xy2traceset, std.ypos, std.xpos, stdset $
             , func = func, ncoeff = ncoeff $
             , invvar = float(mask), yfit = fitstd

tset_slits = tset_old

FOR islit = 0L, nslit-1L DO BEGIN
    tset_slits[0].COEFF[1:*, islit] = stdset.COEFF[1:*]
    tset_slits[1].COEFF[1:*, islit] = stdset.COEFF[1:*]
    djs_iterstat, fitstd-leftedge[*, islit], sigrej = 2.0 $
                  , median = left_med, mask = mask
    djs_iterstat, rightedge[*, islit] - fitstd, sigrej = 2.0 $
                  , median = right_med, mask = mask
    tset_slits[0].coeff[0, islit] = stdset.COEFF[0] - left_med
    tset_slits[1].coeff[0, islit] = stdset.COEFF[0] + right_med 
ENDFOR
traceset2xy, tset_slits[0], yy1, xx1
traceset2xy, tset_slits[1], yy2, xx2
slitmask = long_slits2mask(tset_slits, /silent, nslit = nslit1)

; Write output file
splog, 'Writing output file'
mwrfits, slitmask, slitfile_new, hdr, /create
mwrfits, tset_slits, slitfile_new

splog, 'Number of slits = ', nslit1
splog, 'Elapsed time = ', systime(1)-t0, ' sec'




END
