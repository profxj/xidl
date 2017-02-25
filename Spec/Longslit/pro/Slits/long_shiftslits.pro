FUNCTION LONG_SHIFTSLITS, tset_slits, xshift

nx = tset_slits[0].dims[0]
ny = tset_slits[0].dims[1]

if (size(tset_slits[0].COEFF, /n_dimen) EQ 1) then nslit = 1 $
ELSE nslit = (size(tset_slits[0].COEFF, /dimens))[1]
   
FOR j = 0L, nslit-1L DO BEGIN 
    IF tset_slits[0].coeff[0, j] GT 0.0 THEN $
      tset_slits[0].coeff[0, j] =  tset_slits[0].coeff[0, j] + xshift
    IF tset_slits[1].coeff[0, j] LT (nx-1) THEN $
      tset_slits[1].coeff[0, j] = tset_slits[1].coeff[0, j] + xshift
ENDFOR
IF TAG_EXIST(tset_slits, 'XCORR_COEFF') THEN $
  tset_slits.xcorr_coeff[0, *] = tset_slits.xcorr_coeff[0, *] + xshift

RETURN, tset_slits
END
