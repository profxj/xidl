;+ 
; NAME:
; sdss_newzem
;  Version 1.1
;
; PURPOSE:
;  Returns the new QSO redshift if a new one has been calculated.  The
;  redshifts generally come from the Shin et al. 2007 prescription.
;
; CALLING SEQUENCE:
;   zem = sdss_newzem(zmin, zmax, GZFIL=, XZFIL=)
;
; INPUTS:
;   plate -- Plate
;   fiber -- Fiber
;
; RETURNS:
;   zem -- New QSO emission redshift
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; /FMATCH -- Demand that there is a match
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
;   10-Nov-2005 Written by JXP
;-
;------------------------------------------------------------------------------
function sdss_newzem, plate, fiber, FMATCH=FMATCH, FIL=fil, OFIL=ofil

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'z = pdla_fixz(plate, fiber, /FMATCH, FIL= )  [v1.0]'
    return, -1
  endif 

  if not keyword_set(FIL) then $
    fil = getenv('XIDL_DIR')+'/SDSS/General/sdss_dr7qsozem.fits'

  if keyword_set(OFIL) then $
    fil = getenv('XIDL_DIR')+'/SDSS/General/sdss_newqsozem.fits'

  ;; Read
  newqso = xmrdfits(fil,1,/sile)

  ;; Loop
  nqso = n_elements(plate)
  zval = dblarr(nqso)
  for qq=0L,nqso-1 do begin
      ;; Match
      mt = where(plate[qq] EQ newqso.plate and $
                 fiber[qq] EQ newqso.fiberid, nmt)
      if keyword_set(FMATCH) and nmt LT 1 then stop
      ;; Set
      if nmt NE 0 then zval[qq] = newqso[mt[0]].z
  endfor

  ;; Scalar?
  if nqso EQ 1 then zval = zval[0]

  return, zval
end
