;+ 
; NAME:
; parse_rao06
;
; PURPOSE:
;    Calculates f(N,X) for a set of DLA
;
; CALLING SEQUENCE:
;
; INPUTS:
;  GZSTR= -- Structure summarizing g(z) for the quasars
;  GZFIL= -- Filename of the g(z) file for the quasars
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /DR3  -- Restrict to DR3
;  /NODBL -- Do not solve for double power-law
;  CL=   -- Confidence limits for error calculations
;  ZBIN= -- redshift interval [default: 2.2 to 5.5]
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
;   Jun-2011 Written by JXP  (based on fig_fndla)
;-
;------------------------------------------------------------------------------
; fig_fits
pro parse_rao06, strct 

  ;; Get structure if necessary
  if not keyword_set( DATFIL ) then datfil = $
    getenv('XIDL_DIR')+'DLA/Surveys/rao06_table.dat'
  tmp = {dlanoestruct}

  nlin =numlines(datfil)-110
  close, /all
  openr, 11, datfil

  dla = replicate(tmp,nlin)
  msk = bytarr(nlin)

  lin = ''
  for qq=1L,110 do readf, 11, lin
  for qq=0L,nlin-1 do begin
      readf, 11, lin
      ;; Parse name
      dla[qq].qso = strmid(lin,0,9)

      ;; Parse mag
      dla[qq].qso_mag = float(strmid(lin,10,4))
      
      ;; Parse zem
      dla[qq].qso_zem = float(strmid(lin,15,5))

      ;; zabs
      dla[qq].zabs = float(strmid(lin,21,6))

      ;; NHI
      dla[qq].NHI = float(strmid(lin,92,5))
      dla[qq].sigNHI[0] = float(strmid(lin,100,4))
      dla[qq].sigNHI[1] = float(strmid(lin,105,4))
  endfor

  gd = where(dla.NHI GE 20.299)
  strct = dla[gd]

  close, /all

  
  return
end
      
