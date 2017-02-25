;+ 
; NAME:
; sdss_chkbal
;    Version 1.1
;
; PURPOSE:
;  used to select the QSOs that have already been checked 
;  for BALs (in 'oldqalfil') and update their values (in 'newqalfil').  
;  The 'matchfil' is then fed to sdss_chkbal and those 
;  QSOs are skipped.
; 
; CALLING SEQUENCE:
;  sdss_balupdate, oldqalfil, newqalfil, matchfil, SDSSPATH=sdsspath
;
; INPUTS:
;  oldqalfil  -- Name of FITS file containing the previous QAL structure
;  newqalfil  -- Name of FITS file containing the current QAL structure
;  
; RETURNS:
;
; OUTPUTS:
;  newqalfil --  With the 'flg_bal' info modified accordingly
;  matchfil  --  Name of FITS file containing array of 'newqsofil'
;                indices that matched
;
; OPTIONAL KEYWORDS:
;  SDSSPATH --  maybe already set in .cshrc file, but just in case- a
;               path to the data
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_balupdate, 'DR1_QSO/sdss_DR1_QAL.fits', $
;      'DR2_QSO/sdss_DR2_QAL.fits', ['DR2_QSO/matched.fits']
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Jun-2004 Written by SHF
;   23-Nov-2004 Update by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;

pro sdss_balupdate, oldqalfil, newqalfil, matchfil, SDSSPATH=sdsspath

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'sdss_balupdate, oldqal, newqal, [matchfil], SDSSPATH= [v1.1]'
    return
  endif 

  if not keyword_set( SDSSPATH ) then sdsspath = 'SDSSPATH'

  x=xmrdfits(getenv('SDSSPATH')+oldqalfil, 1, /silent)
  y=xmrdfits(getenv('SDSSPATH')+newqalfil, 1, /silent)
  
  n2=n_elements(y)
  
  samearr=lonarr(n2)
  
  for i=0L, n2-1 do begin
      
      isame=where(abs(y[i].ra - x.ra) LT 0.00001 AND $
                  abs(y[i].dec - x.dec) LT 0.0001 AND $
                  abs(y[i].z_qso - x.z_qso) LT .01, nisame)
      
      case nisame of 
          0:                    ;print, 'no match! NEW to DR2', i
          else: begin
;        printcol, nisame, x[isame].flg_bal, i
              y[i].flg_bal = x[isame[0]].flg_bal
              samearr[i] = i
          end
      endcase
      if i MOD 1000L EQ 0 then print, 'i = ', i
  endfor
  
;  samearr = samearr[where(samearr GT 0)]
  
  mwrfits, y, getenv('SDSSPATH')+newqalfil, /create
  print, 'All done!'
  
;  mwrfits, samearr, getenv('SDSSPATH')+matchfil, /create
  return

end
