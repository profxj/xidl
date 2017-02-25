;+ 
; NAME:
; sdss_fix_rep
;    Version 1.1
;
; PURPOSE:
;   Remove repeat entries from the 'checked' file of QAL.  Defaulted
;  to deal with DR2 using dr2_bd.list
;
; CALLING SEQUENCE:
;  sdss_fix_rep, fil=, bdlst=, new=new
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  FIL= -- File name of file including duplicates
;  BDLST= -- List of indices which are to be removed
;  NEW=   -- Name of new DR2_chkd file [default: 'DR2_chkd.fits']
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Apr-2004 Written by SHF
;-
;------------------------------------------------------------------------------
pro sdss_fix_rep, fil=fil, bdlst=bdlst, new=new

if not keyword_set(fil) then $
  fil = getenv('SDSSPATH')+'/DR2_QSO/DR2_chkd.fits'

;read in chkd file
dlastr = xmrdfits(fil, 1, /silent)
tot = n_elements(dlastr)

;read in bad_i list
if not keyword_set(bdlst) then $
  bdlst = getenv('SDSSPATH')+'/DR2_QSO/dr2_bd.lst'     
readcol, bdlst, bd_i, format = 'I'

;mask all, delete bdlst
 msk = bytarr(tot)
 msk[*] = 1B
 msk[bd_i] = 0B
 gd_i = where(msk EQ 1B)

;new strct
new_dlastr = dlastr[gd_i]

;create new strct
if not keyword_set(new) then new = getenv('SDSSPATH')+'/DR2_QSO/DR2_chkd.fits'
mwrfits, new_dlastr, new, /create, /silent

;stop
end
