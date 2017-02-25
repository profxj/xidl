;+ 
; NAME:
; x_stdcalibfil
;    Version 1.1
;
; PURPOSE:
;    Given an ra and dec, returns calib file and standard star name
;
; CALLING SEQUENCE:
;   x_stdcalibfil, ra, dec, calibfil, std_name
;
; INPUTS:
;   ra   - RA in : format (12:11:23.3)
;   dec  - DEC in : format (+13:32:21.2)
;
; RETURNS:
;   calibfil -- Name of calib file in $XIDL_DIR/Spec/calibs/standards
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
;   12-Jan-2008 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_stdcalibfil, ra, dec, calibfil, std_name, STD_RA=std_ra, STD_DEC=std_dec, $
                   TOLER=toler
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'x_stdcalibfil, ra, dec, calibfil, [std_name] TOLER= [v1.0]'
    return
  endif 

  calibfil = ''
  if not keyword_set(TOLER) then toler = 100.  ; 1 arcmin

  ;; Approved Standard files
  standards = [ $
              ['Feige34',  '10:39:36.74',   '+43:06:09.2',  'feige34_005.fits'], $
              ['Feige110',  '23:19:58.40',   '-05:09:56.2',  'feige110_005.fits'], $
              ['BD+28deg4211', '21:51:11.0',   '+28:51:50',  'bd_28d4211_005.fits'], $
              ['BD+33deg2642', '15:51:59.9',   '+32:56:54.3',  'bd_33d2642_004.fits'], $
;              ['HZ_44', '13:23:35.3',   '+36:07:59.5',
;              'hz44_stis_001.fits'], $
              ['HZ_44', '13:23:35.3',   '+36:07:59.5','hz44_stis_plus_50000k_ext.fits.gz'], $
              ['G191B2B',   '05:05:30.6119','+52:49:51.945','g191b2b_005.fits'] $
;              ['G191B2B',   '05:05:30.6119',
;              '+52:49:51.945','g191b2b_mod_007.fits'] $
;              ['G191B2B',   '05:05:30.6119', '+52:49:51.945','g191b2b_stisnic_003.fits'] $
              ]

  ;; Get RA/DEC
  x_radec, standards[1,*], standards[2,*], std_ra, std_dec
  x_radec, ra, dec, new_ra, new_dec

  ;; Find offset
  gcirc, 1, new_ra/15., new_dec, std_ra/15., std_dec, dis

  a = where(DIS LT TOLER, na)
  if NA EQ 0 then begin
      print, 'x_stdcalibfil: No standard within '+strtrim(round(toler),2)+ $
             ' arcsec.  ', min(dis)
      print, 'x_stdcalibfil: Better luck next time (contact JXP if you must).'
      return
  endif

  ;; Names
  calibfil = standards[3,a[0]]
  std_name = standards[0,a[0]]

  ;; RA and DEC
  std_ra = standards[1, a[0]]
  std_dec = standards[2, a[0]]

  return
end
      
