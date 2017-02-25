;+ 
; NAME:
; cldy_read_cuba
;   Version 1.1
;
; PURPOSE:
;  Read a CUBA file and pass back wave, Jnu
;
; CALLING SEQUENCE:
;
; INPUTS:
;   fil  - CUBA output file
;   z    - Redshift
;
; RETURNS:
; OUTPUTS:
;;  WAVE = wavelength
;  JNU  = radiation field
; OPTIONAL INPUTS:
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; file, '/u/xavier/Cloudy/Spec/Data/CUBA/Q1G0/bkgthick.out',
; 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;  26-May-2010 Written by JFH
;  15-Nov-2016 updated for H&M 2012 CUBA, KLC
;-
;------------------------------------------------------------------------------

PRO CLDY_READ_CUBA, file, z, WAVE = WAVE_OUT, JNU = JNU_OUT, hm2012=hm2012

if  N_params() LT 2  then begin 
   print, 'Syntax -  read_cuba, file, z, [WAVE=,JNU=, /hm2012]'
   return
endif 

COMMON CUBA_SPEC, zval, wave, jnu, file_current, nwave, nblock, nz

;; This is redundant based on N_params()
;IF NOT KEYWORD_SET(FILE) OR NOT KEYWORD_SET(Z) THEN $
;   message, 'Problem with inputs'

IF NOT KEYWORD_SET(FILE_CURRENT) THEN FILE_CURRENT = ''

IF strmatch(file, file_current) EQ 0 THEN BEGIN
    file_current =  file
;   compute the number of total lines in the Madau file ignoring the 
;   comment lines
    readcol, file, junk, comment = '#', FORMAT = 'F'
    
    if keyword_set(hm2012) then begin
       ;; but if H&M2012
       ;; (http://www.ucolick.org/~pmadau/CUBA/DOWNLOADS.html stripped
       ;; of header and format like $XIDL_DIR...cuba*.out): 1 block of
       ;; 575
       nwave = 575L
       nz = 60L
       fmt = '(12x,60f12.5)'
    endif else begin
       ;; For e.g., $XIDL_DIR/Cloudy/cuba_q1g01_bkgthick.out, data should
       ;; be in blocks of 432
       nwave = 432L
       nz = 10L
       fmt = '(11x,60f11.5)'
    endelse
    nblock = n_elements(junk)/nwave
    
;; Open Madau file
    close, /all
    openr, 1, file
    zval = dblarr(nz*nblock)
    zval_temp = fltarr(nz)
    wave = dblarr(nwave)
    jnu = dblarr(nwave, nz*nblock)
    dumf = dblarr(nz+1)
    FOR ii = 0L, nblock-1L DO BEGIN
        readf, 1, zval_temp, FORMAT = fmt 
        zinds = ii*nz + lindgen(nz)
        zval[zinds] = zval_temp
        FOR jj = 0L, nwave-1L DO BEGIN
            readf, 1, dumf          
            jnu[jj, zinds] = dumf[1:*]
            IF ii EQ 0 THEN wave[jj] = dumf[0]
        ENDFOR
    ENDFOR
close, 1
ENDIF

jnu_out = dblarr(nwave)
wave_out = wave

; Spline interpolate onto redshift of interest for each wavelength
FOR ii = 0, nwave-1L DO $
   jnu_out[ii] = interpol(jnu[ii, *], zval, z, /spline)


RETURN
END
