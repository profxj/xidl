;+ 
; NAME:
; esi_echtrcstd   
;     Version 1.1
;
; PURPOSE:
;    Trace a standard star in each ordrer
;
; CALLING SEQUENCE:
;   
;  esi_echtrcstd, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  Image with trcstdered light removed
;
; OPTIONAL KEYWORDS:
;  /NOFND   - Do not repeat step to find object
;  /NOSKY   - Do not repeat sky subtraction
;  /CHK     - Show the final trace
;  MXERR=   - Maximum error in trace centering to include 
;                 (default: 0.3 pix)
;  OFF=     - Offset used in find object routine
;  /CUAR    - Data has only CuAr lamps
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Only setup for 1x1 binning (some hard numbers)
;
; EXAMPLES:
;   esi_echtrcstd, esi, 1.0, /CHK
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Aug-2002 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echtrcstd, esi, slit, NOFND=nofnd, NOSKY=nosky, CHK=chk, MXERR=mxerr $
                   , OFF = off, CUAR = cuar, CBIN = cbin, RBIN = rbin $
                   , FWHM = FWHM, CLOBBER = CLOBBER, BIASFIL=biasfil, $
                   FLATFIL=flatfil, SEDG_FIL=sedg_fil, NOOSCAN=nooscan

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echtrcstd, esi, slit, /NOFND, /NOSKY, /CHK, OFF=, MXERR= '
      print, '     /CUAR, /NOOSCAN [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set( NCOLL ) then ncoll = 10L
  if not keyword_set( MXERR ) then mxerr = 0.3

;;;;;;
;  Find standard star
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.slit EQ slit AND esi.type EQ 'STD', nindx)
  indx2 = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
              esi.type EQ 'STD', nindx2)
  ;; if no standard matching the slit is found, relax and search
  ;; for standard with non-matching slit. If it exists, use. 
  IF nindx EQ 0 AND nindx2 GT 0 THEN BEGIN
     print,'esi_echtrcstrd: No standard star images for current slit'
     print,'esi_echtrcstrd: standard exists for slit=',esi[indx2[0]].slit
     print,'esi_echtrcstrd: Hit continue provided you are okay with using a standard from a different slit'
     nindx=nindx2
     indx=indx2
     slit = esi[indx2[0]].slit
     stop
  ENDIF
  
  IF nindx EQ 0 and nindx2 EQ 0 THEN BEGIN
        print, 'esi_echtrcstrd: No standard star images! Returning..' 
        return
  ENDIF

  IF nindx EQ 1 THEN BEGIN
     print, 'esi_echtrcstd: Tracing standard star image --- ', esi[indx].img_root
  ENDIF

  IF nindx GT 1 THEN BEGIN
          print, 'esi_echtrcstd: Warning -- Multiple standard star images'
          indx = indx[0]
          print, 'esi_echtrcstd: Taking first one ', esi[indx].img_root
  ENDIF

;;;;;;;;;
; Process
  esi_echproc, esi, indx, SUBSCAT = subscat, CLOBBER = CLOBBER, BIASFIL=biasfil,$
               FLATFIL=flatfil, NOOSCAN=nooscan
  obj_id = esi[indx[0]].OBJ_ID

;;;;;;;;;
; FIND Obj + Create Obj Structure
;  if not keyword_set( NOFND ) then $
  esi_echfndobj, esi, indx, NFIND = 1, CBIN = CBIN, /STD, FWHM = FWHM, $
                 sedg_fil=sedg_fil
;;;;;;;;;
; Sky subtract
  if not keyword_set( NOSKY ) then begin
     if keyword_set(CUAR) then ordrs = [0L, 8L] else ordrs = lindgen(10)
     esi_echskysub, esi, indx, /STD, ORDR = ordrs, FCHK = chk, CBIN = cbin, $
                    SEDG_FIL=sedg_fil
  endif
  
;; Run fndobj again on sky-subtracted frame  
  esi_echfndobj, esi, indx, NFIND = 1, CBIN = CBIN, /STD, /SKYSUB, CHK = CHK $
                 , FWHM = FWHM, sedg_fil=sedg_fil

  ;; OBJ
  objfil = esi[indx].obj_fil
  if x_chkfil(objfil+'*') EQ 0 then begin
      print, 'esi_echtrcstd: Object file not found ', objfil
      return
  endif
  objstr = xmrdfits(objfil, 1, STRUCTYP = 'dblsobjstrct', /silent)

  ;; TRACE
  std_trc = esi_getfil('std_trc', SLIT = slit, cbin = cbin, rbin = rbin, /name)
  mwrfits, objstr[0:9].trace, std_trc, /create, /silent
  
  return
end
              
      
      
