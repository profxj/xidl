;+ 
; NAME:
; x_apall   
;   Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   spec = x_apall(ydat, [head])
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   wave       - wavelength array
;   DISPLAY    - Display the sky subtracted image with xatv
;   OVR        - String array for ov region:  '[2050:2100, *]'
;   ERROR      - Variance array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   spec = x_apall('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Nov-2001 Written by JXP
;   05-Dec-2001 Modularized most of the routine
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_apall, img, OVR=ovr, CLINE=cline, APER=aper, $
                  SKYFUNC=skyfunc, SKYNORD=skynord, SKYREG=skyreg, $
                  DISPLAY=DISPLAY, NOREJ=norej, TRACE=trace, NOOV=noov, $
                  NOSKY=nosky, STRCT=strct, RN=rn, GAIN=gain, VAR=var, $
                  CRUDE=crude, TINTER=tinter, ERROR=error, SILENT=silent


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'spec = x_apall(img, OVR=, CLINE=, APER=, SKYFUNC=, SKYNORD=, '
    print, '        SKYREG=, /DISPLAY, /NOREJ, TRACE=, /NOOV, /NOSKY, STRCT=,'
    print, '        RN=, GAIN=, VAR=, /crude, /tinter, ERROR=, /SILENT ) [V1.1]'
    return, -1
  endif 

;  Optional Keywords

  ; gain and read noise
  if not keyword_set(GAIN) then gain=1.
  if not keyword_set(RN) then rn=5.

; Test if it is a file or image

  dat = x_readimg(img, /fscale)
  if n_elements(dat) EQ 1 then return, -1
      
  ; Overscan
  if keyword_set( strct ) then ovr = strct.ovr

  if not keyword_set( NOOV ) then $
    xoverscan, dat, 'GEN', OVIMG=ovdat, OVSEC=ovr, /nofits $
  else ovdat = dat              ; No OV subtraction
      
  sz = size(ovdat, /dimensions)

  ; variance
  if not keyword_set(VAR) then begin
      var = fltarr(sz[0], sz[1])
      var = (ovdat + rn^2 + (ovdat EQ 0)) > (rn^2)
  endif

  ; EXTRACTION STRUCTURE
  if keyword_set( strct ) then begin
      cline = strct.cline
      aper = strct.aper
      trace = *strct.trace
      skyfstr = strct.skyfstr
      skyreg = *strct.skyreg
  endif else begin
      ; CLINE
      if not keyword_set( CLINE ) then  cline = sz[0]/2L
  endelse

; Set Aperture

  if not keyword_set( APER ) then begin
      if not keyword_set( SILENT) then print, 'x_apall: Define the Aperture'
      aper = x_setapergui(ovdat, CLINE=cline)
  endif

; TRACE
  if not keyword_set( TRACE ) then begin
      if not keyword_set( SILENT) then print, 'x_apall: Tracing'
      trace = x_trace(ovdat, aper, cline, VAR=var, TFITSTR=tfitstr, $
                      CRUDE=crude, /ROT, INTER=tinter)
  endif

;;;;; SKY ;;;;;;;;;;
  if not keyword_set( NOSKY ) then begin
      if not keyword_set( SKYFSTR ) then begin
          skyfstr = { fitstrct }
          if not keyword_set( SKYFUNC ) then skyfstr.func = 'POLY'
          if not keyword_set( SKYNORD ) then skyfstr.nord = 1
          if not keyword_set( SKYLOW ) then skyfstr.lsig = 2.
          if not keyword_set( SKYHIGH ) then skyfstr.hsig = 2.
          skyfstr.niter = 2
      endif

      ; FIT Interactively to set the Regions
      if not keyword_set( SKYREG ) then begin
          if not keyword_set( SILENT) then print, 'x_apall: Define sky regions'
          fit = x1dfit(findgen(sz[1]), $
                       djs_median(ovdat[cline-10L:cline+10L,*],1), $
                       FITSTR=skyfstr, $
                       REG=skyreg, /INTER)
      endif
      
      ; SUBTRACT SKY
      if not keyword_set( SILENT) then print, 'x_apall: Subtracting the Sky'
      fimg = x_skysub(ovdat, sky, REG=skyreg, CREG=cline, $
                      TRACE=trace, NOREJ=norej, SKYFSTR=skyfstr, SKYRMS=skyrms) 
  endif else fimg = ovdat ; NO SKY SUBTRACTION
                    
;;;;;;; DISPLAY

  if keyword_set( DISPLAY ) then xatv, fimg, /block

; EXTRACT

  stop
  if not keyword_set( SILENT) then print, 'x_apall: Extracting'
  if arg_present( ERROR ) then begin
      spec = x_extract(fimg, aper, trace, error, sky, CAPER=cline, $
                       SKYRMS=skyrms, RN=rn, GAIN=gain) 
  endif else spec = x_extract(fimg, aper, trace, CAPER=cline)
                     

; Release memory
  delvarx, ovdat, dat, fimg, var
  if keyword_set( SKY ) then delvarx, sky

; Fill in the structure
  if arg_present( STRCT ) then begin
      strct = { extrctstrct }
      if keyword_set(ovr) then strct.ovr = ovr ; OV
      strct.aper = aper
      strct.cline = cline
      if not keyword_set( NOSKY ) then begin ; SKY
          if keyword_set( skyreg ) then strct.skyreg = ptr_new(skyreg)
          strct.skyfstr = skyfstr
      endif
      if keyword_set( trace ) then strct.trace = ptr_new(trace)
  endif

  if not keyword_set( SILENT) then print, 'x_apall: All Done!'
  return, spec
end
