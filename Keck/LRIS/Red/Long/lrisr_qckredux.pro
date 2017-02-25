;+ 
; NAME:
; lrisr_qckredux   
;   Version 1.1
;
; PURPOSE:
;    Perform a quick (and reasonably final) reduction of long
;  slit spectrum taken with LRIS-B
;
; CALLING SEQUENCE:
;  lrisb_qckredux, img, flat, spec, sig, CLINE=, APER=,
;    SKYFUNC=, SKYNORD=, SKYREG=, /DISPLAY, /NOREJ, TRACE=, 
;    /NOSKY, STRCT=, RN=, VAR=, CRUDE=, TINTER=, /SILENT,
;    ARC=, WAVE=, /OPTIMAL, AINTER=, OBJLIN=, /DEBUG, APINTER=,
;    YMODEL=, STD=, GRISM=, CHKTRC=, /FULL
;
; INPUTS:
;   img  -- Image to extract from (raw fits file)
;  flat  -- Name of flat file (fits)
;
; RETURNS:
;  spec  -- 1D spectrum
;   sig  -- Error of 1D spectrum
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CLINE=   -- Row to use to identify object  [default: 1300L]
;  APER=    -- Aperture for boxcar extraction  (output or input)
;  TRACE=   -- Trace of the object (output or input)
;  /CHKTRC  -- Check the trace by eye
;  /DISPLAY -- Display the sky subtracted image with xatv
;  /NOREJ   -- Turn off rejection in sky subtraction (not recommended)
;  /FULL    -- Image was readout in all amps (i.e. long slit mode not enabled)
;  /AINTER  -- Interactively fiddle with the Arc peakup
;  OBJLIN=  -- Row near the object trace [default: 150 or 104] 
;  /OPTIMAL -- Perform optimal extraction
;  YMODEL   -- 2D model of the image created by extract_image
;  RN=      -- Readnoise of the image [default: 3.29]
;  ARC=     -- Arc image (fits file)
;  /DEBUG   -- Debug the program (this one often 'breaks')
;  GRISM    -- Name of GRISM used (e.g. '400')
;
; OPTIONAL OUTPUTS:
;  WAVE=    -- Wavelength array (requires Arc image)
;
; COMMENTS:
;
; EXAMPLES:
;   lris_qckredux, 'b109.ccd', spec, sig, WAVE=wave, ARC='b102.ccd',
;                  FLAT='b100.ccd', /OPTIMAL, /NOSKY 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lrisr_qckredux, img, flat, spec, sig, CLINE=cline, APER=aper, $
                    SKYFUNC=skyfunc, SKYNORD=skynord, SKYREG=skyreg, $
                    DISPLAY=DISPLAY, NOREJ=norej, TRACE=trace, $
                    NOSKY=nosky, STRCT=strct, RN=rn, VAR=var, $
                    CRUDE=crude, TINTER=tinter, SILENT=silent,$
                    ARC=arc, WAVE=wave, OPTIMAL=optimal, AINTER=ainter, $
                    OBJLIN=objlin, DEBUG=debug, APINTER=apinter, $
                    YMODEL=ymodel, STD=std, GRISM=grism, CHKTRC=chktrc, $
                    FULL=full, SAPER=saper, NOFLAT=noflat, NSKY=nsky, $
                    FIMG=fimg


;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'lrisr_qckredux, img, flat, spec, sig, CLINE=, APER=, ' + $
      'SKYFUNC=, SKYNORD=, '
    print, '        SKYREG=, /DISPLAY, /NOREJ, TRACE=, /NOSKY, STRCT=,'
    print, '        VAR=, /crude, /tinter, /SILENT, SIG=, /CHKTRC ) [V1.1]'
    return
  endif 

;  Optional Keywords

  if not keyword_set(SAPER) then saper = 10.
  if not keyword_set(GRISM) then grism = '400'
  if not keyword_set(RN) then rn = 3.29
  if not keyword_set(NSKY) then nsky = [90., 90.]
  if not keyword_set(FULL) then icut = 200L ELSE icut = 275L

; Process
  print, 'lrisb_qckredux: Processing the image'
  lrisr_proc, img, flat, finimg, /LONG, FULL=full, NOFLAT=noflat
  sz = size(finimg, /dimensions)

; OBJLIN
  if not keyword_set(objlin) then begin
      if keyword_set( FULL ) then objlin = 430. else objlin = 590
  endif

  ; variance
  if not keyword_set(VAR) then begin
      var = fltarr(sz[0], sz[1])
      var = (finimg + rn^2 + (finimg EQ 0)) > (rn^2)
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
      if not keyword_set( CLINE ) then  cline = 1024L
  endelse

; Set Aperture

  if not keyword_set( APER ) then begin
      if not keyword_set( SILENT) then $
        print, 'lrisb_qckredux: Define the Aperture'
      apspec = djs_median(finimg[*,cline-15L:cline+15L],2)
      x_fndpeaks, apspec, center, nsig=3., /force
      mn = min(abs(center-objlin),imn)
      objlin = center[imn]
      ;; AUTO
      if not keyword_set( APINTER ) then $
        aper = x_setaper(apspec, objlin, 0.10, RADIUS=20L) $
      else begin ; Interactive
          blah = x_setapergui(apspec, OBJSTR=sobjstr)
          objlin = sobjstr.ycen
          aper = sobjstr.aper
      endelse
      skyaper = aper+saper
      if keyword_set( STD ) then begin
          aper = aper + 15
          skyaper = skyaper + 100
      endif
      print, 'lrisb_qckredux: Aperture = ', aper
      if keyword_set(DEBUG) then stop
  endif else skyaper = aper+saper

; TRACE
  if not keyword_set( TRACE ) then begin
      if not keyword_set( SILENT) then print, 'lrisb_qckredux: Tracing'
      trace = x_trace(finimg, objlin, $
                      cline, VAR=var, TFITSTR=tfitstr, $
                      /CRUDE, INTER=tinter)
      ;; Test
      if keyword_set( DEBUG ) or keyword_set( CHKTRC ) then begin
          rnd_trc2 = lindgen(sz[1])
          trc_msk = round(trace) + rnd_trc2*sz[0]
          tmp = finimg
          tmp[trc_msk] = -10000
          xatv, tmp, /block
      endif
      ;; Debug
      if keyword_set(DEBUG) then stop
  endif

;;;;; SKY ;;;;;;;;;;
  if not keyword_set( NOSKY ) then begin
      if not keyword_set( SKYFSTR ) then begin
          skyfstr = { fitstrct }
          if not keyword_set( SKYFUNC ) then skyfstr.func = 'POLY'
          if not keyword_set( SKYNORD ) then skyfstr.nord = 2
          if not keyword_set( SKYLOW ) then skyfstr.lsig = 2.
          if not keyword_set( SKYHIGH ) then skyfstr.hsig = 2.
          skyfstr.niter = 2
          skyfstr.flg_rej = 1
          skyfstr.maxrej = 5
      endif
      if keyword_set(STD) then skyfstr.nord = 0L

      ;; SUBTRACT SKY
      if not keyword_set( SILENT) then $
        print, 'lrisb_qckredux: Subtracting the Sky'

;      skyreg = [ [30, 120.], [90., 160.] ] $
      skyreg = [ [-nsky[0], skyaper[1]], [-skyaper[0], nsky[1]]] 
      skyreg[*] = skyreg[*] + objlin
          
      fimg = x_skysub(transpose(finimg), sky, REG=skyreg, CREG=cline, $
                      TRACE=trace, NOREJ=norej, SKYFSTR=skyfstr, SKYRMS=skyrms) 
      fimg = transpose(fimg)
  endif else fimg = finimg ; NO SKY SUBTRACTION
                    
;;;;;;; DISPLAY
  if keyword_set( DISPLAY ) or keyword_set( DEBUG ) then $
    xatv, fimg, /block, min=-10, max=100

  if arg_present(YMODEL) and not keyword_set(OPTIMAL) then ymodel=fimg

  if keyword_set( DEBUG ) then stop


; EXTRACT

  if not keyword_set( OPTIMAL ) then begin
      if not keyword_set( SILENT) $
        then print, 'lrisb_qckredux: Extracting'
      if arg_present( SIG ) then begin
          spec = x_extract(transpose(fimg), $
                           [objlin-aper[0],objlin+aper[1]], $
                           trace, fvar, sky, CAPER=cline, $
                           SKYRMS=skyrms, RN=rn, GAIN=1., VAR=transpose(var)) 
          fsig = sqrt(fvar)
      endif else spec = x_extract(transpose(fimg), $
                                  [objlin-aper[0],objlin+aper[1]], $
                                  trace, CAPER=cline)
  endif else begin
      cutimg = fimg[0:icut,*]
      cutvar = var[0:icut,*]
      ivar = 1. / (cutvar > 1.)
      ;; Kludge the trace
      xcen = trace
      kludge = [[xcen], [xcen*0.+9999.]]
      sigma = (total(aper)/4.) > 1.5
      ;; Optimal
      extract_image, cutimg, ivar, kludge, wfixed=[1,1,1], $
        sigma, flux, finv, ymodel=ymodel, nPoly=2
      print, 'lrisb_qckredux:  Optimal Sigma = ', sigma
      spec = flux[*,0]
      fsig = sqrt(1./finv[*,0])
      if keyword_set( DEBUG ) then stop
  endelse
                     
; ARC
  if keyword_set( ARC ) then begin
      if not keyword_set( SILENT) then $
        print, 'lrisb_qckredux: Wavelength calib'
      ;; Assume B1
      linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/lris_blue.lst'
      case strtrim(grism,2) of 
          '400': begin
              templt = $
                getenv('XIDL_DIR')+'/Keck/LRIS/CALIBS/lrisbfit_400.idl'
          end
          '1200': begin
              templt = $
                getenv('XIDL_DIR')+'/Keck/LRIS/CALIBS/lrisbfit_1200.idl'
          end
          else: stop
      endcase
      ;; Find Shift
      restore, templt
      ;; Normalize
      a = where(aspec GT 2000., na)
      if na NE 0 then aspec[a] = 2000.
      ;; Read arc
;      aimg = xmrdfits(arc, /fscale, /silent)
      lrisb_proc, arc, flat, aimg, /LONG, FULL=full
      ;; Extract
      arcsp = x_extract(transpose(aimg), [objlin-7.,objlin+7.], $
                       trace, CAPER=cline)
      ;; Normalize
      tmp = arcsp
      a = where(tmp GT 2000., na)
      if na NE 0 then tmp[a] = 2000.
      ;; Cross-correlate
      step = lindgen(400L) - 200L
      corr = c_correlate(aspec, tmp, step, /double)
      mx = max(corr, imx)
;      stop
      imx = step[imx]
      print, 'lris_qckredux: Offseting ', strtrim(imx,2), ' pixels'

      ;; Line list
      x_arclist, linlist, lines
      ;; ID
      x_templarc, arcsp, lines, calib, SHFT=imx
      ;; Peak up
      x_arcpeakup, arcsp, lines, FFIT=ffit, WV=wave, INTER=ainter, $
        SIG=[2.,2.], NFRST=5L
  endif


; Release memory
;  delvarx, finimg, dat, fimg, var
  if keyword_set( SKY ) then delvarx, sky

; Fill in the structure
  if arg_present( STRCT ) then begin
      strct = { extrctstrct }
;      if keyword_set(ovr) then strct.ovr = ovr ; OV
      strct.aper = aper
      strct.cline = cline
      if not keyword_set( NOSKY ) then begin ; SKY
          if keyword_set( skyreg ) then strct.skyreg = ptr_new(skyreg)
          strct.skyfstr = skyfstr
      endif
      if keyword_set( trace ) then strct.trace = ptr_new(trace)
  endif

  if not keyword_set( SILENT) then print, 'x_apall: All Done!'
  if arg_present(SIG) then sig = fsig
  return
end
