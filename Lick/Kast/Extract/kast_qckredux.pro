;+ 
; NAME:
; kast_qckredux   
;   Version 1.1
;
; PURPOSE:
;   Perform a quick and dirty (and rather accurate) reduction of a
;  Kast spectrum.  Very useful at the telescope.
;
; CALLING SEQUENCE:
;  kast_qckredux, img, flat, spec, sig, ARC=, OVR=, $
;                  CLINE=, APER=, SKYFUNC=, SKYNORD=, SKYREG=, 
;                  /DISPLAY, /NOREJ, TRACE=, /NOSKY=, STRCT=, RN=, GAIN=, 
;                  VAR=, CRUDE=, TINTER=, /SILENT,
;                  WAVE=, /OPTIMAL, /AINTER, OBJLIN=, /DEBUG, /APINTER, 
;                  YMODEL=, SIGMA=, /CHK=
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
;  /CRUDE   -- Use the routine trace_crude to trace [recommended]
;  CLINE=   -- Row to use to identify object  [default: half way on image]
;  APER=    -- Aperture for boxcar extraction  (output or input)
;  TRACE=   -- Trace of the object (output or input)
;  /CHKTRC  -- Check the trace by eye
;  /DISPLAY -- Display the sky subtracted image with xatv
;  /NOREJ   -- Turn off rejection in sky subtraction (not recommended)
;  /AINTER  -- Interactively fiddle with the Arc peakup
;  OBJLIN=  -- Row near the object trace [default: center of image]
;  /OPTIMAL -- Perform optimal extraction
;  YMODEL   -- 2D model of the image created by extract_image
;  RN=      -- Readnoise of the image [default: 3.29]
;  ARC=     -- Arc image (fits file)
;  /DEBUG   -- Debug the program (this one often 'breaks')
;  SIGMA=   -- Approximate sigma of the Gaussian profile (for Optimal)
;
; OPTIONAL OUTPUTS:
;  WAVE=    -- Wavelength array (requires Arc image)
;
; COMMENTS:
;
; EXAMPLES:
;   kast_qckredux, 'b109.ccd', spec, sig, WAVE=wave, ARC='b102.ccd',
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

pro kast_qckredux, img, flat, spec, sig, ARC=arc, OVR=ovr, $
                   CLINE=cline, APER=aper, GRISM=grism, $
                   SKYFUNC=skyfunc, SKYNORD=skynord, SKYREG=skyreg, $
                   DISPLAY=DISPLAY, NOREJ=norej, TRACE=trace, NOOV=noov, $
                   NOSKY=nosky, STRCT=strct, RN=rn, GAIN=gain, VAR=var, $
                   CRUDE=crude, TINTER=tinter, SILENT=silent,$
                   WAVE=wave, OPTIMAL=optimal, AINTER=ainter, $
                   OBJLIN=objlin, DEBUG=debug, APINTER=apinter, $
                   YMODEL=ymodel, SIGMA=sigma, CHK=chk


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'kast_qckredux, img, spec, sig, CLINE=, APER=, SKYFUNC=, SKYNORD=, '
    print, '        SKYREG=, /DISPLAY, /NOREJ, TRACE=, /NOOV, /NOSKY, STRCT=,'
    print, '        RN=, GAIN=, VAR=, /crude, /tinter, ERROR=, /SILENT,'
    print, '        /CHK, SIGMA=, [v1.1]'
    return
  endif 

;  Optional Keywords

  ; gain and read noise
  if not keyword_set(GAIN) then gain=3.9
  if not keyword_set(GRISM) then grism = 'B1'
  if not keyword_set(RN) then rn=6.5 ;blue-- red = 8.3 for slow

; Test if it is a file or image

  ovdat = x_readimg(img, /fscale)
  ovdat = ovdat * gain
  sz = size(ovdat, /dimensions)
  if n_elements(ovdat) EQ 1 then return

; OBJLIN
  if not keyword_set(objlin) then objlin = sz[1]/2.

; FLAT
  if keyword_set( FLAT ) then begin
      print, 'kast_qckredux: Flattening...'
      flat = x_readimg(flat, /fscale)
      ;; Normalize
      medflt = djs_median(flat[*,long(objlin)-20L:long(objlin)+20L],2)
      bset = bspline_iterfit(findgen(n_elements(medflt)),medflt, yfit=yfit,$
                                     everyn=15)
      nrmflat = flat / ( yfit # replicate(1., sz[1]))
      ovdat = ovdat / nrmflat
  endif
      
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
      if not keyword_set( SILENT) then $
        print, 'kast_qckredux: Define the Aperture'
      apspec = djs_median(ovdat[cline-15L:cline+15L, *],1)
      x_fndpeaks, apspec, center, nsig=3., /force
      mn = min(abs(center-objlin),imn)
      objlin = center[imn]
      if not keyword_set( APINTER ) then $
        aper = x_setaper(apspec, objlin, 0.05, RADIUS=20L) else begin
          aper = x_setapergui(apspec)
          aper[0] = objlin - aper[0]
          aper[1] = aper[1] - objlin 
      endelse
      skyaper = aper+3.
      if keyword_set(DEBUG) then stop
      print, 'kast_qckredux: Aper = ', aper
  endif

; TRACE
  if not keyword_set( TRACE ) then begin
      if not keyword_set( SILENT) then print, 'x_apall: Tracing'
      trace = x_trace(ovdat, objlin, $
                      cline, VAR=var, TFITSTR=tfitstr, $
                      /CRUDE, /ROT, INTER=tinter)
      ;; Test
      if keyword_set( DEBUG ) OR keyword_set( CHK ) then begin
          rnd_trc2 = lindgen(sz[0])
          trc_msk = rnd_trc2 + round(trace)*sz[0]
          tmp = ovdat
          tmp[trc_msk] = -10000
          xatv, tmp, /block
      endif
  endif

;;;;; SKY ;;;;;;;;;;
  if not keyword_set( NOSKY ) AND not keyword_set( OPTIMAL ) then begin
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

      ;; SUBTRACT SKY
      if not keyword_set( SILENT) then print, 'x_apall: Subtracting the Sky'
      skyreg = [ [-30., skyaper[1]], [-skyaper[0], 30.]] 
      skyreg[*] = skyreg[*] + objlin
      fimg = x_skysub(ovdat, sky, REG=skyreg, CREG=cline, $
                      TRACE=trace, NOREJ=norej, SKYFSTR=skyfstr, SKYRMS=skyrms) 
  endif else fimg = ovdat ; NO SKY SUBTRACTION
                    
;;;;;;; DISPLAY
  if keyword_set( DISPLAY ) or keyword_set( CHK ) $
    then xatv, fimg, /block, min=-10, max=100


; EXTRACT

  if not keyword_set( OPTIMAL ) then begin
      if not keyword_set( SILENT) then print, 'x_apall: Extracting'
      if arg_present( SIG ) then begin
          spec = x_extract(fimg, [objlin-aper[0],objlin+aper[1]], $
                           trace, spec_var, sky, CAPER=cline, $
                           SKYRMS=skyrms, RN=rn, GAIN=gain, VAR=var) 
          fsig = sqrt(spec_var)
      endif else spec = x_extract(fimg, [objlin-aper[0],objlin+aper[1]], $
                                  trace, CAPER=cline)
  endif else begin
      cutimg = transpose(fimg[*, round(objlin)-20L:round(objlin)+20L])
      cutvar = transpose(var[*, round(objlin)-20L:round(objlin)+20L])
      ivar = 1. / (cutvar > 1.)
      ;; Kludge the trace
      xcen = trace-objlin+20
      kludge = [[xcen], [xcen*0.+9999.]]
      if not keyword_set( SIGMA ) then sigma = (total(aper)/4.) > 1.0
      ;; Optimal
      extract_image, cutimg, ivar, kludge, wfixed=[1,1,1], $
        sigma, flux, finv, ymodel=ymodel, nPoly=2
      spec = flux[*,0]
      fsig = sqrt(1./finv[*,0])
      if keyword_set( DEBUG ) then stop
      if keyword_set( CHK ) then xatv, ymodel, /block
  endelse
                     
; ARC
  if keyword_set( ARC ) then begin
      if not keyword_set( SILENT) then print, 'x_apall: Wavelength calib'
      ;; Assume B1
      case grism of 
          'B1': begin
              templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_g1.idl'
              linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_blue3.lst'
          end
          'B3': begin
              ;templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_g3.idl'
              templt = 'kastfit_g3.idl'
              linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_blue2.lst'
          end 
          'R4': begin
              ;templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_12005000.idl'
              templt = 'kastfit_12005000.idl'
              linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_blue2.lst'
          end 
          else: stop
      endcase
      ;; Find Shift
      restore, templt
      ;; Normalize
      a = where(aspec GT 2000., na)
      if na NE 0 then aspec[a] = 2000.
      ;; Read arc
      arc = x_readimg(arc, /fscale)
      ;; Extract
      arcsp = x_extract(arc, [objlin-7.,objlin+7.], $
                       trace, CAPER=cline)
      ;; Normalize
      tmp = arcsp
      a = where(tmp GT 2000., na)
      if na NE 0 then tmp[a] = 2000.
      ;; Cross-correlate
      step = lindgen(400L) - 200L
      corr = c_correlate(aspec, tmp, step, /double)
      mx = max(corr, imx)
      imx = step[imx]
      print, 'kast_qckredux: Offseting ', strtrim(imx,2), ' pixels'

      ;; Line list
      x_arclist, linlist, lines
      ;; ID
      x_templarc, arcsp, lines, calib, SHFT=imx
      ;; Peak up
      x_arcpeakup, arcsp, lines, FFIT=ffit, WV=wave, INTER=ainter
  endif


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
  if arg_present(SIG) then sig = fsig
  return
end
