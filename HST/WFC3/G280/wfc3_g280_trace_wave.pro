;+ 
; NAME:
;  wfc3_g280_trace_wave
;
; PURPOSE:
;   This code uses the direct image centroid to generate a first guess
;   at the trace.  It then tweaks up this trace using Gaussian
;   centroiding along the columns.  Lastly, it generates a wavelength
;   array based on the AXE team coefficients.
;
; CALLING SEQUENCE:
;  wfc3_g280_trace_wave, wfc3_g280_strct, ii, specim, AXE1=, $
;                          JXP_KLUDGE=,CALIB=, OLD=, $
;                          EXT_BOX=, RANGE=, BEAM=, NOSHIFT=
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;   ii -- the index of the object in the structure
;   specim -- 2D spectral image
;
; RETURNS:
;
; OUTPUTS:
;   updates the structure with the trace and wavelength fits
;
; OPTIONAL KEYWORDS:
;   AXE1= -- Use the old AXE1 data (DEPRICATED)
;   JXP_KLUDGE= -- Offset in x for generating the wavelength solution
;                  [default: 3.0; based on eyeball comparison to SDSS
;                  spectra] (DEPRICATED)
;   CALIB= -- Supply your own calibration (NOT RECOMMENDED)
;   OLD= -- Use the old (circa 2010) code (DEPRICATED)
;   EXT_BOX= -- Extraction box [default:10L]
;   RANGE= -- Range to calculate the tracewidth [default: [2700d,4500]]
;   BEAM= -- Which beam to trace and find the wavelength for
;   NOSHIFT= -- Do not apply a y-pixel shift in fitting the trace 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  wfc3_g280_trace_wave, wfc3_g280_strct, ii, specim, AXE1=axe1, $
;                        JXP_KLUDGE=jxp_kludge,CALIB=calib, OLD=old, $
;                        EXT_BOX=ext_box, RANGE=range, BEAM=beam, $
;                        NOSHIFT=noshift
;
; PROCEDURES CALLED:
;  wfc3_g280_mkcalibstrct
;  gaussfit
;  poly_fit
;
; REVISION HISTORY:
;   23-Dec-2010 Written by JXP/JMO
;   10-Jun-2016 Major update in structure and how it runs by MN
;------------------------------------------------------------------------------
function wfc3_calc_coeff, calib, xcen,ycen, WAVE=wave, BEAM=beam

  case beam of
     0: begin
        if keyword_set(wave) then begin
           ncof=calib.DISP_ORDER_A
           cof=calib.DLDP_A
        endif else begin
           ncof=calib.DYDX_ORDER_A
           cof=calib.DYDX_A
        endelse
     end
     1: begin
        if keyword_set(wave) then begin
           ncof=calib.DISP_ORDER_C
           cof=calib.DLDP_C
        endif else begin
           ncof=calib.DYDX_ORDER_C
           cof=calib.DYDX_C
        endelse
     end
     else: stop
  endcase

  dim=size(cof,/dim)
  coeff=dblarr(ncof+1)
  for i=0L,ncof do begin
     if size(dim,/dim) eq 1 then coeff(i)=cof(i) else $
        case dim(0) of
        3: coeff(i)=cof(0,i)+cof(1,i)*xcen+cof(2,i)*ycen
        6: coeff(i)=cof(0,i)+cof(1,i)*xcen+cof(2,i)*ycen+$
                 cof(3,i)*xcen^2+cof(4,i)*xcen*ycen+$
                 cof(5,i)*ycen^2
        else: stop
     endcase
  endfor
  
  return, coeff 
end


pro wfc3_g280_trace_wave, wfc3_g280_strct, ii, specim, AXE1=axe1, $
                          JXP_KLUDGE=jxp_kludge,CALIB=calib, OLD=old, $
                          EXT_BOX=ext_box, RANGE=range, BEAM=beam, $
                          NOSHIFT=noshift

  if keyword_set(OLD) and size(JXP_KLUDGE,/type) EQ 0 then jxp_kludge = 3.0
  if not keyword_set(EXT_BOX) then ext_box = 10L
  if not keyword_set(RANGE) then range = [2700d,4500]
  
  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
          'trace_strct =  wfc3_g280_trace_wave(x0,y0,specim, JXP_KLUDGE=, WFC3=) [v1.0]'
    return
  endif 
  
  dim=size(specim,/dim)
  
  ;; Read in the calibration file into a structure
  if keyword_set(old) then chip = 3 else chip=wfc3_g280_strct(ii).chip
  calib = wfc3_g280_mkcalibstrct(chip, AXE1=axe1)

  case beam of
     0: begin
        beam_col=calib.BEAMA
        xoff=calib.XOFF_A
        yoff=calib.YOFF_A
        xoff2=calib.XOFF_A2
        yoff2=calib.XOFF_A2
     end
     1: begin
        beam_col=calib.BEAMC
        xoff=calib.XOFF_C
        yoff=calib.YOFF_C
        xoff2=calib.XOFF_C2
        yoff2=calib.XOFF_C2
     end
     else: stop
  endcase
        
  ;; Calculate the coefficients of the aXe Team
  trace_coeff=wfc3_calc_coeff(calib,wfc3_g280_strct(ii).x0+xoff,$
                              wfc3_g280_strct(ii).y0+yoff, BEAM=beam)
  wave_coeff=wfc3_calc_coeff(calib,wfc3_g280_strct(ii).x0+xoff,$
                             wfc3_g280_strct(ii).y0+yoff, BEAM=beam, /wave)

  ;; Theoretical Trace + wavelength + trace offset
  cnt = 0L
  trace_x=lonarr(beam_col[1]-beam_col[0]+1)
  trace_y_orig=dblarr(beam_col[1]-beam_col[0]+1)
  wave=dblarr(beam_col[1]-beam_col[0]+1)
  trace_offset=dblarr(beam_col[1]-beam_col[0]+1)
  sigma=dblarr(beam_col[1]-beam_col[0]+1)  
  for jj=beam_col[0],beam_col[1] do begin
     
     ;; Trace
     dx=float(jj)+round(wfc3_g280_strct(ii).x0)-wfc3_g280_strct(ii).x0
     dy=poly(dx,trace_coeff)
     trace_x(cnt)=wfc3_g280_strct(ii).x0 + dx + xoff + xoff2
     trace_y_orig(cnt)=wfc3_g280_strct(ii).y0 + dy + yoff + yoff2
     
     ;; Wavelength
     if keyword_set(old) then $
        dx=1.0*jj + xoff - xoff2 - JXP_KLUDGE ;; JXP kludge!
     wave(cnt)=poly(dx,wave_coeff)
     
     ;; Trace offset
     a=round(trace_y_orig(cnt))
     ;; take into account edges of chip
     if (a-ext_box ge 0 and a+ext_box lt dim(1)) and $
        (trace_x(cnt) ge 0 and trace_x(cnt) lt dim(0)) then begin
        ;;take ext_box window about trace y posn.
        xpix=indgen(2*ext_box+1)
        ypix=specim[trace_x(cnt),a-ext_box:a+ext_box]
        g=gaussfit(xpix,ypix,coeffs,nterms=4,sigma=coefferrs)
        trace_offset(cnt)=(a-ext_box+coeffs[1]) - trace_y_orig(cnt)
        sigma(cnt)=coeffs[2]
     endif else begin
        trace_offset(cnt)=0D
        sigma(cnt)=0D
     endelse 
     cnt = cnt+1
  endfor
  
  cnt=cnt-1

  ;; Tweak the Trace
  if not keyword_set(noshift) then begin
     ;; Choose wavelength range free of strong absorption
     trace_gwv=where(wave gt range[0] and wave lt range[1],ntgwv)  ;; Avoid LLS
     median_trace_offset=median(trace_offset[trace_gwv])
     trace_y_fit=trace_y_orig+median_trace_offset
     trace_y_fit=trace_y_orig+median_trace_offset
     print, 'wfc3_mktrace: Median offset = ', median_trace_offset
     ;; KLUDGE TO FIX TRACE C AT BLUE END
     a=lindgen(70)
     if beam eq 1 then begin
        b=a(where(abs(trace_offset(a)) lt 2d))
        c=poly_fit(trace_x(b),trace_offset(b),3,meas=sqrt(abs(trace_offset(a))))
        fit=poly(trace_x(a),c)
        trace_y_fit(a)=trace_y_orig(a)+fit
     endif
  endif else begin
     trace_y_fit=trace_y_orig
     median_trace_offset=0D
  endelse
  
  ;now look at the Sigma
  if keyword_set(old) then begin
     sigma_gwv=where(reverse(wave) ge 2700.0,nsgwv)
     spx=findgen(nsgwv)
     px2=findgen(cnt)
     pfy2=poly_fit(spx,sigma[sigma_gwv],3)
     polysigma=pfy2[0]+pfy2[1]*px2+pfy2[2]*px2*px2+pfy2[3]*px2*px2*px2
                                ;again, pull the average polynomial
                                ;from the calib file, but this time
                                ;scale multiplicatively to match the data
     sigmafit=calib.sigmaa_poly[0]+calib.sigmaa_poly[1]*px2+$
           calib.sigmaa_poly[2]*px2*px2+calib.sigmaa_poly[3]*px2*px2*px2
     gwv=where(wave ge 3200.0,ngwv)
     scale=median(sigma[gwv]/reverse(sigmafit[gwv]))
     sigmafit=reverse(sigmafit)*scale
     
  endif else begin
     ;; good wavelength away from LLS and no outlying sigma
     sigma_gd=where(wave ge range[0] and wave lt range[1] and $
                    sigma lt ext_box,nsgwv)
     
     ;; fit a line to the sigma of the width of the trace
     ;; (used to determine the sky region)
     pfy2=poly_fit(wave[sigma_gd],sigma[sigma_gd],0)
     ;; Save
     sigmafit=poly(wave,pfy2)
  endelse

  ;; Save in the correct BEAM:
  case beam of
     0: begin
        wfc3_g280_strct(ii).trace_xa=trace_x
        wfc3_g280_strct(ii).trace_ya_orig=trace_y_orig
        wfc3_g280_strct(ii).trace_ya_fit=trace_y_fit
        wfc3_g280_strct(ii).trace_offseta=trace_offset
        wfc3_g280_strct(ii).trace_sigmaa=sigma
        wfc3_g280_strct(ii).trace_sigma_fita=sigmafit
        wfc3_g280_strct(ii).median_trace_offseta=median_trace_offset
        wfc3_g280_strct(ii).wavea=wave
        wfc3_g280_strct(ii).cnta=cnt
     end
     1: begin
        wfc3_g280_strct(ii).trace_xc=trace_x
        wfc3_g280_strct(ii).trace_yc_orig=trace_y_orig
        wfc3_g280_strct(ii).trace_yc_fit=trace_y_fit
        wfc3_g280_strct(ii).trace_offsetc=trace_offset
        wfc3_g280_strct(ii).trace_sigmac=sigma
        wfc3_g280_strct(ii).trace_sigma_fitc=sigmafit
        wfc3_g280_strct(ii).median_trace_offsetc=median_trace_offset
        wfc3_g280_strct(ii).wavec=wave
        wfc3_g280_strct(ii).cntc=cnt
     end
     else: stop
  endcase
  
end
