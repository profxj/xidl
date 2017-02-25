;+ 
; NAME:
;  wfc3_g280_flux
;
; PURPOSE:
;   Fluxes the WFC3/G280 data using an archived sensitvity function.
;
; CALLING SEQUENCE:
;   
;  wfc3_g280_flux, wfc3_g280_strct, ii, calib, OLD=, BEAM=
;
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;   ii -- the index of the object in the structure
;
; RETURNS:
;
; OUTPUTS:
;   Updates the structure to contain the fluxed data
;
; OPTIONAL KEYWORDS:
;  OLD= -- Use the old (circa 2010) data
;  BEAM= -- Beam to reduce
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  wfc3_g280_flux, wfc3_g280_strct, ii, calib, OLD=old, BEAM=beam
;
; PROCEDURES CALLED:
;  bspline_valu
;
; REVISION HISTORY:
;   23-Dec-2010 Written by JXP/JMO
;   10-Jun-2016 Updated to deal with multibeam+structure by MN
;------------------------------------------------------------------------------

pro wfc3_g280_flux, wfc3_g280_strct, ii, calib, OLD=old, BEAM=beam

  if (N_params() LT 2) then begin 
     print,'Syntax - ' + $
           'wfc3_g280_flux, box_strct, exptime, strct, SENS_FIL= [v1.0]'
     return
  endif

  case beam of
     0: begin
        if keyword_set(OLD) then $
           sens_fil= getenv('XIDL_DIR')+'HST/WFC3/G280/gd71_sensfunc_dec2010.fits' $
        else $
           sens_fil=calib.sens_fila
        wave=wfc3_g280_strct(ii).wavea
        count=wfc3_g280_strct(ii).countsa
        var=wfc3_g280_strct(ii).vara
        sky=wfc3_g280_strct(ii).skya
     end
     1: begin
        sens_fil=calib.sens_filc
        wave=wfc3_g280_strct(ii).wavec
        count=wfc3_g280_strct(ii).countsc
        var=wfc3_g280_strct(ii).varc
        sky=wfc3_g280_strct(ii).skyc
     end
     else: stop
  endcase
  
  if keyword_set(OLD) then begin
     mag_set= xmrdfits(sens_fil,1, /silen)
     
     wave_min = mag_set.WAVE_MIN
     wave_max = mag_set.WAVE_MAX
  
     inds = WHERE(wave GE wave_min AND wave LE wave_max, ninds)
     mag_func = bspline_valu(wave[inds], mag_set)
     sens = 10.0^(0.4D*mag_func)
     scale = fltarr(n_elements(wave))
     scale[inds] = sens
  endif else begin
     sens_set=mrdfits(sens_fil,1)
     inds = WHERE(wave GE min(sens_set.wavelength) AND $
                  wave LE max(sens_set.wavelength), ninds)
     sens=fltarr(n_elements(wave))
     sens[inds]=interpol(sens_set.sensitivity,sens_set.wavelength,$
                         wave[inds])
     dw=fltarr(n_elements(wave))    
     dw[inds]=shift(wave[inds],1)-wave[inds]
     dw(inds[0])=dw(inds[1])
     dw=abs(dw)
     scale=1d/(sens*dw)
     scale(where(~finite(scale)))=0D
  endelse

  flux=(count/wfc3_g280_strct(ii).exptime)*scale
  flux_sig=(sqrt(var)/wfc3_g280_strct(ii).exptime)*scale
  fsky=(sky/wfc3_g280_strct(ii).exptime)*scale
  
  case beam of
     0: begin
        wfc3_g280_strct(ii).fluxa=flux
        wfc3_g280_strct(ii).skya=fsky
        wfc3_g280_strct(ii).flux_siga=flux_sig
     end
     1: begin
        wfc3_g280_strct(ii).fluxc=flux
        wfc3_g280_strct(ii).skyc=fsky
        wfc3_g280_strct(ii).flux_sigc=flux_sig
     end
     else: stop
  endcase
end
