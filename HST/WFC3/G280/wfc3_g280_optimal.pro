;+ 
; NAME:
;  wfc3_g280_optimal 
;
; PURPOSE:
;   Simple optimal extraction (Gaussian profile)
;
; CALLING SEQUENCE:
;   
;  wfc3_g280_optimal, wfc3_g280_strct, ii, specim, var, BEAM=, SIGMA=, $
;                       DEBUG=debug, NCRITER=ncriter
;
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;   ii -- the index of the object in the structure
;   specim -- Spectral image (should be sky subtracted)
;   var    -- Variance image
;
; RETURNS:
;
; OUTPUTS:
;  Updates the wfc3_g280 structure to include the optimal flux extraction.
;
; OPTIONAL KEYWORDS:
;  BEAM= -- Beam to reduce
;  SIGMA= -- Sigma for the Gaussian profile
;  DEBUG= -- Debug the program
;  NCRITER= -- CR iterations [default:1L]
;  REJSIG= -- CR rejection criterion [default:9L]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2-May-2016 Written by JXP -- Taken from x_extobjopt
;------------------------------------------------------------------------------
pro wfc3_g280_optimal, wfc3_g280_strct, ii, specim, var, BEAM=beam, SIGMA=sigma, $
                       DEBUG=debug, NCRITER=ncriter, REJSIG=rejsig

  if (N_params() LT 4) then begin 
    print,'Syntax - ' + $
          'wfc3_g280_optimal, specim, var [v1.0]'
    return
  endif 
  
  if not keyword_set(REJSIG) then rejsig = 9.
  if n_elements(NCRITER) EQ 0 then ncriter = 1L

  case beam of
     0: begin
        cnt=wfc3_g280_strct(ii).cnta
        trace_x=wfc3_g280_strct(ii).trace_xa
        trace_y_fit=wfc3_g280_strct(ii).trace_ya_fit
        tsigma=median(wfc3_g280_strct(ii).trace_sigma_fita(0:cnt-1))*1.<2.
        ;tsig=wfc3_g280_strct(ii).trace_sigmaa(0:cnt-1)
        ;tsig(where(abs(tsig) gt 5*median(tsig)))=median(tsig)
        ;sigma=poly(findgen(cnt),poly_fit(findgen(cnt),tsig,8))
     end
     1: begin
        cnt=wfc3_g280_strct(ii).cntc
        trace_x=wfc3_g280_strct(ii).trace_xc
        trace_y_fit=wfc3_g280_strct(ii).trace_yc_fit
        tsigma=median(wfc3_g280_strct(ii).trace_sigma_fitc(0:cnt-1))*1.2<2.
        ;tsig=wfc3_g280_strct(ii).trace_sigmac(0:cnt-1)
        ;tsig(where(abs(tsig) gt 5*median(tsig)))=median(tsig)
        ;sigma=poly(findgen(cnt),poly_fit(findgen(cnt),tsig,8))
     end
     else: stop
  endcase

  if keyword_set(SIGMA) then sigma=sigma else sigma=tsigma
  
  ;; Cut image
  idx = round(trace_x[0]) + lindgen(cnt)
  fx = specim[idx,*]
  cutvar = var[idx,*]
  
  ;; Generate inverse variance
  sz = size(fx, /dim)
  ivar = dblarr(sz[0], sz[1])
  gdp = where(var > 0.)
  ivar[gdp] = 1./cutvar[gdp]
  
  ;; Profile
  xc = replicate(1.,sz[0])#findgen(sz[1])
  ;; Restrict dx according to trace
  ;tot_trc[idx] = trace_y_fit[0:cnt-1]
  dx = (-1.*double(sz[1])) > $
       (xc - trace_y_fit[0:cnt-1]#replicate(1.d,sz[1])) < double(sz[1])
  
  ;P = fltarr(sz)
  ;for kk=0, sz[0]-1 do $
  ;   P[kk,*] =  exp( -(dx[kk,*])^2 / (2.d*sigma[kk]^2) ) * 1.0 / sqrt(2.0 * !pi) / sigma[kk]
  P =  exp( -(dx)^2 / (2.d*sigma^2) ) * 1.0 / sqrt(2.0 * !pi) / sigma
  
  ;; Truncate
  if keyword_set(TRUNCATE) then begin
     outside = where(abs(dx) GT 2*sigma)
     ; inside = where(abs(dx) LT 3*sigma)
     P[outside] = 0.
     ;; Re-normalize
     norm = total(P, 2)
     P = P / (norm ## replicate(1., sz[1]))
     ;stop
  endif
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  msk = bytarr(sz[0],sz[1])
  P2 = P^2

  ;; Optimal with CR rejection
  for qq=0L,NCRITER-1 do begin
      fopt_num = (P * fx * ivar) 
      fopt_den = (P2 * ivar) 
      dd = total(fopt_den,2,/double) 
      gd = where(dd GT 0.)
      fopt_ini = dblarr(sz[0])
      fopt_ini[gd] = (total(fopt_num,2,/double))[gd] / dd[gd]
      fopt_img = fopt_ini # replicate(1., sz[1])

      ;; New Variances
      nwvar = cutvar 
      ;IF KEYWORD_SET(USE_INPUT_VAR) THEN nwvar = var $
      ;ELSE nwvar = readno^2 + P * fopt_img + sky
      ;if nbad NE 0 then nwvar[bad] = 0.  ; Bad pixels
      
      ;; Check for CR (only 1 iteration for now)
      cr = where((fx - P*fopt_img)^2 GT rejsig^2*nwvar, ncr)
  
      ;; Recalculate ivar
      if ncr NE 0 then ivar[cr] = 0.
   endfor

  if keyword_set(DEBUG) then begin
     xatv, ivar, /bloc, min=-0.001, max=0.004
     stop
  endif

  msk = fltarr(sz[0],sz[1]) + 1.
  if ncr NE 0 then msk[cr] = 0.

  ;; Optimal values now
  fopt_num = float(msk * P * fx * ivar) 
  fopt_den = float(msk * P2 * ivar) 
  fvar_num = float(msk * P ) 
  fvar_den = float(msk * P2 * ivar) 

  ;; 1D spectra
  count = total(fopt_num,2) / total(fopt_den,2)
  tvar = total(fvar_num,2) / total(fvar_den,2)

  ;; Check?
  if keyword_set(DEBUG) then begin
     if beam EQ 0 then begin
        x_splot,  findgen(cnt), wfc3_g280_strct(ii).countsa[0:cnt-1], ytwo=count, /bloc
        x_splot,  findgen(cnt), wfc3_g280_strct(ii).vara[0:cnt-1], ytwo=tvar, /bloc
        stop
     endif else begin
        x_splot,  findgen(cnt), wfc3_g280_strct(ii).countsc[0:cnt-1], ytwo=count, /bloc
        x_splot,  findgen(cnt), wfc3_g280_strct(ii).varc[0:cnt-1], ytwo=tvar, /bloc
        stop
     endelse
  endif
  
  ;; Save
  case beam of
     0: begin
        wfc3_g280_strct(ii).countsa=count
        wfc3_g280_strct(ii).vara=tvar
     end
     1: begin
        wfc3_g280_strct(ii).countsc=count
        wfc3_g280_strct(ii).varc=tvar
     end
     else: stop
  endcase

end
