@civ_calcewn
@civ_sensitivity

function civ_omciv_fit, intlim, coeff, alpha, n_norm, siiv=siiv, $
                        sigomciv=sigomciv,$
                        sigcoeff=sigcoeff, sigalpha=sigalpha, fit_fil=fit_fil,$
                        errstrct=errstrct, cosmology=cosmology
  
  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, cosm_L, cosm_r
  if keyword_set(cosmology) then $
     cosm_common, H0=cosmology[0], Omegavac=cosmology[2], $
                  OmegaDM=cosmology[1],/silent $
  else $
     cosm_common, H0=70., Omegavac=0.7, OmegaDM=0.3,/silent

  ;; Constants
  rhocrit = 1.89e-29 * (cosm_h/100.d)^2      ; g/cm^3
  h_invs = double(cosm_h) / (1.e6 * 3.09e13) ; s^-1
  amu = double(1.66053886e-24)                 ; g
  if keyword_set(siiv) then $
     mciv = 28.0855d * amu $
  else mciv = 12.0107d * amu    ; g 
  sol = double(2.998e10)         ; cm/s
  om_coeff = h_invs * mciv / (sol * rhocrit)

  if keyword_set(fit_fil) then begin
     if size(fit_fil,/type) eq 7 then begin
        fit = xmrdfits(fit_fil,1,/silent) 
        errstrct = xmrdfits(fit_fil,3,/silent) 
     endif else fit = fit_fil

     alpha = fit.alpha
     if not keyword_set(sigalpha) then sigalpha = fit.sigalpha
     coeff = fit.coeff
     if not keyword_set(sigcoeff) then sigcoeff = fit.sigcoeff
     n_norm = fit.datanorm
     
  endif 

  if not keyword_set(intlim) then intlim = [10.^13,10.^15]
  if not keyword_set(n_norm) then begin
     if keyword_set(siiv) then n_norm = 10.d^13.5 $
     else n_norm = 10.d^14. 
  endif 

  alpha2 = alpha + 2.
  deltan = intlim[1]^alpha2 - intlim[0]^alpha2
  omciv = om_coeff * coeff/alpha2 * deltan / n_norm^alpha

  ;; Error analysis (remember: alpha and coeff anti-correlated;
  ;; so must stay on error ellipse)
  sigomciv = dblarr(2)
  if keyword_set(sigcoeff) and keyword_set(sigalpha) then begin
     if keyword_set(errstrct) then begin
        ;; This is best way
        omciv_surf = dblarr(n_elements(errstrct.surf)) ; minimize computation
        tmp = fit
        logL_grid = reform(errstrct.logL,n_elements(errstrct.coeff_grid),$
                           n_elements(errstrct.alpha_grid)) ;  for reference
        for nn=0L,n_elements(errstrct.surf)-1 do begin
           mm = array_indices(logL_grid,errstrct.surf[nn])
           omciv_surf[nn] = civ_omciv_fit(intlim,errstrct.coeff_grid[mm[0]],$
                                          errstrct.alpha_grid[mm[1]],$
                                          fit.datanorm,siiv=siiv)
        endfor
        omciv_mn = min(omciv_surf,imn,max=omciv_mx,subscript_max=imx) ; find extrema
        sigomciv[0] = omciv - omciv_mn
        sigomciv[1] = omciv_mx - omciv
     endif else begin
        ;; This is wrong
        ;; Note: if wanting to reproduce the limits with 
        ;; fistrct.sigcoeff_omciv, must use -shift(fitstrct.sigcoeff_omciv,1)
        ;; Lower
        alpha2 = alpha - sigalpha[0] + 2. ; makes deltan smaller
        deltan = intlim[1]^alpha2 - intlim[0]^alpha2 
        omcivlo = om_coeff * (coeff+sigcoeff[1])/alpha2 * deltan / $
                  n_norm^(alpha - sigalpha[0])
        sigomciv[0] = omciv - omcivlo
        
        ;; Upper
        alpha2 = alpha + sigalpha[1] + 2. ; makes deltan larger
        deltan = intlim[1]^alpha2 - intlim[0]^alpha2
        omcivhi = om_coeff * (coeff-sigcoeff[0])/alpha2 * deltan / $
                  n_norm^(alpha + sigalpha[1])
        sigomciv[1] = omcivhi - omciv
     endelse 
  endif else begin
     if keyword_set(sigcoeff) or keyword_set(sigalpha) then $
        print,'civ_omciv(): sigcoeff and sigalpha must both be set'
  endelse 

  return,omciv

end                             ; civ_omciv_fit()


function civ_omciv_summ,ncolm,xpath,siiv=siiv,$
                        signcolm=signcolm,sigx=sigx,sigomciv=sigomciv,$
                        civstrct=civstrct,sens_fil=sens_fil,nciv=nciv,$
                        cosmology=cosmology,_extra=extra
  ;;;;;;;;;;;;;
  ;; Sum data
  ;;;;;;;;;;;;;
  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, cosm_L, cosm_r
  if keyword_set(cosmology) then $
     cosm_common, H0=cosmology[0], Omegavac=cosmology[2], $
                  OmegaDM=cosmology[1],/silent $
  else $
     cosm_common, H0=70., Omegavac=0.7, OmegaDM=0.3,/silent

  ;; Constants
  rhocrit = 1.89e-29 * (cosm_h/100.)^2 ; g/cm^3
  h_invs = cosm_h / (1.e6 * 3.09e13)   ; s^-1
  amu = 1.66053886e-24                 ; g
  if keyword_set(siiv) then begin
     mciv = 28.0855 * amu 
     dblt_name = 'SiIV'
  endif else begin
     mciv = 12.0107 * amu       ; g 
     dblt_name = 'CIV'
  endelse 
  sol = 2.998e10                ; cm/s
  om_coeff = h_invs * mciv / (sol * rhocrit)

  ;; Check
  if n_params() lt 2 and not keyword_set(civstrct) $
     and not keyword_set(sens_fil) then $
     stop,'civ_omciv_summ(): wrong parameter combination'

  if not keyword_set(ncolm) then begin
     if size(civstrct,/type) eq 7 then civ_group,civcand,civstrct,_extra=extra $
     else civcand = civstrct
     ncolm = civ_calcewn_ndblt(civcand,dblt_name,$
                               signcolm=signcolm,flg_colm=flg_colm,/silent)
  endif 

  if not keyword_set(xpath) then begin
     if not keyword_set(sens_fil) then $
        stop,'civ_omciv_summ(): no way to generate xpath'
;     if size(sens_fil,/type) eq 7 then $
        xpath = civ_sensitivity_x(sens_fil,ncolm=alog10(ncolm),$
                                  signcolm=signcolm,sigx=sigx) ;$
;     else xpath = sens_fil
  endif 

  omciv = om_coeff * total(ncolm/xpath)
  sigomciv = dblarr(2)
  if keyword_set(signcolm) and keyword_set(sigx) then begin
     ;; X anti-correlated
     ;; Brute-force method
     sigomciv[0] = om_coeff * sqrt(total((signcolm/xpath)^2) + $
                                   total((ncolm/xpath^2*sigx[*,1])^2))
     sigomciv[1] = om_coeff * sqrt(total((signcolm/xpath)^2) + $
                                   total((ncolm/xpath^2*sigx[*,0])^2))
  endif 

return,omciv

end                             ; civ_omciv_summ()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_omciv,siiv=siiv,intlim=intlim,compare_lim=compare_lim,_extra=extra
  if not keyword_set(intlim) then intlim = [10.^13,10.^15]
  zmed = 0.6

  fmt = '(a18,1x,f7.5,2x,e9.3,1x,e9.3,1x,e9.3,2x,e9.3,1x,e9.3,1x,e9.3)'


  for ii=0,5 do begin
     ;;;;;;;;;;;;;
     ;; Data
     ;;;;;;;;;;;;;
     case ii of 
        0: begin
           ;; Rating = 6
           fit_fil = getenv('MLSS_DIR')+'/analysis/inputs/fn_rtgeq6_fit.fits'
           sens_fil = getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal.fits'
           civ_group,civcand,rating=[6],flg_sys=384,nciv=nciv,$
                     siiv=siiv,_extra=extra ; unsat
           lbl = 'Rating = 6 ('+strtrim(nciv,2)+'): '
        end
        1: begin
           ;; z < 0.6 
           fit_fil = getenv('MLSS_DIR')+'/analysis/inputs/fn_rtgeq6_zlt0.6fit.fits'
           sens_fil = getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal_zlt0.6.fits'
           civ_group,civcand,rating=[6],flg_sys=384,nciv=nciv,zlim=[0.,zmed],$
                     siiv=siiv,_extra=extra
           lbl = 'z < 0.6 ('+strtrim(nciv,2)+'): '
        end
        2: begin
           ;; z >= 0.6 
           fit_fil = getenv('MLSS_DIR')+'/analysis/inputs/fn_rtgeq6_zgt0.6fit.fits'
           sens_fil = getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal_zgt0.6.fits'
           civ_group,civcand,rating=[6],flg_sys=384,nciv=nciv,zlim=[zmed,10.],$
                     siiv=siiv,_extra=extra
           lbl = 'z >= 0.6 ('+strtrim(nciv,2)+'): '
        end
        3: begin
           ;; Rating >= 5
           fit_fil = getenv('MLSS_DIR')+'/analysis/inputs/fn_rtg56_fit.fits'
           sens_fil = getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal.fits'
           civ_group,civcand,rating=[5,6],flg_sys=384,nciv=nciv,$
                     siiv=siiv,_extra=extra
           lbl = 'Rating >= 5 ('+strtrim(nciv,2)+'): '
        end
        4: begin
           ;; z < 0.6 
           fit_fil = getenv('MLSS_DIR')+'/analysis/inputs/fn_rtg56_zlt0.6fit.fits'
           sens_fil = getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal_zlt0.6.fits'
           civ_group,civcand,rating=[5,6],flg_sys=384,nciv=nciv,zlim=[0.,zmed],$
                     siiv=siiv,_extra=extra
           lbl = 'z < 0.6 ('+strtrim(nciv,2)+'): '
        end
        5: begin
           ;; z >= 0.6 
           fit_fil = getenv('MLSS_DIR')+'/analysis/inputs/fn_rtg56_zgt0.6fit.fits'
           sens_fil = getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal_zgt0.6.fits'
           civ_group,civcand,rating=[5,6],flg_sys=384,nciv=nciv,zlim=[zmed,10.],$
                     siiv=siiv,_extra=extra
           lbl = 'z >= 0.6 ('+strtrim(nciv,2)+'): '
        end
        else: begin
        end
     endcase   
     ncolm = 0                  ; reset
     omciv_summ = civ_omciv_summ(ncolm,siiv=siiv,sigomciv=sigomciv_summ,/signcolm,$
                                 civstrct=civcand,sens_fil=sens_fil,nciv=nciv2,$
                                 _extra=extra) ; unsat, etc

     fitstrct = xmrdfits(fit_fil,1,/silent)
     errstrct = xmrdfits(fit_fil,3,/silent)
     if keyword_set(compare_lim) then begin
        intlim = [min(ncolm),max(ncolm)] ; returned by civ_omciv_summ()
     endif  
     omciv_fit = civ_omciv_fit(intlim,siiv=siiv,fit_fil=fitstrct,$
                               sigomciv=sigomciv_fit,errstrct=errstrct)

     ;;;;;;;;;;;;;
     ;; Print
     ;;;;;;;;;;;;;
     print,lbl,median(civcand.zabs[0]),omciv_summ,sigomciv_summ[0],sigomciv_summ[1],$
           omciv_fit,sigomciv_fit[0],sigomciv_fit[1],format=fmt

  endfor                        ; loop all cases
end
