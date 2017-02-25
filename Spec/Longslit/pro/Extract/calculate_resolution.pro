;+
; NAME:
;   calculate_resolution
;
; PURPOSE:
; calculate_resolution, s, anamorph=anamorph
;
; CALLING SEQUENCE:
;   long_superbias, filenames, outfile, [ sigrej=, maxiter=, /verbose ]
;                 
;
; INPUTS:
;   s -- Spectral structure which contains all the key quantities
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  
; EXAMPLES:
;
; BUGS:
;    
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Mar-2005  Written by J. Hennawi (UCB), D. Schlegel (LBL)
;-
;------------------------------------------------------------------------------
pro calculate_resolution, s, anamorph=anamorph

   if NOT keyword_set(anamorph) then anamorph = 1.0


   nslit =  n_elements(s)

   for i=0, nslit-1 do begin
     spatial_binning    = s[i].binning[0]
     dispersion_binning = s[i].binning[1]
     h = median(s[i].arc_fwhm)
     median_slit_qm     = median(s[i].arc_fwqm)
   
;
;  2/3 is due to weird definition in dierfc, represents erf^{-1}(1/3)
;
     b = ((median_slit_qm - h)) / dierfc(2./3.) / 2.
     b = (b > 0.5) < (0.75 * h)


     splog, 'Slit #', i, ', Half Max: ', Erf(h/b)/2., $
                            ', 1/3 Max: ', Erf((median_slit_qm + h)/2/b)/2. - $
                                           Erf((median_slit_qm - h)/2/b)/2. 

     seeing_vs_slit = s[i].fwhmfit / anamorph * spatial_binning $
                          / dispersion_binning / 2.3548 / ( h / 2.)

     min_w = min(seeing_vs_slit)
     max_w = max(seeing_vs_slit) 
    
     fwhm_temp = fltarr(11)
     this_sigma = fltarr(11)
     x = findgen(5001)*0.002 - 5.0
     slit_start = ((abs(x)+0.005) < 1.0) - ((abs(x) - 0.005) < 1.0)

     bkern = gauss_kernel(100.*b/sqrt(2.)/(h/2.))

     for ii=0,10 do begin
       this_sigma[ii] = min_w + 0.1*ii * (max_w - min_w)
       ss = slit_start * gauss_kernel(100.*this_sigma[ii], hpix=2500)
       result = convol(ss, bkern)

       hm = min(where(result GE 0.5*max(result))) > 0
       fwhm_temp[ii] = 0.01*h*((500 - hm) + (result[hm] - 0.5*max(result))/$
                        (result[hm] -result[hm-1]))
       if ii EQ 0 then splog, 'Min FWHM is ', fwhm_temp[ii]/h, $
                              ' percent of slit width'
       if ii EQ 10 then splog, 'Max FWHM is ', fwhm_temp[ii]/h, $
                               ' percent of slit width'
     endfor

     kk = sort(seeing_vs_slit)
     linterp, this_sigma, fwhm_temp, seeing_vs_slit[kk], ans
     s[i].pix_res[kk] = ans

  endfor

end 
       




