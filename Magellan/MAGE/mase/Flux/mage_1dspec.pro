;+
; NAME:
; mage_1dspec
;    Version 1.1
;
; PURPOSE:
;   Combines orders of a echfspec structure into a 1d spectrum.  
;   Outputs flux and error arrays into separate files (_F.fits, _E.fits)
;
; CALLING SEQUENCE:
;  
;   mage_1dspec, mage, setup, obj_id, side
;
; INPUTS:
;   mage      - MAGE structure
;   obj_id   -  Object ID  (e.g. 0L, 1L, etc)
;   [exp_id] -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;   1d flux      -   (fits file; FSpec/name_ech_F.fits)
;   1d error     -   (fits file; FSpec/name_ech_E.fits)
;
; OPTIONAL KEYWORDS:
;    /SILENT   - No text output
;    OBJ_NM=   - Name of object in slit (a = science)
;    /STD      - Run on a standard star
;    OUTNM=    - Alternative output name for FSpec file
;    ENDTRM=   - Trim order edges by this when calculating ratios
;    MINPIX1=  - Minimum 'good' pixels to calculate fitting profile
;    ORDNM=    - 0th or 1st order fitting (1 = 1st, 2 = 0th (default))
;    SNRMIN=   - Minimum S/N per pixel for consideration in fitting
;    MINPIX2=  - Minimum 'good' pixels to calculate simple 0th order fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   mike_1dspec, mike, setup, obj_id, side
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Jun-2004 Written by GEP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro mage_1dspec,spec,outfile,sigfil,combname, resvel=resvel $
                , REBIN=rebin, hdr=hdr, CHK=chk
  
;;;;;;;;;;;;
;;; begin 1d

  use = where(spec.phys_ordr NE 0,nuse)
  
  gd = where(spec.wave NE 0.)
  minwv = min(spec.wave[gd])
  maxwv = max(spec.wave[gd])

  npix=17720L
  ;npix = 18000L
  velpix=22.0d
  cdelt = velpix/299792.458d/alog(10.0d)
  tot_wave = 10^(alog10(3000.d) + dindgen(npix)*cdelt)
  ;cdelt = alog10(1.0d + velpix/299792.458d)
  ;npix = 50000L
  ;tot_wave = 10^(alog10(3000.0d) + dindgen(npix)*cdelt)
  weight = fltarr(npix)
  weight[*] = 0.
  tot_flux = fltarr(npix)
  tot_flux[*] = 0.
  sig = replicate(0.0,npix)

  if (keyword_set(hdr)) then begin
     head = hdr
  endif else begin
     fxhmake, head
  endelse

  sxaddpar, head, 'CDELT1', cdelt
  sxaddpar, head, 'CRPIX1', 1
  sxaddpar, head, 'CTYPE1', 'LINEAR'
  sxaddpar, head, 'DC-FLAG', 1
  sxaddpar, head, 'BITPIX', -32
  sxaddpar, head, 'NAXIS', 1
  sxaddpar, head, 'NAXIS1', n_elements(tot_flux)
  arr = strsplit(sigfil,"/", /extract)
  sxaddpar, head, 'SIGFILE', arr[n_elements(arr)-1]
  sxdelpar, head, 'NAXIS2'

  if (keyword_set(RESVEL)) then sxaddpar, head, 'RESVEL', resvel
  
  ;; Get S/N for each order overlap, zero out fx, sig in bad blaize areas

  edge_weight = 1.0*(spec.var[*,use] GT 0)

  
  for i=0L,nuse-1L do begin
     cnst = 61750.
                                ;cnst = 61700.
     ordr = spec.phys_ordr[use[i]]
     hup_wv = cnst / (ordr - 0.5)
     hlo_wv = cnst / (ordr + 0.5)
     fsr = hup_wv - hlo_wv
     
     gd = where(spec.wave[*,use[i]] LE hlo_wv and $
                spec.wave[*,use[i]] NE 0. ,ngd)
     if ngd GT 0 then begin
        fixfun = (1.0 - (hlo_wv  - spec.wave[gd,use[i]])/ fsr)  > 0
        edge_weight[gd,i] =  fixfun
     endif
     
     gd = where(spec.wave[*,use[i]] GE hup_wv,ngd)
     if ngd GT 0 then begin
        fixfun = (1.0 - (spec.wave[gd,use[i]] - hup_wv)/ fsr)  > 0
        edge_weight[gd,i] =  fixfun
     endif
     
  endfor
 


  if keyword_set( CHK ) then begin
     for ii=1,nuse-2 do begin
        title = strtrim(string(use[ii]),2)
        u1 = where(spec.wave[*,use[ii]] GT 0)
        x_splot, spec.wave[u1,use[ii]], spec.fx[u1,use[ii]],  ytwo=edge_weight[*,ii], xtwo=spec.wave[u1,use[ii]], psym3=-3,  ythr=spec.var[u1,use[ii]], xthr=spec.wave[u1,use[ii]], /block
                                ;x_splot, spec.wave[u1,use[ii]], spec.fx[u1,use[ii]], psym3=-3, /block
        u2 = where(spec.wave[*,use[ii-1]] GT 0)
        u3 = where(spec.wave[*,use[ii+1]] GT 0)

                                x_splot, spec.wave[u1,use[ii]], spec.fx[u1,use[ii]], ytwo=spec.fx[u2,use[ii-1]], xtwo=spec.wave[u2,use[ii-1]], ythr=spec.fx[u3,use[ii+1]], xthr=spec.wave[u3,use[ii+1]], psym3=-3, /block, TITLE=title
                               ;STOP
     endfor
  endif


  sw = spec.wave[*,use]
  sv = spec.var[*,use]
  sf = spec.fx[*,use]
  
;  This is ONLY to be used for the
;  quick look pipeline which employs
;  boxcar extraction.
  

   if keyword_set(rebin) then begin
      for i=0, nuse-1 do begin

         logwv = alog10(3000.d)
         while (10^logwv LT sw[0,i]) do logwv += cdelt

         indrun = where(sw[*,i] NE 0, ngd)

         sw_rebin = 0.0*sw[*,i]         
         sf_rebin = 0.0*sf[*,i]         
         sv_rebin = 0.0*sv[*,i]         

         sw_rebin[indrun] = 10^(logwv+cdelt*dindgen(ngd))

         sset = bspline_iterfit(sw[indrun,i],sf[indrun,i],everyn=2.5, /silent)
         sf_rebin[indrun] = bspline_valu(sw_rebin[indrun], sset)

         sset = bspline_iterfit(sw[indrun,i],sv[indrun,i],everyn=2.5, /silent)
         sv_rebin[indrun] = bspline_valu(sw_rebin[indrun], sset)

         sw[*,i] = sw_rebin
         sf[*,i] = sf_rebin
         sv[*,i] = sv_rebin

      endfor
   endif
   
   gd = where(sw GT 1000.0 AND sv GT 0 AND sf NE 0.0000, ngd)

   sg = gd[sort(sw[gd])]
   min_ind = min(abs(tot_wave - min(sw[gd])), min_ii)
   max_ind = min(abs(tot_wave - max(sw[gd])), max_ii)
   ;min_ii = (where(abs(tot_wave-min(sw[gd])) LT 0.02))[0]
   ;max_ii = (where(abs(tot_wave-max(sw[gd])) LT 0.02))[0]
   ll = 0L
   for ii=min_ii,max_ii do begin
      uu = (ll + 6L) <  (ngd-1L)
      ;Also, check this for the wave info...
     ;a = where(abs((tot_wave[ii]-sw[sg[ll:uu]])/tot_wave[ii]*299792.458d) LT 3*velpix,na)
     ;a = where(abs((tot_wave[ii]-sw[sg[ll:uu]])/tot_wave[ii]*299792.458d) LE 0.01*velpix,na)
     ;a = where(abs((tot_wave[ii]-sw[sg[ll:uu]])) LT 0.02,na)
     a = where(float(tot_wave[ii]) EQ float(sw[sg[ll:uu]]), na)
;;print, na
      if na EQ 0 then continue
      
      if na EQ 1 then begin
         tot_flux[ii] = sf[sg[a[0]+ll]] 
         sig[ii] = sqrt(sv[sg[a[0]+ll]])
        ; STOP
      endif else begin
         h = sg[a + ll]
         wtmp = edge_weight[h] / sv[h]
         wtot = total(wtmp) 
      ;   STOP
         if wtot GT 0 then begin
            
;  Flux = Sum (W * fx) / Sum(W)
;  Var = Sum(W^2 * var) / Sum(W)^2
 ;are these correct?           
            tot_flux[ii] = total(wtmp * sf[h]) / wtot
            sig[ii] =  sqrt(total(wtmp^2 * sv[h]) / total(wtmp^2))
       ;     STOP

            ;trying a new mask here
            ;outlier = where(sf[h] GT tot_flux[ii] + 8.0*sqrt(sv[h]), $
            ;                complement=keep)
            ;outlier = where(abs(sf[h]) GT abs(median(sf[h]))*30., $
            ;               complement=keep)
            ;keep = INDGEN(N_ELEMENTS(h))
            ;if outlier[0] NE -1 then begin
            ;   if n_elements(keep) EQ 1 then begin
            ;      tot_flux[ii] = sf[h[keep[0]]]
            ;      sig[ii] = sqrt(sv[h[keep[0]]])
            ;   endif else begin
            ;      h = h[keep]
                  wtmp = edge_weight[h] / sv[h]
                  wtot = total(wtmp) 
                  tot_flux[ii] = total(wtmp * sf[h]) / wtot
                  sig[ii] =  sqrt(total(wtmp^2 * sv[h])/ total(wtmp^2))               
          ;        STOP           
;    endelse
            endif
        ; endif
      endelse
      
      ll = max(a) + ll
      
   endfor

    minwv = tot_wave[min_ii]
    tot_pix = max_ii - min_ii + 1L
    a = lindgen(tot_pix) + min_ii
   
    if keyword_set( CHK ) then $
      x_splot, tot_wave[a], tot_flux[a], ytwo=sig[a], /block

    sxaddpar, head, 'CRVAL1', alog10(minwv)
    ;sxaddpar, head, 'CRVAL1', alog10(3000.d)


    
    if not keyword_set( SILENT ) then $
      print, 'mage_1dspec: Writing results to files: ', outfile, ',', sigfil

    mwrfits, tot_flux[a], outfile, head, /create, /silent
    mwrfits, sig[a], sigfil, head, /create, /silent

    mwrfits, tot_flux[a], combname, /create, /silent
    mwrfits, sig[a], combname, head,  /silent

;STOP


    return
end
