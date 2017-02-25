 ;+ 
; NAME:
; sdss_ewciv
;  Version 1.1
;
; PURPOSE:
;  Measures the EW of the CIV lines (rest values) in the SDSS
;
; CALLING SEQUENCE:
;  sdss_ewciv, wave, flux, sig, dblt_name, qstrct=, ZABS=, PLOT=
;
; INPUTS:
;  wave  -- Wavelength array
;  flux  -- Flux array 
;  sig   -- Sigma array
;  conti -- continuum
;  dblt_name -- Doublet name or structure (compatible with
;               dblt_retrieve()); not used with /generic
;  zabs  -- Absorption redshift; unless /generic set, then 
;           expect centroid (in Angstrom)
;
; RETURNS:
;  qstrct -- CIV structure with CIV EW filled up
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  istrct= -- input {qalcharstrct} formatted structure for candidate
;            doublets
;  zlim= -- redshift bounds for transition (helps define boundaries);
;           not used with /generic
;  snr_conv= -- convolved SNR "spectrum" from sdss_fndlin (helps
;               define bounds)
;  ncsiglim= -- number of flux error above continuum to include; clipping
;  /contilim -- cap flux at continuum
;  nzsiglim= -- number of flux error below zero to include; clipping
;  /zerolim -- floor flux at zero
;  /plot -- show some stuff
;  /debug -- print some stuff
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
;   27-Feb-2004 Written by JXP
;   27-May-2011 Make generic, KLC
;   30-May-2011 Work with new sdsscivstrct, KLC
;   06-Jun-2011 Turned into function and allowed input, KLC
;   10-Nov-2011 Added sdss_ewciv_strct, KLC
;   25-Jan-2013 Allow inflation of boundaries for AODM, KLC
;   05-Jun-2013 Force upper limit to flux, no excluding pixels, KLC
;   23-Jul-2015 Added nzsiglim and zerolim, no default clipping, KLC
;-
;------------------------------------------------------------------------------

pro sdss_ewciv_strct, civ_fil, out_fil, dblt_name=dblt_name, clobber=clobber, $
                      _extra=extra
  ;; Measure better EW for all doublets; allow widths to grow
  if n_params() ne 2 then begin
     print,'Syntax - sdss_ewciv_strct, civ_fil, out_fil, [dblt_name=, /clobber, _extra=]'
     return
  endif 

  sdssdir = sdss_getsdssdir()
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 7 then dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name

  if size(civ_fil,/type) eq 7 then $
     civstr = xmrdfits(civ_fil,1,/silent) $
  else civstr = civ_fil
  nciv = (size(civstr,/dim))[0] > 1 ; avoid singularity
  
  ;; Avoid reading in too much
  prev_spec = ''
  for iciv=0L,nciv-1 do begin

     if civstr[iciv].sdss_obs[0] ne prev_spec then begin
        ;; Read spectrum
        parse_sdss,sdssdir+civstr[iciv].sdss_obs[0],flux,wave,npix=npix,sig=sigma
        ;; Read conti
        cstrct = xmrdfits(sdssdir+civstr[iciv].abslin_fil,1,/silent)
        cindx = fix(alog(civstr[iciv].cflg)/alog(2)) ; assume all same cflg in LOS
        sig = sdss_calcnormerr(flux,sigma,cstrct,/unnorm)
     endif                      ; read in
     
     ;; _extra includes /debug, /plot, /keepwvlim, /final, /generic,
     ;; /contilim, ncsiglim=, /zerolim, nzsiglim=
     tmpstr = sdss_ewciv(wave[cstrct.ipix0:cstrct.npix-1], $
                         flux[cstrct.ipix0:cstrct.npix-1], $
                         sig[cstrct.ipix0:cstrct.npix-1], $
                         cstrct.conti[cstrct.ipix0:cstrct.npix-1,cindx],$
                         dblt,civstr[iciv].zabs_orig[0],$
                         snr_conv=cstrct.snr_conv[cstrct.ipix0:cstrct.npix-1,cindx],$
                         istrct=civstr[iciv],_extra=extra)
     civstr[iciv] = tmpstr

     if iciv ne 0 and (iciv mod 500) eq 0 then $
        print,'sdss_ewciv_strct: iciv = ',iciv

  endfor                        ; loop iciv=nciv

  test = file_search(out_fil+'*',count=ntest)
  if ntest eq 0 or keyword_set(clobber) then begin
     mwrfits,civstr,out_fil,/create,/silent
     spawn,'gzip -f ' + out_fil
     print,'sdss_ewciv_strct: created ',out_fil
  endif else $
     stop,'sdss_ewciv_strct stop: will not overwrite file ',out_fil

end                             ; sdss_ewciv_strct

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function sdss_ewciv, wave, flux, sig, conti, dblt_name, zabs, istrct=istrct, $
                     zlim=zlim, plot=plot, debug=debug, snr_conv=snr_conv, $
                     keepwvlim=keepwvlim, final=final, contilim=contilim, $
                     ncsiglim=ncsiglim, zerolim=zerolim, nzsiglim=nzsiglim, $
                     generic=generic

  if  N_params() LT 6  then begin 
     print,'Syntax - ' + $
           'sdss_ewciv(wave, flux, sig, conti, dblt_name, zabs, '
     print,'                  [zlim=, snr_conv=, istrct=, /plot, /debug, '
     print,'                  /contilim, ncsiglim=, /zerolim, nzsiglim=, /generic])'
     return,-1
  endif 

  ;; Structure
  if keyword_set(istrct) then qstrct = istrct $
  else begin
     qstrct = {sdsscivstrct}
     if keyword_set(keepwvlim) then $
        stop,'sdss_ewciv: must set istrct= to /keepwvlim.'
  endelse 
  tags = tag_names(qstrct)
  wrtag = (where(tags eq 'WREST'))[0]
  if wrtag eq -1 then begin
     ;; Using sdsscontistrct 
     wrtag = where(tags eq 'CENTROID') 
     cindx = fix(alog(qstrct.cflg)/alog(2))
  endif else cindx = 0

  if keyword_set(final) and not keyword_set(generic) then begin
     ;; Use final zabs and wvlim
     ztag = (where(tags eq 'ZABS_FINAL'))[0]
     sigztag = (where(tags eq 'SIGZABS_FINAL'))[0]
     ewtag = (where(tags eq 'EW_FINAL'))[0]
     sigewtag = (where(tags eq 'SIGEW_FINAL'))[0]
     wvlimtag = (where(tags eq 'WVLIM_FINAL'))[0]
     ncolmtag = (where(tags eq 'NCOLM_FINAL'))[0]
     signcolmtag = (where(tags eq 'SIGNCOLM_FINAL'))[0]
  endif else begin
     ;; Use original zabs and wvlim
     ztag = (where(tags eq 'ZABS_ORIG'))[0]
     sigztag = (where(tags eq 'SIGZABS_ORIG'))[0]
     ewtag = (where(tags eq 'EW_ORIG'))[0]
     sigewtag = (where(tags eq 'SIGEW_ORIG'))[0]
     wvlimtag = (where(tags eq 'WVLIM_ORIG'))[0]
     ncolmtag = (where(tags eq 'NCOLM_ORIG'))[0]
     signcolmtag = (where(tags eq 'SIGNCOLM_ORIG'))[0]
  endelse 

  ;; dwv
  dwv = wave - shift(wave,1)
  dwv[0] = dwv[1]
  npix = (size(wave,/dim))[0]

  ;; possible sigma clipping
  flux0 = flux                  ; safety; for plotting and 'return'
  if keyword_set(contilim) then $
     flux = flux < conti $      ; aggressive
  else $
     if keyword_set(ncsiglim) then $
        flux = flux < (conti + ncsiglim*abs(sig))
  ;; default is to not change flux
  
  if keyword_set(zerolim) then $
     flux = flux > 0. $         ; aggressive
  else $
     if  keyword_set(nzsiglim) then $
        flux = flux > (0. - abs(nzsiglim*sig)) ; could be confusing whether nzsiglim +/-
  ;; else can keep flux as it is

  
  ;; Wavelengths
  if keyword_set(generic) then begin
     nlin = n_elements(zabs)               ; size() could return 0
     qstrct.(wrtag)[0:nlin-1,cindx] = zabs ; not to cause confusion
     zabs = zabs * 0.                      ; will reset at end
     boxmax = 20.                          ; something large
     obswv = qstrct.(wrtag)[0:nlin-1,cindx]
     
     dblt = dblt_retrieve('')   ; dummy
     if not keyword_set(snr_conv) then $
        print,'sdss_ewciv(): consider setting snr_conv in /generic mode'

  endif else begin
     if size(dblt_name,/type) eq 8 then dblt = dblt_name $
     else  dblt = dblt_retrieve(dblt_name)
     qstrct.(wrtag)[0:1,cindx] = [dblt.wvI,dblt.wvII]
     
     obswv = (1.+zabs)*qstrct.(wrtag)[0:1,cindx]
     boxmax = 0.5*(obswv[1] - obswv[0]) ; can't grow beyond
     nlin = n_elements(obswv)

     ;; Set box size 
     if dblt.ion eq 'CIV' then $
        ;; Just use natural width of doublet
        boxw = boxmax $ 
     else begin
        ;; For anything else,default to 5 pixels centered on absorber
        mn = min(abs(wave-obswv[0]),icen)
        idx = icen - 2 + lindgen(5)
        boxw = 2.5*mean(dwv[idx]) ; could fail if idx out of bounds
     endelse 

     if keyword_set(zlim) then begin
        ;; Do comparison and take whichever is bigger (and use dblt.wvII
        ;; because that will be bigger)
        nwboxw = 0.5 * dblt.wvII*(zlim[1]-zlim[0]) 
        if keyword_set(debug) and nwboxw gt boxw then $
           print,'sdss_ewciv debug: zlim sets wavelength bounds'
        boxw = boxw > nwboxw
     endif 

  endelse 

  if ztag ne -1 then $
     qstrct.(ztag)[0:nlin-1] = zabs

  if keyword_set(plot) or keyword_set(debug) then clr = getcolor(/load)

  ;; Loop
  for ii=0L,nlin-1 do begin

     if (ii gt 1) or keyword_set(generic) then begin
        ;; Need to set bounds, becasue who knows what it is and how
        ;; wide it should be, and qstrct.(wrtag)[ii,cindx] should already
        ;; be set
        if qstrct.(wrtag)[ii,cindx] le 0. then $
           stop,'sdss_ewciv(): must instantiate wrest['+strtrim(ii,2)+']'
        mn = min(abs(wave-obswv[ii]),icen)
        idx = icen - 2 + lindgen(5)
        boxw = 2.5*mean(dwv[idx]) ; could fail if indx out of bounds
     endif 

     if (ii eq 0 or keyword_set(generic)) and $
        (keyword_set(plot) or keyword_set(debug)) and $
        not keyword_set(keepwvlim) then begin
        ;; Set up plot window
        mn = min(abs(wave-(obswv[ii]-boxw)),ilhs)
        if keyword_set(generic) then $
           mn = min(abs(wave-(obswv[ii]+boxw)),irhs) $
        else mn = min(abs(wave-(obswv[ii+1]+boxw)),irhs)
        
        i1a = (ilhs - 15) > 0L
        i2b = (irhs + 20) < (npix - 1)
        mnw_lhs = mean(wave[i1a:ilhs])
        mnw_rhs = mean(wave[irhs:i2b])
        plot, wave, flux0, xrange=[mnw_lhs-10,mnw_rhs+10], $
              yrange = [0., max(flux0[i1a:i2b])*1.1], $
              backgroun=clr.white, color=clr.black, psym=10
        oplot, wave, conti, color=clr.blue, linest=2, thick=3
     endif
     
     
     ;; EW Edges
     if keyword_set(keepwvlim) then begin
        ;; Use limits set in structure already
        mn = min(abs(wave-qstrct.(wvlimtag)[ii,0]),ilhs)
        mn = min(abs(wave-qstrct.(wvlimtag)[ii,1]),irhs)
     endif else begin
        ;; Pre-defined bounds and maybe adjust
        mn = min(abs(wave-obswv[ii]),icen)
        mn = min(abs(wave-(obswv[ii]-boxw)),ilhs)
        mn = min(abs(wave-(obswv[ii]+boxw)),irhs)
        if keyword_set(snr_conv) then begin
           ;; Roll away from default bounds until convolved S/N
           ;; "spectrum" reaches minimum
           ;; Exclude badness  
           ;; Find where LHS pixels are less than the ones more to the
           ;; left 
           gd = where(snr_conv[0:ilhs] lt shift(snr_conv[0:ilhs],1) $
                      and finite(snr_conv[0:ilhs]) eq 1,ngd)
           if ngd gt 0 then $   ; but must be real part of spectrum
              if sig[gd[ngd-1]] gt 0. then ilhs = gd[ngd-1]

           ;; Find where RHS pixels are less than the ones more to the
           ;; right 
           gd = where(snr_conv[irhs:*] lt shift(snr_conv[irhs:*],-1) $
                      and finite(snr_conv[irhs:*]) eq 1,ngd)
           if ngd gt 0 then $
              if sig[irhs + gd[0]] gt 0. then irhs = irhs + gd[0]

           ;; Also set width, which will float for doublets
           boxw = 0.5*(wave[irhs]-wave[ilhs])

           if boxw gt boxmax then begin
              ;; Can't let this get too big, just set it to be
              ;; the difference.
              if keyword_set(generic) then begin
                 mn = min(abs(wave-(obswv[ii]-boxmax)),imn)
                 ilhs = ilhs > imn ; stay tight to centroid
                 mn = min(abs(wave-(obswv[ii]+boxmax)),imn)
                 irhs = irhs < imn
              endif else begin
                 if ii eq 0 then begin
                    ;; Assume dblt.wvI and don't grow to the red
                    mn = min(abs(wave-(obswv[ii]+boxmax)),imn)
                    irhs = irhs < imn
                 endif 
                 if ii eq 1 then begin
                    ;; Assume dblt.wvII and don't grow to blue
                    mn = min(abs(wave-(obswv[ii]-boxmax)),imn)
                    ilhs = ilhs > imn ; stay tight to centroid
                 endif 
              endelse 
              
              if keyword_set(generic) and keyword_set(debug) then $
                 print,'sdss_ewciv() debug: bounds grew beyond 20 Ang'
           endif 

           qstrct.(wvlimtag)[ii,*] = [wave[ilhs],wave[irhs]] ;save
        endif                                                ; snr_conv= set


        ;; Also do a test to see if should extend further (though
        ;; do not want to enclose another line), by comparing flux to
        ;; actual continua
        ;; Including the +/-sigma check makes sure the boundaries
        ;; doesn't beecome ridiculous; and use original flux because
        ;; that's how it's calibrated 
        bd = where(flux0[0:icen]+sig[0:icen] gt conti[0:icen] and conti[0:icen] ne 0 and $
                   finite(conti[0:icen]),nbd)
        if nbd ne 0 then begin
           if sig[bd[nbd-1]] gt 0. then ilhs_test = bd[nbd-1] $
           else ilhs_test = ilhs
        endif else ilhs_test = ilhs
        bd = where(flux0[icen:*]+sig[icen:*] gt conti[icen:*] and conti[icen:*] ne 0 and $
                   finite(conti[icen:*]),nbd)
        if nbd ne 0 then begin
           if sig[icen + bd[0]] gt 0. then irhs_test = icen + bd[0] $
           else irhs_test = irhs
        endif else irhs_test = irhs

        boxw_test = 0.5*(wave[irhs_test]-wave[ilhs_test])
        if boxw lt boxmax then begin
           case ii of           ; set min/max limits based on other detected lines
              0: begin
                 wvlo = wave[0]     
                 if ztag ne -1 then $
                    wvhi = (1.+qstrct.(ztag)[ii+1])*qstrct.(wrtag)[ii+1,cindx] $
                 else wvhi = qstrct.(wrtag)[ii+1,cindx]
              end
              nlin-1: begin
                 wvhi = wave[npix-1]
                 if ztag ne -1 then $
                    wvlo = (1.+qstrct.(ztag)[ii-1])*qstrct.(wrtag)[ii-1,cindx] $
                 else wvlo = qstrct.(wrtag)[ii-1,cindx]
              end
              else: begin
                 if ztag ne -1 then begin
                    wvlo = (1.+qstrct.(ztag)[ii-1])*qstrct.(wrtag)[ii-1,cindx] 
                    wvhi = (1.+qstrct.(ztag)[ii+1])*qstrct.(wrtag)[ii+1,cindx] 
                 endif else begin
                    wvlo = qstrct.(wrtag)[ii-1,cindx]
                    wvhi = qstrct.(wrtag)[ii+1,cindx]
                 endelse 
              end
           endcase
           
           ;; Accept new bounds so long as they do not enclose another
           ;; line and do not shrink the width
           if ilhs_test lt ilhs and wave[ilhs_test] gt wvlo then ilhs = ilhs_test
           if irhs_test gt irhs and wave[irhs_test] lt wvhi then irhs = irhs_test

           ;; Set up necessary numbers
           qstrct.(wvlimtag)[ii,*] = [wave[ilhs],wave[irhs]] ;save
           boxw = 0.5*(wave[irhs]-wave[ilhs])
        endif                   ; boxw < boxmax
        
     endelse                    ; /keepwvlim not set


     ;; (Re)measure redshift; flux-weighted for good (error != 0)
     ;; and absorption (flux < 1)
     wgt = replicate(0.,irhs-ilhs+1)
     gd = where(sig[ilhs:irhs] ne 0.,ngd) ; $
     if gd[0] ne -1 then begin
        ;; Make a new measurement
        wgt[gd] = 1. - flux[ilhs + gd]/conti[ilhs + gd]
        if keyword_set(debug) then $
           print,'sdss_ewciv debug: flux weighting'

        ;; EW measurement
        qstrct.(ewtag)[ii] = total( (1.-flux[ilhs + gd]/conti[ilhs + gd]) $
                                    *dwv[ilhs + gd]) 
        qstrct.(sigewtag)[ii] = sqrt(total( (sig[ilhs + gd]/conti[ilhs + gd] $
                                             *dwv[ilhs + gd])^2)) 
        bd = where(wgt le 0.,nbd)
        if nbd eq irhs-ilhs+1 then begin
           wgt[*] = 1.          ; all portions
           if keyword_set(debug) then $
              print,'sdss_ewciv debug: uniform weighting (flux > conti)'
        endif else $
           wgt = wgt > 0.       ; no negative weights dragging things aside 
     endif else begin
        ;; Just do mean of region  (equal weight)
;        gd = where(sig[ilhs:irhs] ne 0.) ; why the repeating?!
;        if gd[0] ne -1 then $ 
;           wgt[gd] = 1. $       ; only good portion
;        else 
        wgt[*] = 1.             ; all portions
        if keyword_set(debug) then $
           print,'sdss_ewciv debug: uniform weighting (sig == 0)'

        ;; EW default
        qstrct.(ewtag)[ii] = !values.f_nan ; maybe these should be 0?
        qstrct.(sigewtag)[ii] = !values.f_infinity
     endelse 
     obswv[ii] = total(wave[ilhs:irhs]*wgt)/total(wgt)
     if keyword_set(debug) and finite(obswv[ii]) eq 0 then stop 
     if ztag ne -1 then begin
        qstrct.(ztag)[ii] = obswv[ii]/qstrct.(wrtag)[ii,cindx] - 1.
        qstrct.(sigztag)[ii] = sqrt(1./total(1./wgt))/qstrct.(wrtag)[ii,cindx]
     endif 

     ;; Re-evaluate
     if not keyword_set(snr_conv) and not keyword_set(keepwvlim) then begin
        ;; Edges: round 2
        mn = min(abs(wave-(obswv[ii]-boxw)),ilhs)
        mn = min(abs(wave-(obswv[ii]+boxw)),irhs)
        ;; Width
        qstrct.(wvlimtag)[ii,*] = [wave[ilhs],wave[irhs]]
     endif 
     
     ;; Rest EW and message
     if ztag ne -1 then begin
        qstrct.(ewtag)[ii] =  qstrct.(ewtag)[ii] / (1.+qstrct.(ztag)[ii])
        qstrct.(sigewtag)[ii] =  qstrct.(sigewtag)[ii] / (1.+qstrct.(ztag)[ii])

        ;; Set flag;  enables setting of lower and upper limit flags
        ;; (denotes badness)
        qstrct.ewflg[ii] = sdss_getlimflg() ; measured!
        
        if keyword_set(debug) then begin
           if n_elements(zabs) eq 1 then $
              print,'sdss_ewciv debug: changing z'+$
                    strtrim(round(qstrct.(wrtag)[ii,cindx]),2) $
                    + '= '+string(zabs,format='(f8.5)')+$
                    ' to ' + string(qstrct.(ztag)[ii],format='(f8.5)') $
           else $
              print,'sdss_ewciv debug: changing z'+$
                    strtrim(round(qstrct.(wrtag)[ii,cindx]),2) $
                    + '= '+string(zabs[ii],format='(f8.5)')+$
                    ' to ' + string(qstrct.(ztag)[ii],format='(f8.5)')
        endif 
        if (qstrct.wrest[ii]*(1+qstrct.(ztag)[ii]) lt qstrct.(wvlimtag)[ii,0] or $
            qstrct.wrest[ii]*(1+qstrct.(ztag)[ii]) gt qstrct.(wvlimtag)[ii,1]) and $
           keyword_set(debug) then stop,'sdss_ewciv() debug stop: centroid out of wavelength bounds'
     endif 

     if ncolmtag ne -1 then begin
        ;; Column density by AODM
        if irhs-ilhs lt 3 then begin
           ;; Check added for MgII 25 Jan 2013
           print,'sdss_ewciv(): npix < 3, growing just for X_AODM ',$
                 qstrct.qso_name,qstrct.(ztag)[0]
           if ilhs eq 0 then begin
              ilhs_aod = ilhs
              irhs_aod = irhs + 2
           endif else begin
              if irhs eq npix-1 then begin
                 ilhs_aod = ilhs - 2
                 irhs_aod = irhs
              endif else begin
                 ilhs_aod = ilhs-1 > 0
                 irhs_aod = irhs+1 < npix-1
              endelse 
           endelse 
        endif else begin
           ilhs_aod = ilhs
           irhs_aod = irhs
        endelse
        ;; Should this be flux0?
        x_aodm,wave[ilhs_aod:irhs_aod],flux[ilhs_aod:irhs_aod]/conti[ilhs_aod:irhs_aod],$
               sig[ilhs_aod:irhs_aod]/conti[ilhs_aod:irhs_aod],qstrct.(wrtag)[ii,cindx],$
               ncolm,signcolm,/log,flg_sat=flg_sat 
        qstrct.(ncolmtag)[ii] = ncolm
        qstrct.(signcolmtag)[ii] = signcolm

        ;; Set flags; enables setting of lower and upper limit flags
        ;; (denotes badness)
        if qstrct.(ncolmtag)[ii] gt 0. then $
           qstrct.ncolmflg[ii] = sdss_getlimflg() ; measurement
        if flg_sat ne 0 then $
           qstrct.ncolmflg[ii] = sdss_setlimflg(qstrct.ncolmflg[ii],/lower) ; sat!
        if qstrct.(signcolmtag)[ii] gt 1./(3*alog(10.)) then $              ; check 3sig
           qstrct.ncolmflg[ii] = sdss_setlimflg(qstrct.ncolmflg[ii],/upper) ; bad measure
     endif else flg_sat = 0

     ;; Plot
     if (keyword_set(plot) or keyword_set(debug)) and $
        not keyword_set(keepwvlim) then begin
        ;; Centroid
        if ztag ne -1 then $
           oplot, (1+[qstrct.(ztag)[ii],qstrct.(ztag)[ii]]) * $
                  qstrct.(wrtag)[ii,cindx], [-1e5,1e5], color=clr.cyan, $
                  linest=2, thick=2 $
        else oplot, [1.,1.]*qstrct.(wrtag)[ii,cindx], [-1e5,1e5], color=clr.cyan, $
                    linest=2, thick=2
        ;; Bounds
        oplot, [qstrct.(wvlimtag)[ii,0],qstrct.(wvlimtag)[ii,0]], $
               [-1e5,1e5], color=clr.green, $ ;GO
               linest=2, thick=4
        oplot, [qstrct.(wvlimtag)[ii,1],qstrct.(wvlimtag)[ii,1]], $
               [-1e5,1e5], color=clr.red, $ ;STOP
               linest=2, thick=4

;        oplot, [wave[ilhs_test],wave[ilhs_test]],$
;               [-1e5,1e5], color=clr.green, $ ;GO
;               linest=1, thick=3
;        oplot, [wave[irhs_test],wave[irhs_test]], $
;               [-1e5,1e5], color=clr.red, $ ;STOP
;               linest=1, thick=3

        print, qstrct.(wrtag)[ii,cindx], qstrct.(ewtag)[ii]
        if keyword_set(generic) then stop $
        else if ii EQ nlin-1 then begin
           mask = 1
           ;; Since many lines may make the initial velocity cut, have to do
           ;; a sanity check because lots of narrow lines are getting
           ;; chopped up. 
           ;; First, inspect EW because don't accept junk here; true
           ;; broad guys should result in positive EW in both lines
           if qstrct.ew_orig[0] lt 0. or qstrct.ew_orig[1] lt 0. then $
              mask = 0 $
           else begin
              ;; Do shape checks, namely how the flux looks around key portions
              mn = min(wave-qstrct.wvlim_orig[0,1],ihiI,/abs)
              mn = min(wave-qstrct.wvlim_orig[1,0],iloII,/abs)
              mn = min(wave-dblt.wvI*(1.+qstrct.zabs_orig[0]),iwvI,/abs)
              mn = min(wave-dblt.wvII*(1.+qstrct.zabs_orig[1]),iwvII,/abs)
              idx = [-1,0,1]    ; going to take mean of range
              mean_fxI = mean(flux[iwvI+idx])
              mean_fxII = mean(flux[iwvII+idx])
              ;; If the mid-point is lower than either other flux-weighted centroid,
              ;; then it's one line that's been divided
              if mean(flux[ihiI+idx]) lt mean_fxI and $
                 mean(flux[iloII+idx]) lt mean_fxII then $
                    mask = 0 $
              else begin
                 ;; Check strength at center of each bounds, secondary check on above
                 mn = min(wave-mean(qstrct.wvlim_orig[0,*]),icI,/abs)
                 mn = min(wave-mean(qstrct.wvlim_orig[1,*]),icII,/abs)
                 if mean(flux[icI+idx]) gt mean_fxI and $
                    mean(flux[icII+idx]) gt mean_fxII then $
                       mask = 0 $
                 else begin
                    if (abs(ihiI-iwvI) le 1 and abs(iloII-iwvII) le 1) or $
                       (ihiI eq iwvI or iloII eq iwvII) then $
                          mask = 0 $
                    else begin
                       ;; Don't want min flux in range to be at boundary
                       mn = min(wave-qstrct.wvlim_orig[0,0],iloI,/abs)
                       mn = min(wave-qstrct.wvlim_orig[1,1],ihiII,/abs)
                       mn = min(flux[iloI:ihiII],imn)
                       if abs(imn-ihiI) le 1 or abs(imn-iloII) le 1 then $
                          mask = 0
                    endelse     ; absolute min
                 endelse        ; center flux > mean fxI/II

                 ;; If well-aligned in boxes, keep (actually pretty good for
                 ;; picking off MgII); undoes ill effects of previous
                 if abs(icI-iwvI) le 1 and abs(icII-iwvII) le 1 then $
                    mask = 1
              endelse           ; boundary flux min
           endelse              ; EW test
           
           xyouts,0.95,0.1,'Mask = '+strtrim(mask,2),color=clr.purple,$
                  charsize=3.,/normal,alignment=1
           stop
        endif                   ; generic == 0
     endif                      ; /plot or /debug
     
  endfor                        ; loop ii=nlin
  
  ;; don't pass back out wrong
  if keyword_set(generic) then zabs = qstrct.(wrtag)[0:nlin-1,cindx] 

  flux = flux0                  ; don't send back a messed up thing

  ;; restore
  return, qstrct

end 


