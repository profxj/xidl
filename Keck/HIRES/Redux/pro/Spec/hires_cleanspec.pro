pro hires_cleanspec, flux, sig, wav, wavmin = wavmin, wavmax = wavmax, $
                  fluxmin=fluxmin,fluxmax=fluxmax

   
 if not keyword_set(wavmin) then wavmin=3000.0
 if not keyword_set(wavmax) then wavmax=10000.0
 if not keyword_set(fluxmin) then fluxmin=-2.0
 if not keyword_set(fluxmax) then fluxmax=2.0

                                ;strip out NaNs, infs
 inf=where(finite(flux,/infinity),ni)
 if(ni ne -1) then begin
    flux[inf]=0.0
    sig[inf]=-1.0
 endif
 
 nan=where(finite(flux,/nan),nn)
 if(nn ne -1) then begin
    flux[nan]=0.0
    sig[nan]=-1.0
  endif
 
  ;nuke all wave < wavemin > wavemaxflux values
 wmn=where(wav lt wavmin,nwmn)
 if(nwmn gt 0) then begin
    flux[wmn]=0.0
    sig[wmn]=-1.0
 endif
 
 wmx=where(wav gt wavmax,nwmx)
 if(nwmx gt 0) then begin
    flux[wmx]=0.0
    sig[wmx]=-1.0
 endif

  ;nuke all flux > fluxmax flux < fluxmin values
 fm=where(flux lt fluxmin,nfm)
 if(nfm gt 0) then begin
    flux[fm]=0.0
    sig[fm]=-1.0
 endif
     
 fmx=where(flux gt fluxmax,nfmx)
 if(nfmx gt 0) then begin
    flux[fmx]=0.0
    sig[fmx]=-1.0
 endif

 
;  mwrfits,flux,fluxfil,head,/create
;  mwrfits,sig,errfil,errhead,/create



 return
end
