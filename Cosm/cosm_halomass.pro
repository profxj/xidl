;+
;
; Function that returns halo mass information using the ST
; approximation and the formalism summarized in Dekel&Birnboim 2006
; appendix. This is a idl version of the code written by Dekel. 
; This calls a c++ program.
;
; INPUT:
; massmin   --> log10 min halo mass (H0 is applied)
; massmax   --> log10 max halo mass (H0 is applied)
; z      --> redshift used 
; EXTRA  --> all the cosmology supported by common
; H0     --> Redefine H0
; prec   --> precision used in the integration of mass function
;            (default 0.008)
; savehm --> when ndens is set, provide the name of a file where the
;            halo mass function is written
;
; OUTPUTS:
; 
; Return the cumulative halo mass in #/Mpc^3 in the interval set (cosmology applied)
;
; ndens  --> if set, compute the halo mass function in bins of 0.2 dex
; emstar --> the characteristic mass (log10 Msun)
; 
; KEYWORD:
; silent --> Suppress messages
; smart  --> Truncate the halo mass function for less than 10^-20
;-


function cosm_halomass, massmin, massmax, z, _EXTRA=extra , ndens=ndens, $
                   emstar=emstar, silent=silent, H0=H0, prec=prec, savehm=savehm, $
                   smart=smart

;;default
  if ~keyword_set(prec) then prec=0.008
  if ~keyword_set(savehm) then begin
     if(z lt 10) then savehm='halomass_z'+string(z,format="(F5.3)")+'.dat'
     if(z ge 10) then savehm='halomass_z'+string(z,format="(F5.2)")+'.dat'
  endif
 
;;start the cosmology
  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, cosm_L, cosm_r, sigma_8
  cosm_common, H0=h0, _EXTRA=extra, SILENT=silent, /W05MAP
  ;cosm_common, H0=h0, _EXTRA=extra, SILENT=silent, /VANILLA
  littleh=cosm_h/100.
  

;;call the c++ program passing the cosmology
  spawn, [getenv("XIDL_DIR")+"/Cosm/cosm_dekelmiki_psn",string(massmin),string(massmax),string(z),string(cosm_dm),string(littleh),string(sigma_8),string(prec)], output, /noshell
        
  
;;parse output
  p1=strpos(output,"W")
  ends=strpos(output,"Z")
  
  cumhalo=1D*strmid(output,0,p1)
  emstar=alog10(1D11*strmid(output,p1+1,ends))
  

;;halo mass function
  if keyword_set(ndens) then begin
     binmass=0.2
     nbin_mass=floor((massmax-massmin)/binmass)
     
     ;;set the halo mass function
     hmfunct=dblarr(nbin_mass)
     halomass=dblarr(nbin_mass)
     hmfunct[0]=cumhalo
     halomass[0]=massmin+0.5*binmass

     for i=1, nbin_mass-1 do begin
        halomass[i]=massmin+(i+0.5)*binmass
    
        minlimit=massmin+i*binmass
        ;;compute function
        spawn, [getenv("XIDL_DIR")+"/Cosm/cosm_dekelmiki_psn",string(minlimit),string(massmax),string(z),string(cosm_dm),string(littleh),string(sigma_8),string(prec)], output, /noshell
        
        ;;parse output
        p1=strpos(output,"W")
        hmfunct[i]=1D*strmid(output,0,p1)
        if (keyword_set(smart) and hmfunct[i] lt 1D-20) then break

     endfor
     
     ;;make the distribution
     hmfunct=(hmfunct-shift(hmfunct,-1))/binmass
     
     ;;cut the last bin which is not right
     hmfunct=hmfunct[0:nbin_mass-2]
     halomass=halomass[0:nbin_mass-2]
     
     ;;print the halo mass in file
     forprint, halomass, hmfunct, textout=savehm, silent=silent, $
               comment='Log10 Mass (10^11 Msun); halo #/Mpc^3/bin with bin '+string(binmass,'(f4.2)')+' dex'

  endif
  

return, cumhalo




end
