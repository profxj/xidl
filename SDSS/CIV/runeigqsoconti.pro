;+
; NAME:
;   RUNEIGQSOCONTI
;     Version 2.0
;
; AUTHOR:
;   Melodie M. Kao
;   Massachusetts Institute of Technology
;   Kavli Institute for Astrophysics and Space Research
;   77 Massachusetts Avenue, Building 37-287, Cambridge, MA 02139 
;   melodie.kao@alum.mit.edu
;   mkao@caltech.edu
;
; PURPOSE:
;   Wrapper script to show how to run eigqsoconti.pro.  Uses SDSS
;   spectra--user will have to modify it to suit your own spectra.
;
;
; CALLING SEQUENCE:
;    runeigqsoconti, /debug
;
; DESCRIPTION: 
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   A FIT file containing an nx3 array for each spectrum input, where
;   n=number of pixels in spectrum.  
;   output[*,0] = eigencontinuum
;   output[*,1] = mask array (0 = used in fit; 1 = not used in fit)
;   output[*,2] = eigenerror (ONLY error from eigenfitting)
;
; OPTIONAL KEYWORDS:
;   sdsssum= -- DR7_QSO catalog-formatted FITS file to use 
;   eigspecfil= -- eigenspectra to use
;   /debug    -- plots fitting information.  IF turned on, then output
;                files are not written. (If debug > 1,
;                eigqsoconti's debug accessed too.)
;   /clobber -- overwrite files
;   istrt= -- index at which to start
;   processor= -- two element array of [ith, total] processors for a
;                 parallel run with sdss_runparallelsrch.sh script,
;                 where 1 <= ith <= (total <= 4). Overrides istrt=.
;   /tcpu -- print information about CPU time used
;   /help -- print Syntax
;   /fit_full -- extrapolate eigenconti to full spectrum (definitely
;                want /debug to inspect; largely crap)
;   _extra= -- passed to subroutines
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: 
;  
; EXAMPLES:
; 
; PROCEDURES/FUNCTIONS CALLED:
;   eigenrecstrct.pro                
;   eigenrecstrct_mij.pro            
;   getcolor
;   x_splot
;
; MODIFICATION HISTORY:
;   4 JUL 2011  Written by M. Kao
;  28 Jul 2011  Major reorg and overhaul, KLC
;   6 Jun 2017  Enable extrapolation to larger wavelength range, KLC
;
;-
; Copyright (C) 2011, Melodie Kao
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
;------------------------------------------------------------------------------
@eigqsoconti                    ; must compile since it's a function


pro runeigqsoconti, sdsssum=sdsssum, debug=debug, eigspecfil=eigspecfil, $
                    clobber=clobber, processor=processor, list=list, fits=fits, $
                    istrt=istrt, help=help, tcpu=tcpu, fit_full=fit_full, $
                    _extra=extra

  if keyword_set(help) then begin
     print,'Syntax - runeigqsoconti [sdsssum=, debug=, eigspecfil=, /clobber, processor='
     print,'                         /list or /fits, istrt=, /help, /tcpu, /fit_full, _extra=]'
     return
  endif 

  ;; Figure out where all the time goes!
  profiler, /reset
  profiler                      ; profile all user calls
  profiler, /system             ; profile all system calls

  sdssdir = sdss_getsdssdir() ; $SDSSPATH/$SDSSDR

  ;; Driving structure
  if keyword_set(sdsssum) then begin
     if (size(sdsssum,/type))[0] eq 8 then sdsstab = sdsssum $ ; already structure
     else sdsstab = xmrdfits(sdsssum,1,/silent) ; read in
  endif else $
     sdsstab = sdss_getqsostrct() ; default
  nqso = (size(sdsstab,/dim))[0]    ; KLC: I learned size() is faster than n_elements()

  if keyword_set(processor) then begin
     ;; Divide out jobs
     sub = sdss_calcparalleljob(sdsstab, processor)
     istrt = sub[0]
     nqso = sub[1] + 1
     if keyword_set(debug) then $
        print,'runeigqsoconti debug: multi-processor run for just ',istrt,nfil

     ;; Force single-thread
     save_cpu = !cpu
     cpu, tpool_nthreads=1
  endif 
  if not keyword_set(istrt) then istrt = 0L ; let processor= supercede


  ;; Defaults
  ;; Only allow redward of Lya to end of SDSS range
  wavebracket = [1230.0D, 10000.0D] ; match convention in sdss_fndlin_fitspline()
  ;; Defaults passed to eigqsoconti.pro
  if not keyword_set(eigspecfil) then $
     eigenArr = xmrdfits(getenv('SDSSPATH')+'/eigenspectra/eigSpec_qso_all.fit',$
                           0,/silent) $ ; save time and just pass this around
  else begin
     ;; Read input
     if keyword_set(list) then $
        eigenArr = eigenrecstrct_mij_rdfil(eigspecfil) $ ; Construct array
     else begin 
        if keyword_set(fits) then $
           eigenArr = xmrdfits(eigspecfil,0,/silent) $ ; Read fits file
        else begin
           ;; Figure it out
           if size(eigspecfil,/n_dimen) eq 2 then $ 
              eigenArr = eigspecfil $ ; assume input is the desired array
           else begin
              prs = strsplit(eigspecfil[0],'.',/extract,count=nprs)
              if stregex(prs[nprs-1],'fit',/boolean) then $
                 eigenArr = xmrdfits(eigspecfil[0],0,/silent) $ ; Read fits file
              else eigenArr = eigenrecstrct_mij_rdfil(eigspecfil) ; assume need to construct
           endelse                                             ; ndimen != 2
        endelse                                                ; /fits not set
     endelse                                                   ; /list not set
  endelse                                                      ; not using default


  if keyword_set(debug) then $
     if debug gt 1 then debug_sub = 1 $ ; for eigqsoconti
     else debug_sub = 0

  ;; Use function and set up names consistently
  specfile = sdss_getname(sdsstab, dir=spectrodir, root=qso_name)
  specfile = spectrodir + specfile 
  outfile = sdss_getname(sdsstab, /eig, extrap=fit_full, dir=contidir)
  outfile = contidir+outfile

  ;; Timing
  tstart = systime(/seconds)
  tlast = tstart

  FOR i=istrt, nqso-1 DO BEGIN
    
     test = file_search(sdssdir+outfile[i]+'*',count=ntest)
     if ntest ne 0 and not keyword_set(clobber) then begin
        ;; Skipping
        print,'runeigqsoconti: file exits; continuing: ',outfile[i]
        continue
     endif

     ;; Files account keeping
     IF NOT KEYWORD_SET(debug) THEN BEGIN
        test = file_search(sdssdir+contidir[i],count=ntest)
        IF ntest EQ 0 THEN                                                   $
           spawn,'mkdir -p '+sdssdir+contidir[i]
     ENDIF
     
     
     ;; Input spectrum information (KLC: use function already written)
     parse_sdss, sdssdir+specfile[i], flux, waveObs, sig=fluxerr, npix=npix, $
                 head=hdreig
     zQSO     = sdsstab[i].z    ; don't trust header
     waveEmit = waveObs/(1.0+zQSO)
     
     ;; Spectrum portions to actually be fitten and written to file
     finalArr = dblarr(npix, 3)
     
     ;; Remove Lya portion of spectrum.
     keep = WHERE(waveEmit GE wavebracket[0] AND waveEmit LE wavebracket[1], nkeep, $
                  complement=bad, ncomplement=nbad)
     IF keyword_set(fit_full) THEN BEGIN
        waveEmit_full = waveEmit  ; trimmed later
        flux_full = flux
        fluxerr_full = fluxerr
     ENDIF ELSE  waveEmit_full = 0 ; for wv_full keyword
     waveEmit = waveEmit[keep]
     flux     = flux[keep]
     fluxerr  = fluxerr[keep]
     IF nbad NE 0 THEN BEGIN
        finalArr[bad,0]  = !VALUES.F_NAN    ; For arrays to be outputted, 
        finalArr[bad,1]  = 1                ; masked out
        finalArr[bad,2]  = !VALUES.F_NAN    ; make the Lya region NaN rather
     ENDIF
     
     ;; ***********************************************************************
     ;; ******************** BEGIN eigenfitting *******************************
     ;; ******                                                          *******
     
     eigArr = eigqsoconti(waveEmit,flux,fluxerr,eigenArr, $
                          FITMASK=fitmask, finalmask=finalmask, $
                          NITER=niter, FAIL=fail, STATUS=status, CHI_SQR=chi_sqr, $
                          STAT_CTOL=stat_ctol, debug=debug_sub, header=hdreig,$
                          wv_full=waveEmit_full, idx_sub=keep, $
                          _extra=extra) ; /silent


     ;; ******                                                            *****
     ;; ******************** END eigenfitting *********************************
     ;; ***********************************************************************

     ;; previous nkeep ne 0 check was unnecessary; would have crashed earlier
     IF keyword_set(fit_full) THEN BEGIN
        keep = lindgen(npix)    ; all

        if keyword_set(debug) then begin
           nkeep = npix
           waveEmit = waveEmit_full ; for plotting
           flux = flux_full
           fluxerr = fluxerr_full
        endif
     ENDIF
     
     finalArr[keep,0] = eigArr[*,0]
     finalArr[keep,1] = finalmask ; good place to store it
     finalArr[keep,2] = eigArr[*,1]
     
     
     ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ;; %%%%%%%%%%%%%%%%%%%% BEGIN debug plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ;; %%%%%% 
     redchisq = chi_sqr[0]/chi_sqr[1]
     IF KEYWORD_SET(debug) THEN BEGIN
        clr = getcolor(/load) 
        masked = replicate(!VALUES.F_NAN,nkeep)
        bad = where(finalArr[keep,1] eq 1,nbad)
        if nbad ne 0 then masked[bad] = flux[bad]
        fluxerr = sqrt(fluxerr^2 + finalArr[keep,2]^2)
        title =  'Final Continuum for '+qso_name[i]+$
                 ' (#'+ STRTRIM( i,2 ) +       $
                 string(redchisq,chi_sqr[2],$
                        format='(") CHISQR=",f6.3," (P=",f6.3,")")')
        x_splot, waveEmit,                                                $
                 flux,               psym1=10,                            $
                 ytwo=fluxerr, psym2=10,        $
                 ythr=finalArr[keep,0], psym3=-3, $
                 yfou=masked, psym4=6, color4=clr.magenta,               $
                 ymnx=[MIN(flux, /NAN, max=mx), mx],                      $
                 TITLE= title, /BLOCK,                                    $
                 LGND=['flux','flux+conti error','final conti','masked']   
     END
     ;; %%%%%%                                                            %%%%%
     ;; %%%%%%%%%%%%%%%%%%%% END debug plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     print,''
     mwrfits, finalArr, sdssdir+outfile[i], hdreig, /CREATE, /SILENT
     mwrfits, stat_ctol, sdssdir+outfile[i], /SILENT
     spawn,'gzip -f '+sdssdir+outfile[i]
     print,'runeigqsoconti: created ',outfile[i]

     print,qso_name[i]+' SUMMARY:'
     print,'Exit Status = '+status,'; Fail = ',strtrim(fail,2)
     print,'# iter = '+strtrim(niter),$
           '; ctol = ',string(stat_ctol[niter-1],format='(f9.6)')
     print,'Reduced chi^2 = ',string(redchisq,format='(f5.2)'),$
           '; Prob = ',string(chi_sqr[2],format='(f5.3)')
     print,''

     if keyword_set(tcpu) and ((i+1) mod 1000L) eq 0 then begin
        tt = systime(/seconds)
        dt = (tt-tlast)         ; seconds
        print,'runeigqsoconti: Elapsed time (m) = ',dt/60.
        print,'runeigqsoconti: Average time per i (s) = ',dt/1000L
        tlast = tt
     endif 

  ENDFOR                        ; loop i=nsdss

  ;; Final messages
  print, 'runeigqsoconti: All done!'
  tlast = systime(/seconds)
  dt = tlast - tstart
  print,'runeigqsoconti: Elapsed time for '+strtrim(nqso,2)+$
        ' QSOs (m) = ',dt/60.
  print,'runeigqsoconti: Average time per QSO (s) = ',dt/(nqso)

  profiler,/report,filename='runeigqsoconti.profile'
  print,'runeigqsoconti: IDL profiler report saved to runeigqsoconti.profile'

  ;; Revert back to desired thread pool
  if keyword_set(processor) then $
     cpu, restore=save_cpu

end
