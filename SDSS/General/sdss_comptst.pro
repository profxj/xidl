;+ 
; NAME:
; sdss_comptst
;    Version 1.1
;
; PURPOSE:
;  Monte-carlo routine to check the completeness of our automoated
;  algorithm
;   
; CALLING SEQUENCE:
;   sdss_comptst, fil, outfil, NTRIAL=, /CHK, /DEBUG
;
; INPUTS:
;  fil -- SDSS spectrum (without a DLA present)
;
; RETURNS:
;
; OUTPUTS:
;  outfil -- FITS file with the qalstrct
;
; OPTIONAL KEYWORDS:
;  NTRIAL=  -- Number of simulations to run [default: 1000L]
;  /CHK  -- Plot the spectrum
;  /DEBUG -- 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_comptst, '/data5/SDSS/DR1_QSO/spectro/1d/0463/0463-51312-122.fit',
;    'monte_snr10.fits', NTRIAL=10000L
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   2003 Written by SHF
;-
;------------------------------------------------------------------------------
pro sdss_comptst, fil, outfil, NTRIAL=ntrial, CHK=chk, DEBUG=debug

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'sdss_comptst, fil, outfil, NTRIAL=, /DEBUG, /CHK [v1.1]'
    return
  endif 

 qalstrct = {qalcharstrct2}
parse_sdss,fil, flux, wave, conti, SIG=sig, IVAR=ivar, ZQSO=zqso, $
                   NPIX=npix, HEAD=head, CFIL=cfil, CDIR=cdir

a=flux
svsig = sig
s=-1233

  if not keyword_set( NTRIAL ) then NTRIAL=1000L


lines = {abslinstrct}
for i=0,ntrial-1 do begin
    
    NH1=20+randomu(s)
    qalstrct.NH1[i]=NH1
    
    zabs=2.9-.7*randomu(s)
    qalstrct.z_abs[i]=zabs
    
    lines[0] = x_setline(1215.6701d)
    lines.zabs = zabs
    print, i, ' NHI', NH1, ' Zabs',zabs
    lines[0].N = NH1
    lines[0].b = 10.
    
    waveem=1215.6701*(1+zqso)
    emi=where(abs(wave-waveem) lt 0.5)
    
    flux=a
    sig = svsig
    
    abswave=(1+zabs)*1215.6701
    
    dwave=abswave*(.00333)
    
    window=where(wave gt (abswave-dwave) and wave lt (abswave+dwave),nwin)
    noise=fltarr(emi+21)
    noise[window]=randomn(s,nwin)*sig[window]
    flux[0:emi+20]=a[0:emi+20]*(x_voigt(wave[0:emi+20],lines,FWHM=2.3))$
      +noise
    
    if keyword_set( DEBUG ) then stop
                                ;x_splot, wave, flux, /block
    sdss_dla_2, STRCT=qalstrct, ZEM=zqso, flux, wave, sig, i, zabs
    
    if keyword_set(CHK) then begin
        if qalstrct.ndla1[i] EQ 0 and NH1 GT 20.5 then begin
            x_splot, wave, flux, /block
            stop
        endif
    endif
    
    
endfor

mwrfits, qalstrct, outfil, /create

end


