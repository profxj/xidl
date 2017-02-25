;+ 
; NAME:
; sdss_dla
;    Version 1.1
;
; PURPOSE:
;  Main routine to search for DLA in SDSS spectra (or ESI).  This
;  is the key driver of the N(HI) search for DLA.
;   
; CALLING SEQUENCE:
;  sdss_dla, n2=, ESI=, SDSS=, WIN_REST=, /PFLUX, /SMOOTH, /NORM, 
;     SIGV=, STRCT=, ZEM=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;  STRCT=  -- QAL structure created by this program
;
; OPTIONAL KEYWORDS:
;  n2= -- Sigma value above which a pix is bad in the DLA profile [default: 5]
;  SIGV= -- Sigma value to add to flux of all pixels
;  SDSS= -- Name of SDSS file to examine
;  /PFLUX  -- Plot the spectrum
;  /SMOOTH -- Smooth the spectrum by 5 pixels
;  /NORM   -- The flux is normalized
;  WIN_REST -- Width of DLA window (angstroms)  [default: 6]
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
;   Aug2006 -- Modified for PDLA
;-
;------------------------------------------------------------------------------
pro sdss_dla, n2=n2, ESI=esi, SDSS=sdss, WIN_REST=win_rest, $
         PFLUX=pflux, SMOOTH=smooth, NORM=norm, SIGV=sigv, $
         STRCT=qalstrct, MAG=mag

;print, 'begin sdss_dla.pro ...'
qalstrct = {qalcharstrct}
 
if not keyword_set( n2 ) then n2= 5.
if not keyword_set( WIN_REST ) then win_rest = 6.
if not keyword_set(POFF) then poff = 0.1 
;if not keyword_set( SIGV ) then sigv = 5.e-17

;;read in info from SDSS file

if keyword_set(SDSS) then begin
    qalstrct.file_name = sdss
    d1= xmrdfits(SDSS, 0, head, /silent)
    iwave = sxpar(head, 'CRVAL1')
    name = sxpar(head, 'NAME')
    fiberid=strtrim(sxpar(head,'FIBERID'),2)
    qalstrct.qso_name=name+'-'+fiberid
    if not keyword_set(MAG) then mag = sxpar(head, 'MAG_R')
    qalstrct.qso_mag=mag
    flux= d1[*,0]
    pixels= n_elements(flux)
    sig= d1[*,2]
    sdss_wave, wave, pixels=pixels, iwave=iwave
    zem = sxpar(head, 'Z')
    z_err= sxpar(head, 'z_err')
    qalstrct.Z_qso= zem
    RA = sxpar(head, 'RAOBJ')
    qalstrct.RA= RA
    DEC = sxpar(head, 'DECOBJ')
    qalstrct.DEC= DEC
endif

if keyword_set(ESI) then begin
    dla_esitosdss, esi, wave, flux, sig, NORM=norm
endif

pixels= n_elements(wave)

if keyword_set(NORM) then begin
    flux =flux + ( (1/15.)*randomn(s, pixels))
    sig= replicate(1./15.,n_elements(flux))
endif 

if keyword_set(SIGV) then begin
  x_splot, wave, flux, ytwo=sig, /block
    sig = replicate(sigv, pixels)
    flux = flux + (sig * randomn(s,pixels))
endif

index= lindgen( pixels)


if keyword_set( SMOOTH ) then begin
    flux = smooth(flux,5)
    sig  = sig / sqrt(5)
endif

if keyword_set( PFLUX ) then begin
    x_splot, wave, flux, ytwo=sig, /block
                                ; stop 
endif

;;mask redward of zem

o= where(wave GE (zem+1)*1215.6701)
;print, 'z_header:', zem
;print, 'z_err:', z_err
;x_splot, wave, flux
;stop

;;if emission peak is too blue
if o[0] LE 200 then begin
    if o[0] LE 50 then begin
        print, 'sdss_dla.pro: Starting wave too low, returning'   
        return
    endif else begin
        n1=fltarr(pixels)
        ;; Avoid sig = 0
        garr = o[0] + 200L + lindgen(151)
        gd = where(sig[garr] GT 0., ngd)
;        gd = where(sig[o[0]+200:o[0]+350] GT 0., ngd)
        if ngd EQ 0 then n1 = 1.5 else begin
            snr=flux[garr[gd]]/sig[garr[gd]]
            n1=1. > (median(snr))/2.5 < 2.
            qalstrct.SNR=median(snr)
        endelse
    endelse
endif

;;SNR (n1)
if not keyword_set( n1 ) then begin
    if o[0] GT 200 then begin        
        n1=fltarr(pixels)
        garr = o[0] - 200L + lindgen(151)
        gd = where(sig[garr] GT 0.,ngd)
;        gd = where(sig[(o[0]-200)>0:(o[0]-50)>0] GT 0.,ngd)
        if ngd EQ 0 then begin
            garr = o[0] + 200L + lindgen(151)
            gd = where(sig[garr] GT 0., ngd)
            if ngd EQ 0 then n1 = 1.5 else begin
                snr=flux[garr[gd]]/sig[garr[gd]]
                n1=1. > (median(snr))/2.5 < 2.
                qalstrct.SNR=median(snr)
            endelse
        endif else begin
            snr=flux[garr[gd]]/sig[garr[gd]]
            n1=1. > (median(snr))/2.5 < 2.
            qalstrct.SNR=median(snr)
        endelse
    endif
endif

;; Zero out data outside the forest
p= where(wave GE (zem+1+POFF)*1215.6701)
flux[p]=0.
sig[p]=0.

;if qalstrct.SNR ne 0. then print, 'qalstrct.SNR:', qalstrct.SNR 

zarr= (wave/1215.6701d)-1

window= win_rest*(1+zarr)

nindx = n_elements(zarr)

score= fltarr(nindx)
bad= fltarr(nindx)

;s=fltarr(4000)
s=where(flux GT 0. AND sig GT 0.)

;stop

if s[0] eq -1 then return

;;calc scores, badscores
for index= s[0], p[0] do begin
    
    w= window[index]
        
    low= wave[index] - w/2      ; Angstorms
    high= wave[index] + w/2     ; Angstroms
    
    pix= where( wave gt low and wave lt high, npix)
    
    gpix= where( flux[pix] lt n1*sig[pix], ngpix)
        
    bpix= where( flux[pix] gt n2*sig[pix], nbpix)
        
    score[index]= float(ngpix)/ float(npix)
    bad[index]= float(nbpix)/ float(npix)
        
endfor
    
;x_splot, zarr, score, /block, ytwo=bad
    
    
    
medscore= median(score, 10)
    
;x_splot, wave, medscore, /block
    
    
x= fltarr(3900L)
    
;if fiberid EQ 500 then stop
;;lls
;if median(sig[s[0]:s[0]+30]/flux[s[0]:s[0]+30]) GT 1. then s[0]=s[0]+100
if median(medscore[s[0]:s[0]+20]) gt .6 then begin

    for index= s[0], pixels-1 do begin
        
        zindex= index + lindgen(16)
        
        x= where( medscore[zindex] lt 0.5, count )
        
        if (count eq 15) then break
        
    endfor 
    
    if wave[index] LT 920.d*(1+zem) then begin

        zll= (wave[index]/912.)-1
        
        print, 'z l.limit:', zll
        
        qalstrct.flg_LLS= 1
        qalstrct.LLS_zabs= zll
                                ;qalstrct.LLS_score= ?
    endif
    
endif

;;start looking for dla
flux[0:s[0]]=0.
sig[0:s[0]]=0.

;; Zero out start of sig array as necessary -- New!!
if median(flux[0:9]/sig[0:9]) LT 4. then flux[0:9] = 0.

if qalstrct.flg_LLS EQ 1 then $
  start=where(median(flux/sig, 20) GT 4. AND sig NE 0. AND wave GT (zll+1)*912.) $
else start=where(median(flux/sig, 20) GT 4. AND sig NE 0.)

if start[0] EQ -1 then begin
    
    print, 'poor data!'
    qalstrct.nDLA1=0
    qalstrct.DLA_zabs1=-1.
                                ;qalstrct.DLA_z=-1.
    qalstrct.DLA_score1=-1.
    qalstrct.start_wave=-1.
endif else begin
    
    qalstrct.start_wave= wave[start[0]]
    
    score[0:start[0]]=0
    
    dv=fltarr(3900)
    gi=where(score gt 0.6)
    
;;if no dla  
    
    if gi[0] EQ -1 then begin
        print, 'no DLAs!!'
        
        qalstrct.nDLA1=0.
        qalstrct.DLA_zabs1=0.
                                ;qalstrct.DLA_z=0.
        qalstrct.DLA_score1=0.
        
    endif else begin
        
;;create gaps, find dla in each gap, give z, score for each
        
        gwave= wave[gi]
        dgwave= gwave-shift(gwave,1)
        dgwave[0]= 0.
        dv=(dgwave/gwave)*3.E5
        gapindex= where(dv GT 2000., ngindex) 
        
        
        
        if ngindex EQ 0 then begin
            z= median(zarr[gi])
            dlascore=max(score[gi])
        endif else begin
            
            gapwave= wave[gi[gapindex]]
            goodz=zarr[gi[gapindex]]
            z=fltarr(ngindex+1)
            dlascore=fltarr(ngindex+1)
            
            if ngindex EQ 1 then begin
                window1=where(gwave LT gapwave[0])
                window2=where(gwave GT gapwave[0])
                z[0]=median(zarr[gi[window1]])
                dlascore[0]=max(score[gi[window1]])
                if window2[0] NE -1 then begin
                    z[1]=median(zarr[gi[window2]])
                    dlascore[1]=max(score[gi[window2]])
                endif
            endif else begin
                
                for j=0, ngindex-1 do begin
                    
                    if j EQ 0 then begin 
                        window=where(gwave LT gapwave[0])
                        z[j]=median(zarr[gi[window]])
                        dlascore[j]=max(score[gi[window]])
                    endif
                    
                    if j GT 0 AND j LT ngindex-1 then begin
                        window=where(gwave GE gapwave[j-1] $
                                     AND gwave LT gapwave[j])
                        z[j]=median(zarr[gi[window]])
                        dlascore[j]=max(score[gi[window]])
                    endif
                    
                    if j EQ ngindex-1 then begin
                        window1=where(gwave LT gapwave[j] AND gwave GE gapwave[j-1])
                        z[ngindex-1]=median(zarr[gi[window1]])
                        dlascore[ngindex-1]=max(score[gi[window1]])
                        window2=where(gwave GE gapwave[ngindex-1])
                        z[ngindex]=median(zarr[gi[window2]])
                        dlascore[ngindex]=max(score[gi[window2]])
                    endif
                    
                endfor
            endelse
        endelse
        
;;fill strct 
;stop
        qalstrct.nDLA1= ngindex+1
        
        print, 'zdla:', z
        print, 'dlascore:', dlascore
        qalstrct.DLA_zabs1[0:ngindex]= z
                                ;qalstrct.DLA_z[0:ngindex]= z
        qalstrct.DLA_score1[0:ngindex]= dlascore 
    endelse
endelse

end








