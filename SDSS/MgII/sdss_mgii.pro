;+ 
; NAME:
; sdss_mgii
;    Version 1.0
;
; PURPOSE:
;    Brute force algorithm which examines the absorption lines
;    detectedin a SDSS spectrum and searches for matches with MgII
;
; CALLING SEQUENCE:
;  sdss_fndciv, abswav
;
; INPUTS:
;   abswav - Array of observed wavelengths of detected absorption lines
;
; RETURNS:
;   
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------
;takes wavelength, wavelength-to-pixel, absorption wavelength arrays,
;passes back score array of matched lines for each z in array of z's, zarr

pro sdss_mgii, abswave, ZEM=zem, QSTRCT=strct


  strct={qalcharstrct}

;define score and z arrays

  nsearch = 50000L

  score= fltarr(nsearch)
  zarr= dblarr(nsearch) 
  sc=fltarr(nsearch)

  ;;define rest wavelengths of lines to be searched for
  rwave= [2600.1729d, 2796.3520d, 2803.5310d]
;  rwave= [2796.3520d, 2803.5310d]
  nrwave = n_elements(rwave)
  maxabswave=max(abswave)
  if maxabswave LE 8000. then maxout=8000. else maxout=9200. 

  ;;  start/end wave
  strct.start_wave = (1+zem)*1230.d > 3820.
  strct.DLA_quality[0] = maxout < (1+zem)*2803.5310d


;loop over z values from 0 to 5 in .001 increments, at each z
;calculate array of shifted wavelengths (abswavearr) corresponding to each rest wavelength

for j=0L, nsearch-1 do begin 
    zarr[j]= j*.0001d
    if zarr[j] LT 0.45 then continue
    if zarr[j] GT 2.4 then continue
    
    
    abswavearr= rwave*(1.+zarr[j])
    goodlines= where((abswavearr GT (1+zem)*1230.d) AND $
                     (abs(abswavearr-5579.d) GT 5.d) AND $
                     (abs(abswavearr-6302.d) GT 5.d) AND $
                     abswavearr LT maxout, nposs) 
    nposs = nposs < 2

    sc[j]=0
    
;at each z, search for a match of each shifted rest wavelength with a
;wavelength in the absorption wavelength array, within a certain
;tolerance, and count how many matches are found 
 

    for i=0, nrwave-1 do begin
        
        x= where(abs(abswavearr[i]- abswave) $
                 LT abswavearr[i]*(1./4000.), count) 
        ;;if any lines match, add 1 to score for each match
        if (count NE 0) then sc[j]= sc[j]+1
      
;        if zarr[j] EQ 2.0677d then stop

    endfor
   

;define score array as scores for each .001 increment in z
;define z arr as z values for each increment
    
    score[j]= sc[j]/float(nposs)
    
    ;if sc[j] EQ 2 and score[j] NE 1. then score[j] = 0
    if sc[j] LT 2 then score[j] = 0.
    
;   if sc[j] GE 3 AND score[j] GT 0.6 then begin
 ;      x_specplot, flux, fltarr(n_elements(flux)), wave=wave, inflg=4,$
                                ;       /block, zin=zarr[j], /qal
                                ;      stop
                                ;  endif

endfor

z=dblarr(100)
dlasc=fltarr(100)
hits=fltarr(100)

iscore= where(score GT 0.9999, ngscore)

if iscore[0] NE -1 then begin

    ghits=sc[iscore]
    gscore=score[iscore]
    gz= zarr[iscore]
    dgz= gz-shift(gz,1)
    dgz[0]= 0.
                                ;gwave=abswave[iscore]
    gapzindex= where(dgz GT .01, ngap)
    
    
    if ngap EQ 0 then begin
        z=max(gz)
        dlasc=max(gscore)
        hits=max(ghits)
    endif else begin
        
        gapz= zarr[iscore[gapzindex]]
             
        if ngap EQ 1 then begin
            window1=where(gz LT gapz[0])
            window2=where(gz GE gapz[0])
            z[0]=max(zarr[iscore[window1]])
            dlasc[0]=max(score[iscore[window1]])
            hits[0]=max(sc[iscore[window1]])
            z[1]=max(zarr[iscore[window2]])
            dlasc[1]=max(score[iscore[window2]])
            hits[1]=max(sc[iscore[window2]])
        endif else begin
            
            for j=0, ngap-1 do begin
                
                if j EQ 0 then begin 
                    window=where(gz LT gapz[0])
                    z[j]=max(zarr[iscore[window]])
                    dlasc[j]=max(score[iscore[window]])
                    hits[j]=max(sc[iscore[window]])
                endif
                
                if j GT 0 AND j LT ngap-1 then begin
                    window=where(gz GE gapz[j-1] $
                                 AND gz LT gapz[j])
                    z[j]=max(zarr[iscore[window]])
                    dlasc[j]=max(score[iscore[window]])
                    hits[j]=max(sc[iscore[window]])
                endif
                
                if j EQ ngap-1 then begin
                    window1=where(gz LT gapz[j] AND gz GE gapz[j-1])
                    z[ngap-1]=max(zarr[iscore[window1]])
                    dlasc[ngap-1]=max(score[iscore[window1]])
                    hits[ngap-1]=max(sc[iscore[window1]])
                    window2=where(gz GE gapz[j])
                    z[ngap]=max(zarr[iscore[window2]])
                    dlasc[ngap]=max(score[iscore[window2]])
                    hits[ngap]=max(sc[iscore[window2]])
                endif
                
            endfor
            
        endelse
    endelse
    no0 = where(z NE 0)
    z = z[no0]
    
    strct.nDLA2=ngap+1
    
    strct.DLA_zabs2=z
    strct.DLA_score2=dlasc
    strct.DLA_hits=hits
    
    print, 'z: '
    print, z
    
    
endif
;stop
end 


