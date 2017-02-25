;+ 
; NAME:
; sdss_metals
;    Version 1.1
;
; PURPOSE:
;   Identify metal-line systems given a set of absorption lines
;
; CALLING SEQUENCE:
;  sdss_metals, abswave, zem, WAVE=, FLUX=, QSTRCT=
;
; INPUTS:
;  abswave -- Absorption line array
;  zem     -- Emission redshift of the QSO
;
; RETURNS:
;
; OUTPUTS:
;  QSTRCT= -- Structure containing the key output of this routine
;
; OPTIONAL KEYWORDS:
;  WAVE= -- Wavelength array for plotting (not required)
;  FLUX= -- Flux array for plutting (not required)
;  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Apr-2004 Written by SHF
;-
;------------------------------------------------------------------------------
;takes wavelength, wavelength-to-pixel, absorption wavelength arrays,
;passes back score array of matched lines for each z in array of z's, zarr
pro sdss_metals, abswave, zem, WAVE=wave, FLUX=flux, QSTRCT=strct, ZMIN=zmin

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'sdss_metals, abswave, zem, WAVE=, FLUX=, QSTRCT= [v1.1]'
    return
  endif 

  if not keyword_set( ZMIN ) then zmin = 1.6

  strct={qalcharstrct}
  
;define score and z arrays
  
  nsearch = 50000L
  
  score= fltarr(nsearch)
  sc=fltarr(nsearch)

;define rest wavelengths of lines to be searched for

rwave= [1260.4221d, 1302.1685d, 1304.3702d, 1334.5323d, 1526.7066d, $
        1608.4449d, 1670.7874d, 1808.0126d, $
        2344.2140d, 2382.7650d, 2586.6500d, 2600.1729d, 2796.3520d, 2803.5310d]

nrwave = n_elements(rwave)

maxabswave=max(abswave)
if maxabswave LE 8000. then maxout=8000. else maxout=9200. 

;loop over z values from 0 to 5 in .0001 increments, at each z

;calculate array of shifted wavelengths (abswavearr) corresponding to each rest wavelength

  azarr= dindgen(nsearch) * 0.0001d
  gd = where(azarr GE ZMIN and azarr LT ZEM, ngd)
  for k=0L, ngd-1 do begin 
    
      ;; Indices
      j = gd[k]
      zarr = azarr[j]
    
      abswavearr= rwave*(1.+zarr)
      goodlines= where((abswavearr GT (1+zem)*1230.d) AND $
                       (abs(abswavearr-5579.d) GT 5.d) AND $
                       (abs(abswavearr-6302.d) GT 5.d) AND $
                       abswavearr LT maxout, nposs) 
      
; Subtract 1 from nposs if 1808 is valid
      l1808 = 1808.0126d * (1+zarr)
      if (l1808 GT (1+zem)*1230.d) AND $
        (abs(l1808-5579.d) GT 5.d) AND $
        (abs(l1808-6302.d) GT 5.d) AND $
        (l1808 LT maxout) then nposs=nposs-1
      
      
;at each z, set score, sc, initially to 0
      
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
      
      
;define score array as scores for each .0001 increment in z
;define z arr as z values for each increment
      
      score[j]= sc[j]/float(nposs)
      
      if sc[j] LT 2 then score[j] = 0.
      
  endfor
  
  z=dblarr(100)
  dlasc=fltarr(100)
  hits=fltarr(100)
  
  iscore= where(score GE 0.6, ngscore)
  
  if iscore[0] EQ -1 then begin
      print, 'no hits!'
      
      strct.nDLA2=0    
      strct.DLA_zabs2=0.
      strct.DLA_score2=0.
      strct.DLA_hits=0.
      
      
      
  endif else begin
      ghits=sc[iscore]
      gscore=score[iscore]
      gz= azarr[iscore]
      dgz= gz-shift(gz,1)
      dgz[0]= 0.
                                ;gwave=abswave[iscore]
      gapzindex= where(dgz GT .01, ngap)
      
      
      if ngap EQ 0 then begin
          z=max(gz)
          dlasc=max(gscore)
          hits=max(ghits)
          
    endif else begin
        
        gapz= azarr[iscore[gapzindex]]
        
             
        if ngap EQ 1 then begin
            window1=where(gz LT gapz[0])
            window2=where(gz GE gapz[0])
            z[0]=max(azarr[iscore[window1]])
            dlasc[0]=max(score[iscore[window1]])
            hits[0]=max(sc[iscore[window1]])
            z[1]=max(azarr[iscore[window2]])
            dlasc[1]=max(score[iscore[window2]])
            hits[1]=max(sc[iscore[window2]])
        endif else begin
            
            for j=0, ngap-1 do begin
                
                if j EQ 0 then begin 
                    window=where(gz LT gapz[0])
                    z[j]=max(azarr[iscore[window]])
                    dlasc[j]=max(score[iscore[window]])
                    hits[j]=max(sc[iscore[window]])
                endif
                
                if j GT 0 AND j LT ngap-1 then begin
                    window=where(gz GE gapz[j-1] $
                                 AND gz LT gapz[j])
                    z[j]=max(azarr[iscore[window]])
                    dlasc[j]=max(score[iscore[window]])
                    hits[j]=max(sc[iscore[window]])
                endif
                
                if j EQ ngap-1 then begin
                    window1=where(gz LT gapz[j] AND gz GE gapz[j-1])
                    z[ngap-1]=max(azarr[iscore[window1]])
                    dlasc[ngap-1]=max(score[iscore[window1]])
                    hits[ngap-1]=max(sc[iscore[window1]])
                    window2=where(gz GE gapz[j])
                    z[ngap]=max(azarr[iscore[window2]])
                    dlasc[ngap]=max(score[iscore[window2]])
                    hits[ngap]=max(sc[iscore[window2]])
                endif
                
            endfor
            
        endelse
    endelse
    
;stop
    strct.nDLA2=ngap+1
    
;print, 'z:', 'wave:', 'hits:', 'score'
                                ;printcol, hitz, hitwave, hitnumber, hitscore
    strct.DLA_zabs2=z
    strct.DLA_score2=dlasc
    strct.DLA_hits=hits
    
    ;print, 'z, hits, score:'
    ;printcol, z, hits, dlasc
;x_specplot, flux, fltarr(3900), inflg=4, wave=wave, /qal, /block 
;stop
    
    
endelse
;stop
end 


