;takes wavelength, wavelength-to-pixel, absorption wavelength arrays,
;passes back score array of matched lines for each z in array of z's, zarr

pro sdss_fndna1, abswave, ZEM=zem, QSTRCT=strct


  strct={qalcharstrct}

;define score and z arrays



  ;;define rest wavelengths of lines to be searched for
  rwave= [5891.5833d, 5897.5581]
  nrwave = n_elements(rwave)
  maxabswave=max(abswave)
  if maxabswave LE 8000. then maxout=8000. else maxout=9200. 

  ;;  start/end wave
  strct.start_wave = (1+zem)*1230.d > 3820.
  strct.DLA_quality[0] = maxout < (1+zem+0.1)*5897.5581


 ;loop over z values from 1.45 to 5 in .001 increments, at each z
 ;calculate array of shifted wavelengths (abswavearr) corresponding to each rest wavelength
  nsearch = 6500L
  zarr = 0. + 0.0001*dindgen(nsearch)
  score= fltarr(nsearch)
  sc=fltarr(nsearch)

  for j=0L, n_elements(zarr)-1 do begin 
      if zarr[j] GT (zem+0.1) then continue
      abswavearr= rwave*(1.+zarr[j])
      goodlines= where((abswavearr GT (1+zem)*1230.d) AND $
                       (abs(abswavearr-5579.d) GT 5.d) AND $
                       (abs(abswavearr-6302.d) GT 5.d) AND $
                       abswavearr LT maxout, nposs) 
      if nposs NE 2 then continue
      
      sc[j]=0
    
      ;;at each z, search for a match of each shifted rest wavelength with a
      ;;wavelength in the absorption wavelength array, within a certain
      ;;tolerance, and count how many matches are found 

      for i=0, nrwave-1 do begin
          x= where(abs(abswavearr[i]- abswave) $
                   LT abswavearr[i]*(1./4000.), count) 
          if (count NE 0) then sc[j]= sc[j]+1
      endfor
   
      ;;define score array as scores for each .001 increment in z
      ;;define z arr as z values for each increment
      score[j]= sc[j]/float(nposs)
  endfor

  z=dblarr(100)
  dlasc=fltarr(100)
  hits=fltarr(100)

  iscore= where(score EQ 1.0, ngscore)
  
  if iscore[0] NE -1 then begin
      
      ghits=sc[iscore]
      gscore=score[iscore]
      gz= zarr[iscore]
      dgz= gz-shift(gz,1)
      dgz[0]= 0.
                                ;gwave=abswave[iscore]
      gapzindex= where(dgz GT .01, ngap)
      
      
      if ngap EQ 0 then begin
          z=median(gz)
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
      stop
      
  endif
end 


