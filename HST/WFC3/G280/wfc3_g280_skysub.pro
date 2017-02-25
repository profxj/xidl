;+ 
; NAME:
;  wfc3_g280_skysub
;
; PURPOSE:
;   This code performs simple sky subtraction on the spectral image.
;
; CALLING SEQUENCE:
;   
;  specim_sky = wfc3_g280_skysub(wfc3_g280_strct, ii, specim, BEAM=beam, SKYDIR=, $
;                                savesky=, useuppersky=, $
;                                uselowersky=)
;
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;   ii -- the index of the object in the structure
;   specim -- 2D spectral image
;
; RETURNS:
;  specim_sky -- Sky subtracted 2D image 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  BEAM= -- which beam to sky subtract
;  SKY_DIR= -- Directory to store sky image
;  SAVESKY= -- saves the sky image
;  USEUPPERSKY= -- Uses the upper sky band for sky subtraction
;  USEUPPERSKY= -- Uses the upper sky band for sky subtraction
;
; COMMENTS:
;
; EXAMPLES:
;  specim_sky = wfc3_g280_skysub(wfc3_g280_strct, ii, specim, BEAM=beam, SKYDIR=skydir, $
;                                SAVESKY=savesky, USEUPPERSKY=useuppersky, $
;                                USELOWERSKY=uselowersky)
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   23-Dec-2010 Written by JXP/JMO
;   10-Jun-2016 Updated to deal with structure and multiple beams by MN
;------------------------------------------------------------------------------
function wfc3_g280_skysub, wfc3_g280_strct, ii, specim, BEAM=beam, SKYDIR=skydir, $
                           SAVESKY=savesky, USEUPPERSKY=useuppersky, $
                           USELOWERSKY=uselowersky

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
          'trace_strct =  wfc3_g280_skysub(specim, trace, SKY_IMG=) [v1.0]'
    return, -1
  endif 

  if not keyword_set(specim) then $
     specim=xmrdfits(wfc3_g280_strct(ii).spec_fil, $
                     7-3*wfc3_g280_strct(ii).chip)

  case beam of
     0: begin
        cnt=wfc3_g280_strct(ii).cnta
        trace_x=wfc3_g280_strct(ii).trace_xa
        trace_y_fit=wfc3_g280_strct(ii).trace_ya_fit
        trace_sigma_fit=wfc3_g280_strct(ii).trace_sigma_fita
        bm='A'
     end
     1: begin
        cnt=wfc3_g280_strct(ii).cntc
        trace_x=wfc3_g280_strct(ii).trace_xc
        trace_y_fit=wfc3_g280_strct(ii).trace_yc_fit
        trace_sigma_fit=wfc3_g280_strct(ii).trace_sigma_fitc
        bm='C'
     end
     else: stop
  endcase
  
  ;; create arrays
  sz = size(specim, /dimen)
  nax=cnt
  uppersky=fltarr(nax,20)
  lowersky=fltarr(nax,20)
  fullsky=fltarr(nax,40)
  upperskymed=fltarr(nax)
  lowerskymed=fltarr(nax)
  fullskymed=fltarr(nax)
  repsky=fltarr(nax,20)
  repskymed=fltarr(nax)
  
  ;; Create sky array
  for jj=0L,nax-1 do begin
     if round(trace_x[jj]) lt 0 or round(trace_x[jj]) ge sz(0) or $
        round(trace_y_fit[jj]+6*trace_sigma_fit[jj]+19) ge sz(1) or $
        round(trace_y_fit[jj]-6*trace_sigma_fit[jj]-19) lt 0 then continue
     uppersky[jj,*]=specim[round(trace_x[jj]),$
                           round(trace_y_fit[jj]+6*trace_sigma_fit[jj])+lindgen(20)]
     lowersky[jj,*]=specim[round(trace_x[jj]),$
                                                      round(trace_y_fit[jj]-6*trace_sigma_fit[jj])-lindgen(20)]
     if keyword_set(useuppersky) then lowersky = uppersky
     if keyword_set(uselowersky) then uppersky = lowersky
     fullsky[jj,0:19]=lowersky[jj,*]
     fullsky[jj,20:39]=uppersky[jj,*]
     repsky[jj,*]=[reform(specim[round(trace_x[jj]),round(trace_y_fit[jj]+6*trace_sigma_fit[jj])-1-lindgen(5)]), $
                   reform(specim[round(trace_x[jj]),round(trace_y_fit[jj]+6*trace_sigma_fit[jj])+20+lindgen(5)]), $
                   reform(specim[round(trace_x[jj]),round(trace_y_fit[jj]-6*trace_sigma_fit[jj])+1+lindgen(5)]), $
                   reform(specim[round(trace_x[jj]),round(trace_y_fit[jj]+6*trace_sigma_fit[jj])-20-lindgen(5)])]
  endfor
  
  ;; Reject
  djs_iterstat, fullsky, sigrej=2.5, mask=badpix, median=med_full

  ;; Grow 3
  grow = 1L
  szb = size(badpix, /dimens)
  cridx = where(badpix LE 0., ncr)
  if ncr GT 0 then begin
     xcr = cridx mod szb[0]
     ycr = cridx / szb[0]
     
     ximg = lindgen(szb[0]) # replicate(1L,szb[1])
     yimg = replicate(1L,szb[1]) # lindgen(szb[1])
     x0 = (xcr - grow) > 0
     x1 = (xcr + grow) < (szb[0]-1)
     y0 = (ycr - grow) > 0
     y1 = (ycr + grow) < (szb[1]-1)
  endif

  ;; Loop
  for qq=0L,ncr-1 do begin
      badpix[x0[qq]:x1[qq],y0[qq]:y1[qq]] = 0
  endfor

  ;; Simple box car median
  for jj=0,nax-1 do begin
     a = where(badpix[jj,*], na)
     ;; Enough for stats?
     if na LT 10 then fullskymed[jj]=med_full else $
        fullskymed[jj]=median(fullsky[jj,a])
     repskymed[jj]=median(repsky[jj,*])
  endfor

  ;; Filter
  svfilt = savgol(10, 10, 0, 2)
  sky_model = convol(fullskymed, svfilt, /edge_trunc)
  sky_img = sky_model # replicate(1., sz[1])
  
  ;; Subtract it
  sky_sub = specim
  sky_sub[(sz(0)-1)<round(trace_x[0])>0:(sz(0)-1)<round(trace_x[0]+nax-1)>0,*] = $
     specim[(sz(0)-1)<round(trace_x[0])>0:(sz(0)-1)<round(trace_x[0]+nax-1)>0,*] - sky_img

  ;; save the sky image
  if keyword_set(savesky) then begin
     skyfile=SKYDIR+wfc3_g280_strct(ii).name+'_sky'+bm+'.fits'
     mwrfits, sky_img, skyfile, /create
  endif
  
  case beam of
     0: wfc3_g280_strct(ii).skya=abs(repskymed-sky_model)
     1: wfc3_g280_strct(ii).skyc=abs(repskymed-sky_model)
     else: stop
  endcase
  
  ;; Return
  return, sky_sub
end
