;+ 
; NAME:
; x_traceslit   
;    Version 1.0
;
; PURPOSE:
;    This is a failed attempt at a tracing program.  Do not use!
;
; CALLING SEQUENCE:
;   
;   x_traceslit, img, pos, top, bottom
;   FFIT=, /INTER)
;
; INPUTS:
;   img       - Input image
;   pos       - [x,y] position on the CCD
;                      (default=1/3 max + min)
;
; RETURNS:
;
; OUTPUTS:
;   top - Trace of the top of the slit
;   bottom - Trace of the bottom of the slit
;
; OPTIONAL KEYWORDS:
;   GAIN - Gain of the CCD
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_traceslit, img, [252.3, 235.5], top, bottom
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   16-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_traceslit, img, pos, xtrc, top, bot, GAIN=gain


;  Error catching
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'x_traceslit, img, pos, xtrc, top, bottom [V1.0]'
    return
  endif 

  sz = size(img)
  if sz[0] NE 2 then begin
      print, 'x_traceslit: need a 2-D image'
      return
  endif

  if n_elements(pos) NE 2 then begin
      print, 'x_traceslit: Need position in [x,y]'
      return
  endif

;  Optional Keywords

  if not keyword_set( max_dslit ) then max_dslit = 10.  ; Max change in slitwid
  if not keyword_set( GAIN ) then gain = 1.
  if not keyword_set( NCOLM ) then ncolm = 11   ; Number of columns to collapse (odd)
  hfwdth = (ncolm-1)/2

  ipos = round(pos)

  ; Number of trace Regions (max) + Initialize the arrays

  nreg = long(sz[1] / float(ncolm)) + 1
  xpnt = dblarr(nreg)
  ytop = dblarr(nreg)
  ybot = dblarr(nreg)
  flgtrc = intarr(nreg)

; Start at the center

  ; Smash ncolm columns
  smash = djs_median(img[ipos[0]-hfwdth:ipos[0]+hfwdth, *], 1)
  ; Get rough value / 3

  ; Find crude edge
      for i=ipos[1],sz[2]-3 do begin  ; top
          sig = sqrt(smash[i]*gain/ncolm)
          if(smash[i+1] LE smash[i] - 10*sig AND $  ; 10 sigma criterion
             smash[i+2] LE smash[i+1] ) then break
      endfor
      if i GE sz[2]-3 then begin  ; Break if we hit the top of the CCD
          print, 'x_traceslit: Cant have this bad a problem!'
          stop
      endif
      ymax = i
      for i=ipos[1],2,-1 do begin ; bottom
          sig = sqrt(smash[i]*gain/ncolm)
          if(smash[i-1] LE smash[i] - 10*sig AND $ ; 10 sigma criterion
             smash[i-2] LE smash[i-1] ) then break
      endfor
      if i LE 0 then begin  ; Break if we hit the bottom of the CCD
          print, 'x_traceslit: Cant have this bad a problem!'
          stop
      endif
      ymin = i

  ; Rough slit width
  swidth = ymax-ymin+1

  ; Fit a Polynomial to the 'good' region
  ytmp = lindgen(ymax-ymin-1)+ymin+1
  fit = x_fitrej(ytmp, $
                 smash[ymin+1:ymax-1], $
                 'POLY', 2, LSIGMA=2.0, $
                 HSIGMA=2., NITER=3, MAXREJ=5, GDPT=gdpt, $
                 RMS=rms)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; TOP EDGE
  ; Find the points to fit in the top edge
  tmp = min(smash[ymax:(ymax+5)<sz[2]-1], i1)  ; Minimum near the edge
  i1 = i1 + ymax
  tmp = max( ytmp[gdpt], imx) ; First point within 2 sigma of the fit
  ymatch = 0.8* smash[ytmp[gdpt[imx]]] ; Value we are trying to match
  yfit = smash[ytmp[gdpt[imx]]-2:i1]  ; Sub 2 off imx to include 2 extra 
  tmp = lindgen(n_elements(yfit)) + ytmp[gdpt[imx]] - 2 ; Dummy array

  ; Fit the edge
  fit = x_fit(tmp, yfit, 'BSPLIN', n_elements(yfit)-2, FFIT=ffit)
  ; Find 0.8 edge
  y1 = x_fndfitval(tmp, yfit, 'BSPLIN', ffit, ymatch)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; BOTTOM EDGE
  ; Find the points to fit in the top edge
  tmp = min(smash[(ymin-5)>0:ymin], i1)  ; Minimum near the edge
  i1 = i1 + ( (ymin-5)>0 ) ; offset
  mxytmp = min( ytmp[gdpt], imn) ; First point within 2 sigma of the crude fit
  ymatch = 0.8* smash[ytmp[gdpt[imn]]] ; Value we are trying to match
  yfit = smash[i1:ytmp[gdpt[imn]]+2]  ; Add 2 to imn to include 2 extra pts
  tmp = lindgen(n_elements(yfit)) + i1 ; Dummy array

  ; Fit the edge
  fit = x_fit(tmp, yfit, 'BSPLIN', n_elements(yfit)-2, FFIT=ffit)
  ; Find 2/3 edge
  y2 = x_fndfitval(tmp, yfit, 'BSPLIN', ffit, ymatch)

  ; Save the values
  xpnt[0] = double(ipos[0])
  ytop[0] = y1
  ybot[0] = y2
  flgtrc[0] = 1


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; STEP TO LEFT

  qq = 1
  flg_trc = 1
  qqgd = 0
  while( flg_trc NE 0) do begin
      ; Get edges of the region
      xmin = ipos[0] - ncolm*qq - hfwdth
      if xmin LT 0 then break ; Dont bother with the left edge of the CCD
      xmax = ipos[0] - ncolm*qq + hfwdth

      ; Previous center 
      cent = round( (ytop[qqgd] + ybot[qqgd])/2.)
      ; Smash
      smash = djs_median(img[xmin:xmax,*], 1)
      ; Find crude edge
      for i=cent,sz[2]-3 do begin  ; top
          sig = sqrt(smash[i]*gain/ncolm)
          if(smash[i+1] LE smash[i] - 10*sig AND $  ; 10 sigma criterion
             smash[i+2] LE smash[i+1] ) then break
      endfor
      if i GE sz[2]-3 then begin  ; Break if we hit the top of the CCD
          flgtrc[qq] = 0
          qq=qq+1
          continue
      endif
      ymax = i
      for i=cent,2,-1 do begin ; bottom
          sig = sqrt(smash[i]*gain/ncolm)
          if(smash[i-1] LE smash[i] - 10*sig AND $ ; 10 sigma criterion
             smash[i-2] LE smash[i-1] ) then break
      endfor
      if i LE 0 then begin  ; Break if we hit the bottom of the CCD
          flgtrc[qq] = 0
          qq=qq+1
          continue
      endif
      ymin = i

      ; Quick check on slit width
      if (ymax-ymin+1) GT swidth+max_dslit then begin
          print, 'x_traceslit: Slit has gotten way too wide', xmin, xmax
          flgtrc[qq] = 0
          qq=qq+1
          continue
      endif

          
      ; Fit a Polynomial to the 'good' region
      ytmp = lindgen(ymax-ymin-1)+ymin+1
      svfit = x_fitrej(ytmp, $
                     smash[ymin+1:ymax-1], $
                     'POLY', 2, LSIGMA=2.0, $
                     HSIGMA=2., NITER=3, MAXREJ=5, GDPT=gdpt, $
                     RMS=rms)

   ;;;;
   ; TOP EDGE
      
       ; Find the points to fit in the top edge
      tmp = min(smash[ymax:(ymax+5)<sz[2]-1], i1) ; Minimum near the edge
      i1 = i1 + ymax
      tmp = max( ytmp[gdpt], imx) ; First point within 2 sigma of the fit
;      ymatch = svfit[gdpt[imx]]-10*sqrt(rms)  ; Find 10 sigma departure
      ymatch = 0.8* smash[ytmp[gdpt[imx]]] ; Value we are trying to match
      yfit = smash[ytmp[gdpt[imx]]-2:i1] ; Subtract 2 off imx to include 2 extra pts
      tmp = lindgen(n_elements(yfit)) + ytmp[gdpt[imx]] - 2 ; Dummy array

      ; Fit the edge
      fit = x_fit(tmp, yfit, 'BSPLIN', n_elements(yfit)-2, FFIT=ffit)
      ; Find 2/3 edge
      y1 = x_fndfitval(tmp, yfit, 'BSPLIN', ffit, ymatch)

      ; Bad fit?
      if y1 LT 0 then flg_top[qq] = 0
      
   ;;;;
   ; BOTTOM EDGE
      ; Find the points to fit in the top edge
      tmp = min(smash[(ymin-5)>0:ymin], i1) ; Minimum near the edge
      i1 = i1 + ( (ymin-5)>0 )  ; offset
      mxytmp = min( ytmp[gdpt], imn) ; First point within 2 sigma of the crude fit
      ymatch = 0.8* smash[ytmp[gdpt[imn]]] ; Value we are trying to match
;      ymatch = svfit[gdpt[imn]]-10*sqrt(rms)  ; Find 10 sigma departure
      yfit = smash[i1:ytmp[gdpt[imn]]+2] ; Add 2 to imn to include 2 extra pts
      tmp = lindgen(n_elements(yfit)) + i1 ; Dummy array
      
      ; Fit the edge
      fit = x_fit(tmp, yfit, 'BSPLIN', n_elements(yfit)-2, FFIT=ffit)
      ; Find 2/3 edge
      y2 = x_fndfitval(tmp, yfit, 'BSPLIN', ffit, ymatch)
      ; Bad fit?
      if y2 LT 0 then flg_bot[qq] = 0

      xpnt[qq] = double( (xmin+xmax)/2.)
      ytop[qq] = y1
      ybot[qq] = y2
      qqgd = qq
      qq = qq+1
  endwhile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; STEP TO RIGHT

  flg_trc = 1
  tt = 0
  nlft = qq ; Number of points to the left
  qqgd = 0

  while( flg_trc NE 0) do begin
      tt = qq-nlft  ; Index
      ; Get edges of the region
      xmax = ipos[0] + ncolm*tt + hfwdth
      if xmax GT sz[1] then break ; Dont bother with the right edge of the CCD
      xmin = ipos[0] + ncolm*tt - hfwdth

      ; Previous center 
      cent = round( (ytop[qqgd] + ybot[qqgd])/2.) 

      ; Smash
      smash = djs_median(img[xmin:xmax,*], 1)
      ; Find crude edge
      for i=cent,sz[2]-3 do begin  ; top
          sig = sqrt(smash[i]*gain/ncolm)
          if(smash[i+1] LE smash[i] - 10*sig AND $  ; 10 sigma criterion
             smash[i+2] LE smash[i+1] ) then break
      endfor
      if i GE sz[2]-3 then begin  ; Break if we hit the top of the CCD
          flgtrc[qq] = 0
          qq=qq+1
          continue
      endif
      ymax = i
      for i=cent,2,-1 do begin ; bottom
          sig = sqrt(smash[i]*gain/ncolm)
          if(smash[i-1] LE smash[i] - 10*sig AND $ ; 10 sigma criterion
             smash[i-2] LE smash[i-1] ) then break
      endfor
      if i LE 0 then begin  ; Break if we hit the bottom of the CCD
          flgtrc[qq] = 0
          qq=qq+1
          continue
      endif
      ymin = i

      ; Quick check on slit width
      if (ymax-ymin+1) GT swidth+max_dslit then begin
          print, 'x_traceslit: Slit has gotten way too wide', xmin, xmax
          flgtrc[qq] = 0
          qq=qq+1
          continue
      endif

          
      ; Fit a Polynomial to the 'good' region
      ytmp = lindgen(ymax-ymin-1)+ymin+1
      svfit = x_fitrej(ytmp, $
                     smash[ymin+1:ymax-1], $
                     'POLY', 2, LSIGMA=2.0, $
                     HSIGMA=2., NITER=3, MAXREJ=5, GDPT=gdpt, $
                     RMS=rms)

   ;;;;
   ; TOP EDGE
      
       ; Find the points to fit in the top edge
      tmp = min(smash[ymax:(ymax+5)<sz[2]-1], i1) ; Minimum near the edge
      i1 = i1 + ymax
      tmp = max( ytmp[gdpt], imx) ; First point within 2 sigma of the fit
;      ymatch = svfit[gdpt[imx]]-10*sqrt(rms)  ; Find 10 sigma departure
      ymatch = 0.8* smash[ytmp[gdpt[imx]]] ; Value we are trying to match
      yfit = smash[ytmp[gdpt[imx]]-2:i1] ; Subtract 2 off imx to include 2 extra pts
      tmp = lindgen(n_elements(yfit)) + ytmp[gdpt[imx]] - 2 ; Dummy array

      ; Fit the edge
      fit = x_fit(tmp, yfit, 'BSPLIN', n_elements(yfit)-2, FFIT=ffit)
      ; Find 2/3 edge
      y1 = x_fndfitval(tmp, yfit, 'BSPLIN', ffit, ymatch)
      
   ;;;;
   ; BOTTOM EDGE
      ; Find the points to fit in the top edge
      tmp = min(smash[(ymin-5)>0:ymin], i1) ; Minimum near the edge
      i1 = i1 + ( (ymin-5)>0 )  ; offset
      mxytmp = min( ytmp[gdpt], imn) ; First point within 2 sigma of the crude fit
      ymatch = 0.8* smash[ytmp[gdpt[imn]]] ; Value we are trying to match
;      ymatch = svfit[gdpt[imn]]-10*sqrt(rms)  ; Find 10 sigma departure
      yfit = smash[i1:ytmp[gdpt[imn]]+2] ; Add 2 to imn to include 2 extra pts
      tmp = lindgen(n_elements(yfit)) + i1 ; Dummy array
      
      ; Fit the edge
      fit = x_fit(tmp, yfit, 'BSPLIN', n_elements(yfit)-2, FFIT=ffit)
      ; Find 2/3 edge
      y2 = x_fndfitval(tmp, yfit, 'BSPLIN', ffit, ymatch)
      
      ; Bad fit?
      if (y1 LT 0) OR (y2 LT 0) then begin
          flgtrc[qq] = 0
          qq=qq+1
          continue
      endif

      xpnt[qq] = double( (xmin+xmax)/2.)
      ytop[qq] = y1
      ybot[qq] = y2
      flgtrc[qq] = 1
      qqgd = qq
      qq = qq+1
  endwhile
  
  gdpt = where(flgtrc NE 0)
  xtrc = temporary(xpnt[gdpt])
  top = temporary(ytop[gdpt])
  bot = temporary(ybot[gdpt])

end
