;+ 
; NAME:
; x_traceflat   
;    Version 1.1
;
; PURPOSE:
;    Traces the y-distortion of a flat of a number of slits
;
; CALLING SEQUENCE:
;  x_traceflat, img, tracestrct, VAR=, GAIN=, RN=, SAW=, $
;                ystrt=, /NOTRANS
;
; INPUTS:
;   img       - Input image
;
; RETURNS:
; trcstr -- A trace structure describing the curvature of individual
;           slits or orders
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  YSTRT    -- COLUMN/ROW to start at
;  /NOTRANS  -- Do not transpose original image
;  VAR=     -- Variance image input
;  GAIN,RN= -- Values to create a VAR array
;
; OPTIONAL OUTPUTS:
;  SAW=  -- The sawtooth image
;
; COMMENTS:
;
; EXAMPLES:
;   x_traceflat, img, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_traceflat, img, tracestrct, VAR=var, GAIN=gain, RN=rn, SAW=saw, $
                 ystrt=ystrt, NOTRANS=notrans



;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_traceflat, img, tracestrct, VAR=, GAIN=, RN=, SAW= '
    print,   '           ystrt= [v1.1]'
    return
  endif 

  sz = size(img)
  if sz[0] NE 2 then begin
      print, 'x_traceflat: need a 2-D image'
      return
  endif

;  Optional Keywords

  if not keyword_set( YSTRT ) then YSTRT = 400L
  if keyword_set( GAIN ) then tmpimg = img*gain else tmpimg=img
  if not keyword_set( RN ) then rn = 5.  ; Readnoise

;  Create inv variance as necessary

  if not keyword_set( VAR ) then var = (tmpimg + (tmpimg EQ 0) + RN^2)

;  Assume transpose

  if not keyword_set( NOTRANS ) then begin
      tmpimg = transpose(tmpimg)
      tmpvar = transpose(var)
  endif else tmpvar = var

;;;;;;;;;;;;;
;  Sawtooth and Variance

  ; Standard  -- This seems to be working best at the moment
  saw = shift(tmpimg,1) - shift(tmpimg,-1)   ; 1/2 of the edges
;  saw = abs(shift(tmpimg,1) - shift(tmpimg,-1))   ; 1/2 of the edges

  ; Gauss smoothing
;  g = gauss_kernel(1., nsub=5)
;  saw_kern = shift(g,-1) - shift(g,1)
;  saw = convol(tmpimg, saw_kern, /edge_truncate)

  sigsaw = sqrt( (shift(tmpvar,1) + shift(tmpvar,-1)) > 0 )
  ivar = 1./( shift(tmpvar,1) + shift(tmpvar,-1) )
  gdpt = where(abs(saw) GT 8.*sigsaw AND sigsaw NE 0., ngd)

  sz_saw = size(saw,/dimensions)
  mask = intarr(sz_saw[0],sz_saw[1])
  mask[gdpt] = 1

  ; Threshold -- Taking 1/2 of the median of all significant values
  thresh = 0.5 * median( abs(saw[gdpt]) )

;;;;;;;;;;;;;
;  Trace crude :: Positive Half
  if not keyword_set( SILENT ) then print, 'x_traceflat: Tracing'

  ; Take ystart=YSTRT for now.  Could vary this by slit position
  ;   Also, set all negative values to 0
;  xcen= trace_crude(saw, ivar, yset=ycen, thresh=thresh, $
;                     radius=1., ystart=YSTRT, xerr=xerr)
  xcen_pos = trace_crude((saw>0), ivar, yset=ycen_pos, thresh=thresh, $
                     radius=3., ystart=YSTRT, xerr=xerr_pos)

  ; Take ystart=YSTRT for now.  Could vary this by slit position
  ;   Also, set all negative values to 0
  xcen_neg = trace_crude(((-saw)>0), ivar, yset=ycen_neg, thresh=thresh, $
                     radius=3., ystart=YSTRT, xerr=xerr_neg)

  ; Put it together
  sz_p = size(xcen_pos,/dimensions)
  sz_n = size(xcen_neg,/dimensions)

  xcen = fltarr(sz_p[0],sz_n[1]+sz_p[1])
  xerr = fltarr(sz_p[0],sz_n[1]+sz_p[1])
  ycen = lindgen(sz_p[0]) # replicate(1,sz_n[1]+sz_p[1])

  txcen = fltarr(sz_p[0],sz_n[1]+sz_p[1]) ; temporary
  txerr = fltarr(sz_p[0],sz_n[1]+sz_p[1]) ; temporary

  ; Copy them in
  txcen[*,0:sz_p[1]-1] = temporary(xcen_pos)
  txcen[*,sz_p[1]:sz_p[1]+sz_n[1]-1] = temporary(xcen_neg)
  txerr[*,0:sz_p[1]-1] = temporary(xerr_pos)
  txerr[*,sz_p[1]:sz_p[1]+sz_n[1]-1] = temporary(xerr_neg)

  ; Sort

  srt = sort(txcen[YSTRT,*])
  for i=0L,sz_p[1]+sz_n[1]-1 do begin
      xcen[*,i] = txcen[*,srt[i]]
      xerr[*,i] = txerr[*,srt[i]]
  endfor
  delvarx, txcen, txerr
  

;;;;;;;;;;;;;
; Mask out bad points

  if not keyword_set( SILENT ) then print, 'x_traceflat: Masking'

  sz_cen = size(xcen,/dimensions)
  msk_cen = intarr(sz_cen[0], sz_cen[1])
  for i=0L,sz_cen[0]-1 do begin
      for j=0L,sz_cen[1]-1 do begin
          if((mask[long(xcen[i,j]),ycen[i,j]] EQ 1 OR $
              mask[long(xcen[i,j])+1,ycen[i,j]] EQ 1) AND $
             xerr[i,j] LT 0.25) then msk_cen[i,j] = 1
          
          ; Do not allow trace to skip > 50 pixels
          if i GT YSTRT then begin
              if max(msk_cen[i-50:i-1,j]) EQ 0 then msk_cen[i,j] = 0
          endif
      endfor
  endfor

;;;;;;;;;;;;;;;;
;  Find min, max of good regions

;
  miny = lonarr(sz_cen[1])
  maxy = lonarr(sz_cen[1])

  for i=0L,sz_cen[1]-1 do begin
      miny[i] = min( where(msk_cen[*,i] EQ 1) )
      maxy[i] = max( where(msk_cen[*,i] EQ 1) )
  endfor

  ; Create structure and transpose x,y by the variable name
  tracestrct = { $
                 xstrt: long(ystrt), $
                 leny: sz[2], $
                 xcen: temporary(ycen), $
                 ycen: temporary(xcen), $
                 yerr: temporary(xerr), $
                 msk_cen: temporary(msk_cen), $
                 minx: temporary(miny), $
                 maxx: temporary(maxy) }
                 

  delvarx, tmpimg, tmpvar, ivar

  ; Transpose saw
  if not keyword_set( NOTRANS ) AND arg_present(SAW) then saw = transpose(saw)

  return
end
