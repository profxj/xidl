;+ 
; NAME:
; x_ffttrace   
;    Version 1.1
;
; PURPOSE:
;    Traces the y-distortion of a flat using FFTs.  It is unlikely
;    this program works well.  I highly recommend tracing the
;    'sawtooth' image of a flat instead.
;
; CALLING SEQUENCE:
;   x_ffttrace, img
;
; INPUTS:
;   img    - Input image (flat)
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
;   x_ffttrace, img, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   29-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_ffttrace, img, xval, yval, dyval, PLT=plt, XSTRT=xstrt

;  Error catching
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_ffttrace, img, [V1.0]' 
    return
  endif 

  sz_img = size(img)
  if sz_img[0] NE 2 then begin
      print, 'x_ffttrace: need a 2-D image'
      return
  endif

;  Optional Keywords

  if not keyword_set( XSTRT ) then xstrt = 400
  if not keyword_set( XSTEP ) then xstep = 32
  if not keyword_set( YSTEP ) then ystep = 128
  if not keyword_set( MXDY ) then mxdy = 5.
  mxshft = long(mxdy)*2

;  Begin

  ; Create save arrays
  ;   Extra one is for zeropt
  ny = sz_img[2]/ystep
  dyval = fltarr( (sz_img[1]/xstep) + (sz_img[1] MOD xstep NE 0)  + 1, ny)
  xval = fltarr( n_elements(dyval), ny) - 100.
  yval = fltarr( n_elements(dyval), ny)

  
  ref = fltarr(ystep)
  smsh = fltarr(ystep)
;  for j=ystep, sz_img[2]-1, ystep do begin
  jcnt = 0
  for j=ystep, 14*ystep, ystep do begin
      ; Reference 
      x1 = (xstrt/xstep) * xstep
      x2 = x1 + xstep
      ; Smash + Truncate at the end
      tmp = djs_median(img[x1:x2,j:(j+ystep-1)<(sz_img[2]-1)],1) ; Smash
      ntmp = n_elements(tmp)
      ref[0:ntmp-1] = tmp
      if ntmp NE ystep then ref[ntmp:ystep-1] = median(tmp) ; Pad with median
      
      ; FFT
      fft_ref = fft(ref)

      ; yval
      xval[0,jcnt] = (x1+x2)/2.
      dyval[0,jcnt] = 0.
      yval[0,jcnt] = (2*j+ystep)/2.

      cnt = 1L
      clr = getcolor(/load)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Head left
      for i=x1, 1, -xstep do begin
          xmn = (i-xstep) > 0
          xmx = i
          ; Smash
          tmp = djs_median(img[xmn:xmx,j:(j+ystep-1)<(sz_img[2]-1)],1) 
          ntmp = n_elements(tmp)
          smsh[0:ntmp-1] = tmp
          ; Pad with the median
          if ntmp NE ystep then smsh[ntmp:ystep-1] = median(tmp) 

          ; FFT
          fft_smsh = fft(smsh)

          ; Correlate
          corr = fft( fft_ref * conj(fft_smsh), /inverse)

          ; Shift to make things easy
          if cnt GT 1 then ishft = long(mxshft + dyval[cnt-1]) else ishft=mxshft
          if ishft MOD 2 EQ 1 then ishft=ishft+1 ; Make sure it is even
          
          ; Shift

          shft = shift(double(corr),ishft)
          dyval[cnt,jcnt] = x_maxspln( findgen(ishft+1), shft[ishft/2:3*ishft/2] ) $
            - float(ishft/2)
          xval[cnt,jcnt] = (xmx+xmn)/2.
          yval[cnt,jcnt] = yval[0,jcnt] - dyval[cnt,jcnt]
          
          ; Plot
          if keyword_set( PLT ) then begin
              plot, ref
              oplot, smsh, color=clr.red
              stop
          endif
          
          cnt = cnt + 1
      endfor

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Head Right
      for i=x2, sz_img[1]-1, xstep do begin
          xmn = i
          xmx = (i+xstep) < (sz_img[1]-1)
          ; Smash
          tmp = djs_median(img[xmn:xmx,j:(j+ystep-1)<(sz_img[2]-1)],1) 
          ntmp = n_elements(tmp)
          smsh[0:ntmp-1] = tmp
          ; Pad with the median
          if ntmp NE ystep then smsh[ntmp:ystep-1] = median(tmp) 

          ; FFT
          fft_smsh = fft(smsh)

          ; Correlate
          corr = fft( fft_ref * conj(fft_smsh), /inverse)

          ; Shift to make things easy
          if i NE x2 then ishft = long(mxshft - dyval[cnt-1]) else ishft=mxshft
          if ishft MOD 2 EQ 1 then ishft=ishft+1 ; Make sure it is even
          
          ; Shift

          shft = shift(double(corr),ishft)
          dyval[cnt,jcnt] = x_maxspln( findgen(ishft+1), shft[ishft/2:3*ishft/2] ) $
            - float(ishft/2)
          xval[cnt,jcnt] = (xmx+xmn)/2.
          yval[cnt,jcnt] = yval[0,jcnt] - dyval[cnt,jcnt]

          ; Plot
          if keyword_set( PLT ) then begin
              plot, ref
              oplot, smsh, color=clr.red
              stop
          endif

          cnt = cnt + 1
      endfor
      jcnt = jcnt + 1
  endfor

end
          
