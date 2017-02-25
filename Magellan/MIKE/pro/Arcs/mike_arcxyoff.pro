;+ 
; NAME:
; mike_arcxyoff   
;     Version 1.2
;
; PURPOSE:
;    This program is nearly obsolete.  Use mike_arcalign
;
;    Determine the pixel offset between a pair of Arc images
;    (typically the Template Arc and a different Arc) due to thermal
;    expansion of the instrument.
;
;    The program calls mike_arcxyoff_work which resamples a sub-region
;    of each Arc and then cross-correlates the images in 2D using a 2D
;    FFT.  Optionally, one can 'brute-force' the cross-correlation
;    analysis.  We have found both approaches give the same answer.
;   
;
; CALLING SEQUENCE:
;   
;  mike_arcxyoff, mike, side, 
;
; INPUTS:
;   mike     -  MIKE structure
;   side     -  Blue (1) or Red (2)
;   id1      -  Index of the first Arc image (e.g. Template Arc)
;   id2      -  Index of the second Arc image
;
; RETURNS:
;   xyoff -- The xy offsets
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /CHK  -- Shows a visual assessment of the cross-correlation
;  REGION -- Image sub-region for preforming the cross-correlation
;  MAGF  -- Level of resampling (default: 8L corresponding to 0.125
;           pixels)
;  /SLOW -- Turns on brute-force mode (NOT recommended)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_arcxyoff, mike, 1, 45L, 48L, xyoff
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_arcxyoff_work
;
; REVISION HISTORY:
;   13-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_arcxyoff_work, arc1, arc2, side, Magf=Magf, $
                  SLOW=slow, REGION=region, CHK=chk, XYINIT=xyinit

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'xyoff = mike_arcxyoff_work( arc1, arc2, side ) [v1.1]'
      return, -1
  endif 
  if not keyword_set( XYINIT ) then xyinit = [0., 0.]
  if not keyword_set( XOFF ) then xoff = 0
  if not keyword_set( yOFF ) then yoff = 0
  if not keyword_set( MAGF ) then Magf = 8L

  stop  ;; Odds are you do not want to be running this!!
  print, 'mike_arcxyoff: Correlating: ', arc1, ' ', arc2

  ;; Read images
  i1 = xmrdfits(arc1, /fscale, /silent)
  i2 = xmrdfits(arc2, /fscale, /silent)
  szw = size(i1, /dimensions)

  cbin = round(2048. / szw[0])
  rbin = round(4096. / szw[1])

  ;; Region
  if not keyword_set( REGION ) then begin
      if side EQ 1 then region = [ 575, 830, 1015L, 1270] $  ; Blue
      else region = [ 475, 730, 1000L, 1255]  ; Red
      ;; Binning
      region[0] =  round(region[0] * (2./cbin))
      region[1] =  round(region[1] * (2./cbin))
      region[2] =  round(region[2] * (2./rbin))
      region[3] =  round(region[3] * (2./rbin))
  endif

  print, 'mike_arcxyoff: Using region: ', region

  ;; Take sub-region
  i1 = i1[region[0]:region[1],region[2]:region[3]]
  i2 = i2[region[0]:region[1],region[2]:region[3]]

  sz = size(i1, /dimensions)


  ;; Cross correlate
  if keyword_set( SLOW ) then begin  ;; Brute force

     ;; REBIN
      r1 = rebin(i1, sz[0]*Magf, sz[1]*Magf)
      r2 = rebin(i2, sz[0]*Magf, sz[1]*Magf)
      ;; LOOP
      corr_val = fltarr(4*Magf, 2*Magf)
      
      x1_1 = 40 
      x2_1 = sz[0]*Magf - 1 - 40
      y1_1 = 20
      y2_1 = sz[1]*Magf - 1 - 20
      
      sub_r1 = r1[x1_1:x2_1,y1_1:y2_1]
      tot11 = total( sub_r1 * sub_r1 )
      
      for i=0L, 4*Magf-1 do begin  ;; Loop on 2 pix offset in x
          for j=0L, 2*Magf-1 do begin  ;; Loop on 1 pix offset in y
              
              x1_2 = 40 + xoff*Magf + (i - 2*Magf)
              x2_2 = sz[0]*Magf - 1 - 40 + xoff*Magf + (i - 2*Magf)
              y1_2 = 20 + yoff*Magf + (j - Magf)
              y2_2 = sz[1]*Magf - 1 - 20 + yoff*Magf + (j - Magf)
              
              ;; Auto Corr
              sub_r2 = r2[x1_2:x2_2,y1_2:y2_2]
              tot22 = total( sub_r2 * sub_r2 )
             
              ;; Cross-Corr
              corr_val[i,j] = total( sub_r1 * sub_r2 ) / sqrt( tot11 * tot22)
          endfor
      endfor
          
      if keyword_set( CHK ) then xatv, corr_val, /block
      mx = max( corr_val, imx)
      xopt = xoff + ((imx mod (4*Magf)) - 2*Magf) * (1./float(Magf))
      yopt = yoff + (imx / (4*Magf) - Magf) * (1./float(Magf))
      print, 'Offset: ', xopt, yopt
      stop ;; I am worried about a sign error

  endif else begin  ;; FFT

      sh1 = sz/2

      s1 = sqrt(i1>1)
      s2 = sqrt(i2>1)

      ;; Calculate FFT
      g1 = fft(s1-mean(s1))
      g2 = fft(s2-mean(s2))

      h1 = complexarr(sz[0]*MagF,sz[1]*MagF) 
      h1[0:sz[0]-1,0:sz[1]-1] = shift(g1,sh1[0],sh1[1])  
      h1 = shift(h1, -1 * sh1[0], -1 * sh1[1])

      h2 = complexarr(sz[0]*MagF,sz[1]*MagF) 
      h2[0:sz[0]-1,0:sz[1]-1] = shift(g2,sh1[0],sh1[1])  
      h2 = shift(h2, -1 * sh1[0], -1 * sh1[1])

      corr = fft( h1 * conj(h2), /inverse)


      ans = shift(float(corr), sh1[0]*MagF, sh1[1]*MagF)

      if keyword_set( CHK ) then $
        xatv, ans[64-2*Magf:64+2*Magf,64-Magf:64+Magf], /block

      ;; Find max
      mx = max( ans, imx)

      ym1 = ans[imx-sz[0]*MagF]
      yp1 = ans[imx+sz[0]*MagF]
      y_denom = (yp1 + ym1) - 2*mx
      ydel = y_denom LT 0 ? (imx / (sz[0] * MagF) - sh1[1]*MagF) + $ 
            (ym1-yp1)/(2*y_denom) : 0
      ydel = ydel / MagF
      yopt = yoff + ydel

      xm1 = ans[imx-1]
      xp1 = ans[imx+1]
      x_denom = (xp1 + xm1) - 2*mx
      xdel = x_denom LT 0 ? (imx mod (sz[0] * MagF) - sh1[0]*MagF) + $ 
          (xm1-xp1)/(2*x_denom) : 0
      xdel = xdel / MagF
      xopt = xoff + xdel

      ;; Find offset
      ;; xopt = xoff + ((imx mod (4*Magf + 1)) - 2*Magf) * (1./float(Magf))
      ;; yopt = yoff + (imx / (4*Magf+1) - Magf) * (1./float(Magf))

      ;; Sign flip (not exactly sure why)
      xyoff = [-xopt, -yopt] + xyinit

      print, 'mike_arcxyoff: Offsetting ', xyoff

  endelse

  return, xyoff
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_arcxyoff, mike, side, id1, id2, xyoff, Magf=Magf, $
                  SLOW=slow, REGION=region, CHK=chk

  if  N_params() LT 5  then begin 
      print,'Syntax - ' + $
        'mike_arcxyoff, mike, side, id1, id2, xyoff, MAGF=, REGION=, /CHK' + $
        '/SLOW, [v1.1]'
      return
  endif 
  ;; Image names
  i1nm = mike[id1].rootpth+mike[id1].img_root
  i2nm = mike[id2].rootpth+mike[id2].img_root

  ;; Here we go
  xyoff = mike_arcxyoff_work(i1nm, i2nm, side, XYINIT=mike[id1].arc_xyoff, $
                            MAGF=magf, SLOW=slow, CHK=chk, REGION=region)

  return
end
