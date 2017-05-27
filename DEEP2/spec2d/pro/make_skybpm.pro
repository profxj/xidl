
;+
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; make_skybpm.pro
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; PURPOSE
;       The procedure make_skybpm.pro 
;
;
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; SYNTAX
;       make_skybpm
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; INPUTS
;       None
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; KEYWORDS
;       None
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; OUTPUTS
;       skybpm.fits.gz = a file called 
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; PROCEDURES CALLED 
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; EXAMPLES
;       
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; HISTORY
;       Created June 12, 2002 by jnewman and mcc.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-

PRO make_skybpm, maskfile

if n_elements(maskfile) eq 0 then maskfile = getenv('D2_RESULTS')+$
  '/1HSmask/1HSmask.4207.fits'

objdbase = mrdfits(maskfile, 2)

wh=where(objdbase.OBJNO MOD 1E6 GT 5E5 AND objdbase.OBJNO MOD 1E6 LT 6E5)

skyslits = objdbase[wh].slitn

;;;READ IN SPECTRA FROM .fits FILES
path = '/home/marc/d2_results/beta/4207/2002dec30/g600'
filestem = 'spSlit.4207.'
median_scale = 201.
smooth_scale = 5.
full_smooth = 5
thresh = 50.
thresh2 = thresh*5.

madearr = 0

skyslits = skyslits[[2, 4, 5,6]]

FOR k=0, n_elements(skyslits)-1 DO BEGIN
   
   IF skyslits[k] LT 10 THEN $
      stringnum = '00'+string(skyslits[k], format='(i1)') $
   ELSE IF skyslits[k] LT 100 THEN $
      stringnum = '0'+string(skyslits[k], format='(i2)') $
   ELSE stringnum = string(skyslits[k], format='(i3)') 

   fileR = filestem+stringnum+'R.fits.gz'
   fileB = filestem+stringnum+'B.fits.gz'
   tmp = findfile(path+fileR, count=ct) 
   isgood = ct NE 0

   IF isgood THEN BEGIN

      print, filer
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;READ-IN FILE fileR
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%
      s = MRDFITS(path + fileR, 2)
;;;DETERMINE ARRAY DIMENSION
      N = N_ELEMENTS(s.fullbkpt)
;;;EXTRACT THE SKY SPECTRUM FROM 5 COLUMNS
;      cols = findgen(5) / 4. * (N-1)
;      sky = s.skymodel[*,cols]
;      lambda = lambda_eval(s)
;      lambda = lambda[*, cols]
;      sortvec = SORT(lambda)
      lambda = s.fullbkpt
      sky = bspline_valu(lambda, s)
;;;SMOOTH SKY SPECTRUM TO DETERMINE SKY CONTINUUM
      sky = median(sky, 3)
      smooth_sky = djs_MEDIAN(sky, width=median_scale,boundary='reflect') 
;;;SUBTRACT SKY CONTINUUM FROM SKY SPECTRUM
      resids = sky - smooth_sky

;;;DEFINE THRESHOLD LEVEL
;thresh = 50.
;;;CREATE skybpm ARRAY
      M = N_ELEMENTS(sky)
      skybpmR1 = fltarr(M,2)
; skybpm2 = on bright sky lines (instead of off)
      skybpmR2 = skybpmR1
;;;FILL skybpm ARRAY ALONG DIMENSION skybpm[*,1]
      skybpmR1[(where((resids) GE thresh)), 1] = 0.
      skybpmR1[(where((resids) LT thresh)), 1] = 1.
      skybpmR2[(where((resids) GE thresh2)), 1] = 1.
      skybpmR2[(where((resids) LT thresh2)), 1] = 0.
;;;SMOOTH LINE DETECTIONS USING BOXCAR AND ROUND
;;;TO INTEGER VALUES (0,1)
      skybpmR1[*,1] = FLOOR(SMOOTH(skybpmR1[*,1], smooth_scale))
      skybpmR1[*,0] = lambda
;      skybpmR2[*,1] = FLOOR(SMOOTH(skybpmR2[*,1], 5))
      skybpmR2[*,0] = lambda

      IF madearr EQ 0 THEN lambdaR = lambda $
         ELSE lambdaR = [lambdaR, lambda]
      IF madearr EQ 0 THEN skyR = sky $
         ELSE skyR = [skyR, sky]
      IF madearr EQ 0 THEN smooth_skyR = smooth_sky $
         ELSE smooth_skyR = [smooth_skyR, smooth_sky]
      IF madearr EQ 0 THEN skybpmR = skybpmR1 $
         ELSE skybpmR = [skybpmR, skybpmR1]
      IF madearr EQ 0 THEN skybpmRs = skybpmR2 $
         ELSE skybpmRs = [skybpmRs, skybpmR2]

;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;;READ-IN FILE fileB
;;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      s = MRDFITS(path + fileB, 2)
;;;DETERMINE ARRAY DIMENSION
      N = N_ELEMENTS(s.fullbkpt)
;;;EXTRACT THE SKY SPECTRUM FROM 5 COLUMNS
;      cols = findgen(5) / 4. * (N-1)
;      sky = s.skymodel[*,cols]
;      lambda = lambda_eval(s)
;      lambda = lambda[*, cols]
;      sortvec = SORT(lambda)
      lambda = s.fullbkpt
      sky = bspline_valu(lambda, s)

;      s = MRDFITS(path + fileB, 1)
;;;DETERMINE ARRAY DIMENSION
;      N = N_ELEMENTS(s.lambda[0,*])
;;;EXTRACT THE SKY SPECTRUM FROM 5 COLUMNS
;      cols = findgen(5) / 4. * (N-1)
;      sky = s.skymodel[*,cols]
;      lambda = lambda_eval(s)
;      lambda = lambda[*,cols]
;      sortvec = SORT(lambda)
;      lambda = lambda[sortvec]
;      sky = sky[sortvec]
;;;SMOOTH SKY SPECTRUM TO DETERMINE SKY CONTINUUM
      sky = median(sky, 3)
      smooth_sky = djs_MEDIAN(sky, width=median_scale,boundary='reflect') 
;;;SUBTRACT SKY CONTINUUM FROM SKY SPECTRUM
      resids = sky - smooth_sky
;;;CREATE skybpm ARRAY
      M = N_ELEMENTS(sky)
      skybpmB1 = fltarr(M,2)
      skybpmB2 = skybpmB1
;;;FILL skybpm ARRAY ALONG DIMENSION skybpm[*,1]
      skybpmB1[(where((resids) GE thresh)), 1] = 0.
      skybpmB1[(where((resids) LT thresh)), 1] = 1.

      skybpmB2[(where((resids) GE thresh2)), 1] = 1.
      skybpmB2[(where((resids) LT thresh2)), 1] = 0.
;;;SMOOTH LINE DETECTIONS USING BOXCAR AND ROUND
;;;TO INTEGER VALUES (0,1)

      skybpmB1[*,1] = FLOOR(SMOOTH(skybpmB1[*,1], smooth_scale))
      skybpmB1[*,0] = lambda
;      skybpmB2[*,1] = FLOOR(SMOOTH(skybpmB2[*,1], 5))
      skybpmB2[*,0] = lambda



      IF madearr EQ 0 THEN lambdaB = lambda $
         ELSE lambdaB = [lambdaB, lambda]
      IF madearr EQ 0 THEN skyB = sky $
         ELSE skyB = [skyB, sky]
      IF madearr EQ 0 THEN smooth_skyB = smooth_sky $
         ELSE smooth_skyB = [smooth_skyB, smooth_sky]
      IF madearr EQ 0 THEN skybpmB = skybpmB1 $
         ELSE skybpmB = [skybpmB, skybpmB1]
      IF madearr EQ 0 THEN skybpmBs = skybpmB2 $
         ELSE skybpmBs = [skybpmBs, skybpmB2]

      IF madearr EQ 0 THEN madearr = 1
      
   ENDIF

ENDFOR

;;; FINALLY, COMBINE THE ARRAYS TO MAKE ONE
;;; MONSTER skybpm
skybpm = [skybpmB, skybpmR]
skybpms = [skybpmBs, skybpmRs]

index = SORT(skybpm[*,0])
skybpm[*,0] = skybpm[index,0]
skybpm[*,1] = skybpm[index,1]

index = SORT(skybpms[*,0])
skybpms[*,0] = skybpms[index,0]
skybpms[*,1] = skybpms[index,1]



range = minmax(skybpm[*, 0])
lambdaout = MIN(range) + 0.01 + 0.5 * $
            FINDGEN( (MAX(range)-MIN(range)) *2.) < (MAX(range)-0.01)
intersky = FLOOR( SMOOTH( INTERPOL(skybpm[*, 1], skybpm[*, 0], lambdaout), $
                          full_smooth) ) > 0
interskys = CEIL( SMOOTH( INTERPOL(skybpms[*, 1], skybpms[*, 0], lambdaout),$
                          3) ) > 0

;whred = where(lambdaout gt 9075., redct)
;if redct gt 0 then intersky(whred) =  0

skybpout = [[lambdaout], [intersky]]
skybpsout = [[lambdaout], [interskys]]

skybpm = skybpout

MWRFITS, skybpm, 'skyfreebpm.fits', /CREATE
MWRFITS, skybpsout, 'brightskybpm.fits', /CREATE
SPAWN, 'gzip -f skyfreebpm.fits'
SPAWN, 'gzip -f brightskybpm.fits'
sky = [skyB, skyR]
lambda = [lambdaB, lambdaR]
smooth_sky = [smooth_skyB, smooth_skyR]

index= SORT(lambda)
lambda = lambda[index]
sky = sky[index]

smooth_sky = smooth_sky[index]

resids = sky-smooth_sky


print, mean(skybpm[*, 1])
skybpmuse = FLOOR( INTERPOL(skybpm[*, 1], skybpm[*, 0], lambda) ) > 0
print, minmax(skybpmuse)

save, lambda, sky, f='tempsky.sav'
splot, lambda, resids, xr=[6300, 9100], yr=[-500,500], /xsty,/ysty, $
xtitle="wavelength (angstroms)", title="sky mask"
soplot, skybpm[*,0], -300.*(skybpm[*,1]-0.5), psym=4, color=1
soplot, lambda, resids*skybpmuse, color=200*255

END

