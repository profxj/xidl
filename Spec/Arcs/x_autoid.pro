;+ 
; NAME:
; x_autoid   
;    Version 1.0
;
; PURPOSE:
;    Given a 1D arc-line spectrum from the WFCCD or ESI-LowD, this
;    code will automatically determine the wavelength solution.
;
; CALLING SEQUENCE:
;   x_autoid, img, [extrct], instr=, IDLIST=, GUESS=, LINELIST=,
;                       TOLER=
;
; INPUTS:
;   img       - Input arc image or spectrum
;   [extrct]    - Slit edges (default to image edges)
;
; RETURNS:
;
; OUTPUTS:
;   FITSTR -  Fit structure;  Set values prior to passing for 
;
; OPTIONAL KEYWORDS:
;   CPROG  -  Run the monte-carlo in the c program (advised)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_autoid, img, instr='WFCCDB'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_autoid, img, extrct, INSTR=instr, GUESS=guess, NSIG=nsig, $
              LINELIST=linelist, LINES=lines, DEBUG=debug, $
              FINFIT=finfit, CPROG=cprog, SILENT=silent, NGUESS=nguess, $
              FLG_AID=flg_aid, ZEROPT=zeropt, LBUFFER=lbuffer, $
              RBUFFER=rbuffer, STATUS=status

status=0

;  Error catching
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_autoid, img, [extrct],  INSTR=, GUESS=, TOLER=, IDLIST='
    print,   '           LINELIST=, /SILENT, ZEROPT= [v1.0]'
    status=1
    return
  endif 

if(n_elements(lbuffer) eq 0) then lbuffer=50L
if(n_elements(rbuffer) eq 0) then rbuffer=500L

;; Optional Keywords
  flg_aid = 1
;; pushed from 50 sigma to 30 sigma --- necessary to catch some very
;;                                      blue slits
  if not keyword_set( NSIG ) then nsig = 30.

;  Read in image if necessary
  dat = x_readimg(img, /fscale)

;  Extract 1D if necessary 
  
  sz = size(dat)
  if sz[0] EQ 2 then begin
      if not keyword_set( extrct ) then begin
          print, 'x_autoid: Must provide extraction structure'
          flg_aid = 0
          status=2
          return
      endif
      ;; EXTRACT (Assume already ov subtracted)
      spec = x_apall(dat, STRCT=extrct, /NOOV, /NOSKY)
  endif else spec = dat
  delvarx, dat


;  Setup Guess

  if keyword_set( INSTR ) then begin
      case instr of
          'WFCCDB': begin
              linelist = $
                getenv('IDLUTILS_DIR')+'/xidl/Spec/Arcs/Lists/wfccdB_HeNe.lst'
              nguess = [1500L, 40L, 30L]
              guess = dblarr(2,4)
              szguess = 4L
              guess[0,0] = 6000.
              guess[1,0] = 11000.
              guess[0,1] = -3.20
              guess[1,1] = -3.10
              guess[0,2] = -0.00015
              guess[1,2] =  0.00015
              guess[0,3] = 1.058D-7
          end
          'ESI_LWD': begin
              linelist = $
                getenv('IDLUTILS_DIR')+'/xidl/Spec/Arcs/Lists/XeHgNe.lst'
              nguess = [50L, 20L, 20L]
              guess = dblarr(2,5)
              szguess = 5L
              guess[0,0] = 3890
              guess[1,0] = 3920
              guess[0,1] = -0.800
              guess[1,1] = -0.815
              guess[0,2] = 0.00315
              guess[1,2] =  0.00325
              guess[0,3] = -2.4159573D-6
              guess[0,4] = 8.1675796D-10
          end
          else:
      endcase
  endif

  if not keyword_set( SZGUESS ) then szguess = 4L



  ;; Set number of pixels
  npix = n_elements(spec)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LINELIST

  if not keyword_set( lines ) then x_arclist, linelist, lines $
  else lines.flg_plt=0

  ;; Parse the 'good' lines
  gdlin = where(lines.flg_qual GE 5, ngd)
  if ngd EQ 0 then begin
      print, 'x_autoid: No lines for template!'
      flg_aid = 0
      status=3
      return
  endif

  gdwv = lines[gdlin].wave
  weight = float(lines[gdlin].flg_qual)
  mxgd = max(gdwv, min=mngd)

  ;; Final arrays
  gdpx = dblarr(ngd)
  gdflg = lonarr(ngd)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find all 50 sigma features in the Arc

  x_fndpeaks, spec, center, NSIG=nsig, /silent
  npk = n_elements(center)

  if keyword_set( DEBUG ) then $
    x_splot, spec, XTWO=center, YTWO=fltarr(npk), PSYM_Y2=1, /block

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Min Chi^2
  
  all_chi = fltarr(nguess[0],nguess[1],nguess[2])
  all_nsub = lonarr(nguess[0],nguess[1],nguess[2])

  if not keyword_set( SILENT ) then print, 'x_autoid: Monte carlo'
  
  if not keyword_set( CPROG ) then begin
      cent2 = center^2
      cent3 = center^3
      if szguess GT 4 then cent4 = center^4
      for q1=0L,nguess[0]-1 do begin
          Zro = guess[0,0] + q1*(guess[1,0]-guess[0,0])/float(nguess[0])
          for q2=0L,nguess[1]-1 do begin
              dlmb = guess[0,1] + q2*(guess[1,1]-guess[0,1])/float(nguess[1])
              for q3=0L,nguess[2]-1 do begin
                  nonl = guess[0,2] + q3*(guess[1,2]-guess[0,2])/float(nguess[2])
                  
                  ;; Calculate lambda for the peaks 
                  lambda = Zro + center*dlmb + cent2*nonl + cent3*guess[0,3]
                  if szguess GT 4 then lambda = lambda + cent4*guess[0,4]
                      
                  ;; Get the subset of model lines
; ??? hack: add buffer at edge of chip for vignetting probs? (mg&mb 022104)
;     ( this appears to work, but we have only tested the C version )
;t1 = Zro
;t2 = Zro + npix*dlmb + (npix^2)*nonl + (npix^3)*guess[0,3] 
                  buffer=50.
                  t1 = Zro + buffer*dlmb + (buffer^2)*nonl + $
                    (buffer^3)*guess[0,3] 
                  t2 = Zro + (npix-buffer)*dlmb + ((npix-buffer)^2)*nonl + $
                    ((npix-buffer)^3)*guess[0,3] 
                  if szguess GT 4 then t1 = t1 + (buffer^4)*guess[0,4]
                  if szguess GT 4 then t2 = t2 + ((npix-buffer)^4)*guess[0,4]
                  mxwv = t1 > t2
                  mnwv = t1 < t2
                  
                  subgd = where(gdwv GT mnwv AND gdwv LT mxwv, nsubgd)
                  subwv = gdwv[subgd]

                  ;; Find closest for each template line and calulate chisq
                  chisq = 0.d
                  for i=0L,nsubgd-1 do begin
                      min_sep = min(abs(lambda-subwv[i]))
                      chisq = chisq + min_sep^2/(dlmb/4.)^2
                  endfor
                  ;; Save chisq
                  if nsubgd LE 1 then begin
                      splog, 'BOMBING on x_autoid (0)'
                      status=4
                      return
                  endif
                  all_chi[q1,q2,q3] = chisq/float(nsubgd-1.)
                  ;; Save nsub
                  all_nsub[q1,q2,q3] = nsubgd
              endfor
          endfor
      endfor
  endif else begin
      ndim = 3L
      soname = filepath('libxmath.' + idlutils_so_ext(), $
                        root_dir=getenv('XIDL_DIR'), subdirectory='lib')
      retval = call_external(soname, 'arclincorr', $
                             ndim, all_chi, all_nsub, nguess, szguess, guess, $
                             long(npk), center, long(ngd), gdwv, weight, $
                             long(npix), long(lbuffer), long(rbuffer))
;                             long(npix), /UNLOAD)
  endelse

;;;;;;;;;;;;;;;;
  ;; MIN CHI
  min_chi = min(all_chi, imin)
          
  ;; Find the parameters
  q1 = imin MOD nguess[0]
  q2 = (imin MOD (nguess[0]*nguess[1]))/nguess[0]
  q3 = imin/(nguess[1]*nguess[0])
  Zro = guess[0,0] + q1*(guess[1,0]-guess[0,0])/float(nguess[0])
  dlmb = guess[0,1] + q2*(guess[1,1]-guess[0,1])/float(nguess[1])
  nonl = guess[0,2] + q3*(guess[1,2]-guess[0,2])/float(nguess[2])

;;;;;;; CHOOSE LINES FOR REFERENCE ;;;;;;;;;
  ;; Calculate Wavelengths
  if keyword_set( DEBUG ) then stop
  lambda = Zro + center*dlmb + (center^2)*nonl + (center^3)*guess[0,3] 
  if szguess GT 4 then lambda = lambda + (center^4)*guess[0,4]

  ;; Calculate tru_dlmb
  cent_plus = center+1.
  lmb_plus = Zro + cent_plus*dlmb + (cent_plus^2)*nonl + (cent_plus^3)*guess[0,3] 
  if szguess GT 4 then lmb_plus = lmb_plus + (cent_plus^4)*guess[0,4]
  tru_dlmb = abs(lambda-lmb_plus)
  
  ;; Keep when peak is within 5 pix 
  for i=0L,ngd-1 do begin
      min_sep = min(abs(lambda-gdwv[i]), imin)
      ;; Keep it? (5 pix separation)
      if min_sep LT abs(5*tru_dlmb[imin]) then begin
          gdpx[i] = center[imin]
          gdflg[i] = 1
      endif
  endfor

  fingd = where(gdflg EQ 1)
  if n_elements(fingd) LT 3 then begin
      print, 'not enough good lines'
      print, '(prochaska should document his own damn stops)'
      flg_aid = 0
      status=5
      return
  endif
;  x_splot, spec, XTWO=gdpx[fingd], YTWO=fltarr(n_elements(fingd)), PSYM_Y2=1, /block

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Fit 'reference' lines

  fitref = { fitstrct }
  fitref.func = 'POLY'
; order of the fit can't be greater than number of wavelengths MINUS ONE
  fitref.nord = (szguess-1) < (n_elements(fingd)-1L)
  fitref.lsig = 3.
  fitref.hsig = 3.
  fitref.niter = 3
  fitref.minpt = 5
  fitref.maxrej = 5
  
  pixfit = gdpx[fingd]
  wvfit = gdwv[fingd]
  ;; Add ZEROPT
  if keyword_set( ZEROPT) then begin
      pixfit = [pixfit, 0.]
      wvfit = [wvfit, zeropt]
  endif

  fit = x_fitrej(pixfit, wvfit, FITSTR=fitref)
  if fit[0] EQ -1 then begin
      print, 'x_autoid: Problem with fit to ref wave! x_autoid has failed!'
      flg_aid = 0
      status=7
      return
  endif
  
  wave = x_calcfit(findgen(npix), FITSTR=fitref)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Auto ID Remaining lines

  x_fndpeaks, spec, center, NSIG=3., PEAK=peak, /silent, /ALL

  xdat = findgen(npix)
  npk = n_elements(center)

  ; Add to lines
  for i=0L,npk-1 do begin
      ;; Require that a calibration line lies within +/- 5 pixels
      subwv = wave[(peak[i]-5)>0:(peak[i]+5)<(npix-1)]
      ;; Get min, max of range
      minwv = min(subwv, max=maxwv)
      gdln = where((lines.wave-minwv)*(lines.wave-maxwv) LE 0., nmatch)
      if nmatch GT 0 then begin
          ; Find the closest calibration line
          sep = abs(lines[gdln].wave - wave[peak[i]])
          minsep = min(sep, imin)
          if center[i] NE -1 then begin
              lines[gdln[imin]].pix = center[i]
              lines[gdln[imin]].flg_plt = 1 
          endif else lines[gdln[imin]].flg_plt = 0
      endif
  endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FIT THE GOOD LINES

  ; Just the good lines
  fitlin = where(lines.flg_plt EQ 1)

  ; Fit structure
  if not keyword_set( FINFIT ) then begin
      finfit = { fitstrct }
      finfit.func = 'POLY'
      finfit.nord = 4
      finfit.lsig = 3.
      finfit.hsig = 3.
      finfit.niter = 3
      finfit.minpt = 5
      finfit.maxrej = 5
  endif else begin ; Check it is a fit structure
      if tag_names(finfit, /structure_name) NE 'FITSTRCT' then begin 
          splog, 'finfit needs to be a fit structure!'
          status=4
          return
      endif
  endelse

  ;; Adjust max order of fit
  finfit.nord = finfit.nord < (n_elements(fitlin)-2)
  
  ;; Add ZEROPT
  pixfit = lines[fitlin].pix
  wvfit = lines[fitlin].wave
  if keyword_set( ZEROPT) then begin
      pixfit = [pixfit, 0.]
      wvfit = [wvfit, zeropt]
  endif
  ; FIT
  fit = x_fitrej(pixfit, wvfit, FITSTR=finfit)

  status=0
  return
end
