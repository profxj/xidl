;+ 
; NAME:
; x_autoid   
;    Version 1.0
;
; PURPOSE:
;    Traces a series of arc lines
;
; CALLING SEQUENCE:
;   
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

pro x_newauto, img, extrct, INSTR=instr, GUESS=guess, NSIG=nsig, $
              LINELIST=linelist, LINES=lines, DEBUG=debug, $
              FINFIT=finfit, CPROG=cprog, SILENT=silent, NGUESS=nguess, $
              FLG_AID=flg_aid, ZEROPT=zeropt


;  Error catching
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_autoid, img, [extrct],  INSTR=, GUESS=, TOLER=, IDLIST='
    print,   '           LINELIST=, /SILENT, ZEROPT= [v1.0]'
    return
  endif 

;; Optional Keywords
  flg_aid = 1
  if not keyword_set( NSIG ) then nsig = 50.

;  Read in image if necessary
  dat = x_readimg(img, /fscale)

;  Extract 1D if necessary 
  
  sz = size(dat)
  if sz[0] EQ 2 then begin
      if not keyword_set( extrct ) then begin
          print, 'x_autoid: Must provide extraction structure'
          flg_aid = 0
          return
      endif
      ;; EXTRACT (Assume already ov subtracted)
      spec = x_apall(dat, STRCT=extrct, /NOOV, /NOSKY)
  endif else spec = dat
  delvarx, dat


  ;; Set number of pixels
  npix = n_elements(spec)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LINELIST

  linelist = '/home/xavier/idl/xidl/Spec/Arcs/Lists/kast_blue3.lst'
  if not keyword_set( lines ) then x_arclist, linelist, lines $
  else lines.flg_plt=0

  ;; Parse the 'good' lines
  gdlin = where(lines.flg_qual GE 5, ngd)
  if ngd EQ 0 then begin
      print, 'x_autoid: No lines for template!'
      flg_aid = 0
      return
  endif

  gdwv = lines[gdlin].wave
  weight = float(lines[gdlin].flg_qual)
  alllin = lines[where(lines.flg_qual NE 0)].wave
  nall = n_elements(alllin)
  mxgd = max(gdwv, min=mngd)
  ngd = n_elements(gdwv)
  if ngd NE 4 then stop

  ;; Final arrays
  gdpx = dblarr(ngd)
  gdflg = lonarr(ngd)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find all 50 sigma features in the Arc

  x_fndpeaks, spec, center, NSIG=nsig, /silent
  npk = n_elements(center)
  mn = dblarr(npk)
  wvsv = dblarr(npk)

  if keyword_set( DEBUG ) then $
    x_splot, spec, XTWO=center, YTWO=fltarr(npk), PSYM_Y2=1, /block

  ;; Find the 4 peaks
  gdlin = dblarr(4)
  for i=0L,3 do begin
      a = where(center GT float(i)*npix/4. AND $
                center LT float(i+1)*npix/4., na)
      if na EQ 0 then stop
      mx = max( center[a], imx)
      gdlin[i] = center[a[imx]]
  endfor
  print, gdlin
  ngd = 4L

  all_chi = fltarr(nall,nall,nall,nall)
  for i=0L,nall-ngd-1 do begin
      for j=i+1,nall-ngd do begin
          for k=j+1,nall-ngd+1 do begin
              for l=j+1,nall-ngd+2 do begin
                  xval = [alllin[i], alllin[j], alllin[k], alllin[l]]
;                  ffit = poly_fit(xval, gdwv, 2, chisq=chisq)
                  ffit = poly_fit(gdlin, xval, 2, chisq=chisq)
;                  if chisq LT 100. then begin
;                      ffit = poly_fit(xval, gdwv, 3)
;                      ffit = poly_fit(gdlin, xval, 3)
;                      wvcen = poly(center,ffit)
;                      scre = 0.
;                      for qq=0L,npk-1 do begin
;                          mn[qq] = min( abs(alllin-wvcen[qq]), imn)
;                          if mn[qq] LT 5. then begin
;                              scre = scre + 1.
;                              wvsv[qq] = alllin[imn]
;;                              msk[qq] = 1B
;                          endif
;                      endfor
;                      fit2 = poly_fit(center[where(msk EQ 1B)], $
;                                      wvsv[where(msk EQ 1B)], 5)
;                      tmp = poly([0.,1.], fit2)
;                      if scre Ge npk/2 then begin
;                          print, chisq, scre;, tmp[0], tmp[1]-tmp[0]
;;                          print, mn
;                      endif
;                  endif 
                  all_chi[i,j,k,l] = chisq
              endfor
          endfor
      endfor
  endfor
  stop

  fingd = where(gdflg EQ 1)
  if n_elements(fingd) LT 3 then begin
      stop
      flg_aid = 0
      return
  endif
;  x_splot, spec, XTWO=gdpx[fingd], YTWO=fltarr(n_elements(fingd)), PSYM_Y2=1, /block

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Fit 'reference' lines

  fitref = { fitstrct }
  fitref.func = 'POLY'
  fitref.nord = (szguess-1) < n_elements(fingd)
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
      stop
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
          print, 'x_autoid: finfit needs to be a fit structure!'
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

  return
end

;  wvsv = fltarr(npk)
;  msk = bytarr(npk)
;  mn = fltarr(npk)
;
;  all_chi = fltarr(npk,npk,npk,npk)
;  for i=0L,npk-ngd-1 do begin
;      for j=i+1,npk-ngd do begin
;          for k=j+1,npk-ngd+1 do begin
;              for l=j+1,npk-ngd+2 do begin
;                  xval = [center[i], center[j], center[k], center[l]]
;;                  ffit = poly_fit(xval, gdwv, 2, chisq=chisq)
;                  ffit = poly_fit(gdwv, xval, 2, chisq=chisq)
;                  if chisq LT 100. then begin
;;                      ffit = poly_fit(xval, gdwv, 3)
;                      ffit = poly_fit(gdwv, xval, 3)
;                      allpix = poly(alllin,ffit)
;                      scre = 0.
;                      msk[*] = 0B
;                      for qq=0L,npk-1 do begin
;                          mn[qq] = min( abs(allpix-center[qq]), imn)
;                          if mn[qq] LT 2. then begin
;                              scre = scre + 1.
;                              wvsv[qq] = alllin[imn]
;                              msk[qq] = 1B
;                          endif
;                      endfor
;                      fit2 = poly_fit(center[where(msk EQ 1B)], $
;                                      wvsv[where(msk EQ 1B)], 5)
;                      tmp = poly([0.,1.], fit2)
;                      if scre Ge npk/2 then begin
;                          print, chisq, scre, tmp[0], tmp[1]-tmp[0]
;;                          print, mn
;                      endif
;                  endif 
;                  all_chi[i,j,k,l] = chisq
;              endfor

;;          endfor
;      endfor
;  endfor
;  stop

