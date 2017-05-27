;+
; NAME:
;   deimos_arcfit
;
; PURPOSE:
;   Use preliminary solution from lin_arcfit_guess to trace arc lines
;
; CALLING SEQUENCE:
;   sigma = deimos_arcfit( slitarc, arcline_x, lamps, [wset, polyx=
;   ,s_coeff=, ncoeff=, slitwid= , polyflag=, /plot,  maxerr=maxerr)
;
; INPUTS:
;   slitarc    - extracted 2-D spectrum for a single slitlet
;   arcsat     - saturation mask (1=bad) for slitarc;
;                             set 2=bad for pixels >0% vignetted
;   arcline_x  - x position (from lin_arcfit_guess) of line
;   lamps      - lamp structure (see lamps.pro)
;
; OUTPUTS:
;   wset       - traceset (pix -> lambda)
;   sigma       - RMS uncertainty in wavelength solution
;   polyx      - optional output- polynomial expansion for
;                lambda(x,yc)
;   s_coeff    - polynomial expansion of slope of lines vs x
;   
;
; OPTIONAL KEYWORD:
;   lamdif     - list-fit residual, in Angstroms
;   ncoeff     - number of coefficients in fit along each line (default=2)
;   slitwid    - slit width [pix]
;   arcivar    - inverse variance of slitarc
;   polyflag   - flag to fit solution to tailored polynomial expansion
;   plot       - make atv plots
;   dlam       - mean lambda error in each row, averaged over lines
;   maxerr     - maximum error in line position allowed before a line is
;                thrown out
;   anamorph   - anamorphic factor (used in setting maximum row-to-row shifts)
;   ymid       - row to treat as the middle row
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   traceset2xy()
;   xy2traceset
;
; REVISION HISTORY:
;   Written by Finkbeiner, hacked by MD, jan02
;   2002-Jun-05  unhacked by DPF - Waimea
;   2002-Oct-17  nsmooth keyword in peaklines fixes alignment boxes - DPF
;   - this code is still in ugly form, and much longer than it should
;     be - DPF
;-
;------------------------------------------------------------------------------
function deimos_arcfit, slitarc, arcsat, arcline_x, lamps_in, wset, $
        lamdif=lamdif, ncoeff = ncoeff,  slitwid=slitwid,  arcivar= arcivar, $
               plot=plot,  polyflag=polyflag, polyx=polyx, s_coeff=s_coeff, $
               dlam= dlam, qaparams=qaparams, dline=dline, dlsig=dlsig, $
                        maxerr=maxerr, anamorph=anamorph, ymid=ymid


; -------- preliminaries
  npix = (size(slitarc, /dimens))[0]
  nrow = (size(slitarc, /dimens))[1]
  if n_elements(ymid) eq 0 then ymid = nrow/2

  if (NOT keyword_set(ncoeff)) then ncoeff = 4
  if (NOT keyword_set(slitwid)) then slitwid = 10. ; .75 asec 
  if n_elements(maxerr) eq 0 then maxerr=0.15
  if n_elements(anamorph) eq 0 then anamorph=1.6

  anascale=1.6/anamorph

  maxerr=maxerr < slitwid/8.

  errscale=maxerr/0.15

; separate out vignetted pixels from saturated

  nearvig = (arcsat AND 4b) eq 4b

  arcsat = arcsat AND 1b

; -------- only consider lines on the chip according to input guess
  onchip = (arcline_x LT npix) AND (arcline_x GT 0)

  fairlines=strpos(lamps_in.quality,'FAIR') GE 0
  
  w = where((lamps_in.good OR fairlines) AND onchip, nlamp)
  if nlamp gt 0 then begin
    lamps = lamps_in[w]
    arcline_xon = arcline_x[w]
    fairlines=fairlines[w]
  endif else begin
        message, ' no good or fair lines found!', /info
        print, ' no good or fair lines found!'
        return, 0
  endelse       

;  print, 'new arcfit'

; -------- only use nonvignetted / nonsaturated regions for onchip lines
  partvignette = total(nearvig, 2) GT 1. ;partially-vignetted regions

; expand by 30 pix
  partvignette = partvignette OR dilate(partvignette, intarr(61)+1)
; define fully unvignetted region
  nonvignette = partvignette EQ 0

; dilate is faster than smoothing floats...

; do not only look on center of slit for saturation - we could have a
;                                                     bad column there!

  censat = arcsat[*,ymid]

  unsat = (censat OR dilate(censat, intarr(31)+1)) EQ 0
 
  notsatlines = (unsat)[fix(arcline_xon)] gt 0
  notviglines = ( (nonvignette)[fix(arcline_xon)] gt 0) OR $
    (lamps.lambda lt 5500.)
  
  w = where(notsatlines AND notviglines, nlamp) 

  if nlamp gt 0 then begin
      lamps = lamps[w]
      arcline_xon = arcline_xon[w]
      fairlines=fairlines[w]
  endif


  vigstring='nonvignette / unsat mask: '+ $
    strcompress(strjoin(string( notsatlines AND notviglines,/PRINT)))

; -------- center up with peaklines and toss any faint lines
  arcline_xon_save = arcline_xon
  arcline_xon = peaklines(slitarc[*, ymid], arcline_xon, $ 
                width=12*(slitwid/5.33 < 1.), nsmooth=round(slitwid))

;keep only lines detected strongly
  goodline=where(arcline_xon gt 0, goodcnt,complement=badline) 

       if total(badline) gt -1 then $
		deletedcount=total(fairlines[badline] eq 0) $
	     else deletedcount=0

  if goodcnt gt 0 then begin
    arcline_xon = arcline_xon[goodline]
    lamps = lamps[goodline]
    nlamp = n_elements(lamps) 
    fairlines=fairlines[goodline]
  endif else nlamp = 0

  print, 'Lines selected: ', nlamp, '  # deleted: ', $
    deletedcount, format='(A,I4,A,I4)'
  if nlamp LE 0 then begin
        message, ' no good lines found!', /info
        print, ' no good lines found!'
        return, 0
     endif 



;---------------------------------------------------------------------------
; Trace arc lines on the 2D image
;---------------------------------------------------------------------------

; Allow for a shift of up to 2 pixels in the initial centers,
; but only 0.5 pixels while tracing

; xcen is [nrow,nline]

; here we add two to slitwid/2 because convergence is poor with
; nothing added, yet adding much more than 2 just gets you into
; neighboring lines.  - DPF
  rad = 2*anascale+(slitwid > 4.)/2.
  xcen = trace_crude(slitarc, arcivar, xstart=arcline_xon, ystart=ymid, $
                      radius=rad, yset=ycen, nave=1, nmed=1, $
            maxerr=0.5*anascale, maxshifte=.9*anascale, $
                     maxshift0=2.*anascale, xerr=xerr)

; Now fit traceset to xcen, this will "interpolate" over bad rows
; as well as extend trace off CCD if need be.

; TEST 2/25/02

;  xy2traceset, ycen, xcen, crudeset, yfit=xcrudefit,  ncoeff=2, /silent

  invvar=1.-(.999999d0*(xerr gt 10.))
  xy2traceset, ycen, xcen, crudeset, yfit=xcrudefit,  ncoeff=2, $
    invvar=invvar,/silent


; x1 is the output from tracecrude
; x2 is the recentered version of x1
; x3 is a linear fit to x2
; we use peaklines as a first rejection of weak lines, in table but
; not detected.


  x1 = xcrudefit


; median filter spectrum along central row


  slitarcsub = slitarc-median(slitarc)



; iterate this twice, for safety
  xx = trace_fweight(slitarcsub, xcrudefit, ycen, $
        radius=rad, xerr=xerrx,  invvar=arcivar) 

  x2 = trace_fweight(slitarcsub, xx, ycen, $
        radius=rad, xerr=xerr2,  invvar=arcivar) 

  x3 = x2*0.

  rr = findgen(nrow)

  npad = 6
  if nrow lt 21 then npad = (nrow-7)/2

  xerr2fit=xerr2	
  xerr2fit[0:npad-1, *] = 1e6
  xerr2fit[nrow-npad:nrow-1, *] = 1e6
  chifit=fltarr(nlamp)

  for i=0, nlamp-1 do begin

    fitline = svdfit(rr, x2[*, i], 2, measure_errors =xerr2fit[*, i], $
        chisq= chisq_line, yfit= lineout, /double)
    chifit[i]=chisq_line
    x3[*, i] =  float(lineout)
  endfor

; toss anything if you can't centroid to .2 pixel (weak lines)
  
; however, do not blindly toss lines that just have bad columns
; in the detector, etc. at their center
  badrow=xerr2fit gt 99
  weak = xerr2fit[ymid, *] GT maxerr/(1.+3.*fairlines) $
    OR chifit gt 20.*median(chifit)  ; This checks for CR-distorted lines
                                ;will average position along line
                                ;be more restrictive for 'fair' lines

  
  for i=0,nlamp-1 do if (weak[i] AND badrow[ymid,i]) then begin
      whok=where(badrow EQ 0 AND lindgen(nrow) ge npad AND $
                 lindgen(nrow) le nrow-npad-1,okct)
      if okct gt 5 then weak[i]=median(xerr2fit[whok,i]) GT maxerr
   endif


  wweak = where(weak, nweak)
  if nweak GT 0 then begin 

     whcare=where(fairlines[wweak] eq 0,carect)
     if carect gt 0 then begin
	     print, 'WARNING:  weak or CR-distorted lines found'
	     print, 'elements:   ',lamps[wweak[whcare]].element
	     print, 'wavelengths:',lamps[wweak[whcare]].lambda
;	     print, carect,' GOOD lines deleted'
     endif

     w = where(weak eq 0, nlamp)
     if nlamp lt 2 then begin 
        message, 'Oops - I just tossed all the lines!', /info
        print, 'Oops - I just tossed all the lines!'
        return, 0
     endif 
     lamps = lamps[w]
     fairlines=fairlines[w]
     x1 = x1[*, w]
     x2 = x2[*, w]
     x3 = x3[*, w]
     ycen = ycen[*, w]
     xerr2 = xerr2[*, w]
     xerr2fit = xerr2fit[*, w]
  endif 
;  xerror = reform(xerr2[nrow/2, *])
  xerror = reform(xerr2fit[ymid, *])

; try to avoid putting too much weight in any point
  xerror=sqrt(xerror^2+0.003^2*errscale^2)

;make pretty plots
  if keyword_set(plot) then begin 
     atv, slitarc, max=10000, min=0
     for i=0, nlamp-1 do atvplot, x1[*, i], ycen[*, i], color=1
     for i=0, nlamp-1 do atvplot, x3[*, i], ycen[*, i], color=3
  endif 

spline = 0  ; DPF - June 5 (DISABLED CODE REMOVED BY JAN 8/02)


;polyflag = 1  ;MD June 25


; BEGIN POLYFLAG METHOD

  if keyword_set(polyflag) then begin ;try special functional expansion
;     lambda(x,y)=P(x)[1+ (y-yc)*[s_0 + s_1*x  + s_2*x^2]
; setup unchanging parameters
    xx = findgen(npix)/(npix/2.) -1 ;need range -1,1 for orthogonal functions
    dx = 2./npix
    nrow2 = nrow/2
    ulim = 3*nrow/4
    llim = nrow/4
    dy = float(ulim-llim) ;number of rows

; fit lines, delete outliers in loop
    repeat begin
      degree = (6 < (nlamp - 2)) > 2
      t_order = 3 <  (nlamp -1) ;quadratic order if enough lines

      xx3 = reform(x3[nrow2, *])/(npix/2.) -1. ;position of line at central row

; double errors for FAIR lines in the fit

      if nlamp gt 1 then $
         polyx = svdfit(xx3, lamps.lambda, degree, /double, /legendre, $
            yfit=pxx, measure_errors = xerror*(1.+fairlines), $
                        chisq= chisq_polyx  ) $
      else polyx = lamps.lambda
;    print, 'fit to P(x): ', polyx
;    print, 'chisq: ', chisq_polyx 
      dldx = .5*(polyleg(xx3+dx, polyx) - polyleg(xx3-dx, polyx)) ;dlambda/dx
      slope = reform(x3[ulim, *] - x3[llim, *])/dy
      slopeerror=0.2*djsig(slope)*(1.+fairlines)*xerror/median(xerror)

      if nlamp gt 1 then  $
        s_coeff = svdfit(xx3, -slope*dldx/pxx, t_order, /double, $
            yfit= yout, measure_error =  slopeerror, $
	/legendre) $
      else s_coeff =  -slope*dldx/polyx

; check to see if slopes are correct
      if keyword_set(plot) then begin
        dxx = yout*nrow2*pxx/dldx
        for i=0, nlamp-1 do atvplot, [x3[nrow2,i]+dxx[i], $ 
           x3[nrow2, i]- dxx[i]], [0, nrow-1],  color=2
      endif
;fit slope of line vs x
;    print, 'slope coeff: ', s_coeff
;      wave = slitarc*0. 
      Px = polyleg(xx, polyx)
      fx = polyleg(xx, s_coeff)
      
      if nlamp gt 1 then dev = lamps.lambda-pxx $
	else dev = lamps.lambda-polyx
      bad = where(abs(dev) gt maxerr/(1.+3.*fairlines), nbad) ;flag bad lines
      if nbad gt 0 then begin
        worst = max(abs(dev), ind)

        if round(total(fairlines[ind])) NE 0 then $
	   fairstring = ' FAIR' else fairstring = ' GOOD'
        print, 'deleting bad line:  ', lamps[ind].element, $
          lamps[ind].lambda,fairstring

        keep = where(abs(dev) lt worst)
        nlamp = n_elements(keep)
        lamps = lamps[keep]
        x3 = x3[*, keep]
        x2 = x2[*, keep]
        xerr2 = xerr2[*, keep]
        xerr2fit = xerr2fit[*, keep]
        xerror = xerror[keep]
        fairlines=fairlines[keep]
      endif   
    endrep until nbad EQ 0 
;    print, 'finished lambda fits with ', nlamp, ' lines'

    if keyword_set(plot) then begin 
       for i=0, nlamp-1 do atvplot, x3[*, i], ycen[*, i], color=4
    endif

; do a test to check residuals of lines versus row number
;  perhaps this can show evidence of lambda shifts due to notches in
;  masks
    dlam = fltarr(nrow)
    dlamerr=dlam
    baddlam=dlam*0
    wave_est = fltarr(nlamp, nrow)
;    wave_raw=wave_est
    for i=0, nrow-1 do begin
;evaluate fit lambda at individual line centers
        whok=where(xerr2[i,*] lt 999,okct)
       if okct gt 0 then errlimit= (median(xerr2[i,whok]) > maxerr)*2. $
         else errlimit=maxerr*2.
       lli = polyleg(reform(2*x2[i, *]/npix-1.), polyx)*(1.+ (i-nrow2)* $
                     polyleg(reform(2*x2[i, *]/npix-1.), s_coeff))
; keep problematic pixels out of the solution
       goodpix=where(xerr2[i,*] lt errlimit AND (fairlines eq 0),goodct)
       if goodct gt 1 then errors=sqrt(xerr2[i,goodpix]^2+.003^2)

;       if goodct gt 1 then dlam[i] = median((lamps.lambda-lli)[goodpix]) $
       if goodct gt 1 then dlam[i] = $
         total(((lamps.lambda-lli)[goodpix]/errors^2)) / $
            total(1/errors^2) $
         else baddlam[i]=1
       if goodct gt 1 then dlamerr[i] = sqrt(1/total(1/errors^2)) $
         else baddlam[i]=1
       
;total(lamps.lambda - lli)/nlamp
; QA array - DPF
       wave_est[*, i] = lli + dlam[i]
;       wave_raw[*,i]=lli
    endfor

;tmplamps=lamps
;save,tmplamps,lli,x2,npix,polyx,nrow2,s_coeff,goodpix,xerr2,dlam,maxerr,errlimit,wave_est,wave_raw,baddlam,f='temp.sav'
    
    whok=where(baddlam eq 0,okct)
    if okct gt 0 then meddlam=median(dlam[whok]) else meddlam=0.

    whbad=where(baddlam ne 0,badct)
    if badct gt 0 then dlam[whbad]=meddlam

    dlamp = fltarr(nlamp)
    dlampsig = fltarr(nlamp)
    for i=0, nlamp-1 do dlamp[i]=median(lamps[i].lambda-wave_est[i, *])
    for i=0, nlamp-1 do dlampsig[i]=stdev(lamps[i].lambda-wave_est[i, *])

    tmp=wave_est
    for i=0,nlamp-1 do tmp[i,*]=tmp[i,*]-lamps[i].lambda
    for i=0,nrow-1 do tmp[*,i]=tmp[*,i]-dlam[i]

; end QA

;    for i=0, nrow-1 do $
;      wave[*, i] = Px*(1+ (i-nrow2)*fx)
    
  endif  else begin 

; END POLYFLAG METHOD

; alternative scheme:
;--------------------------------------------------------------------------
;   do final traceset fit
;--------------------------------------------------------------------------
; BEGIN TRACESET METHOD

     if (nlamp EQ 1) then ytmp = transpose(lamps.lambda * (dblarr(nrow)+1)) $
       else ytmp = lamps.lambda # (dblarr(nrow)+1)
     xerrmask = (xerr2 NE 999)  ; 0 means trace_fweight found no fit at all
     xivar = transpose((1./xerr2^2) * xerrmask)

     xy2traceset, transpose(double(x3)), ytmp, $
       wset, func=func, ncoeff=ncoeff, $
       maxiter=nlamp, maxrej=1, /sticky, $
       xmin=0, xmax=npix-1, invvar=xivar, yfit=yfit, /silent

     lamdif = (ytmp - yfit)*transpose(xerrmask)
     lam2 = total(lamdif, 2)/nrow
   
; check for bad line

     repeat begin 
        dev = max(abs(lam2), ind)
        lmask = bytarr(nlamp)+1B
        if dev gt 0.15 then begin 
           print, 'Bad line - off by ', dev, ' Angstrom', lamps[ind].lambda, $
             lamps[ind].element, format='(A,F8.4,A,F10.4,A3)'
           lmask[ind] = 0B
        endif 
        w = where(lmask, nlamp)
        if nlamp lt 5 then message, 'not enough lines - initial guess probably wrong', /info
     if nlamp lt 5 then print, 'not enough lines - initial guess probably wrong'

        lamps = lamps[w]
        xerr2 = xerr2[*, w]
        xerrmask = (xerr2 NE 999) 
        xivar = transpose((1./xerr2^2) * xerrmask)

        x2 = x2[*, w]
        x3 = x3[*, w]
      
       if (nlamp EQ 1) then ytmp = transpose(lamps.lambda * (dblarr(nrow)+1)) $
          else ytmp = lamps.lambda # (dblarr(nrow)+1)

ncoeff = 5
        xy2traceset, transpose(double(x2)), ytmp, $
          wset, func=func, ncoeff=ncoeff, $
          maxiter=nlamp, maxrej=8, /sticky, $
          xmin=0, xmax=npix-1, invvar=xivar, yfit=yfit, /silent

        zset = wset
        zset.coeff[0:ncoeff-1, *] = 0

        term4 = median(wset.coeff[4, *])
        term3 = median(wset.coeff[3, *])

        zset.coeff[4, *] = term4
        zset.coeff[3, *] = term3
        traceset2xy, zset, transpose(x2), dx3

; fix 3rd order term
        xy2traceset, transpose(double(x2)), ytmp-dx3, $
        wset, func=func, ncoeff=ncoeff-2, $
        maxiter=nlamp, maxrej=8, /sticky, $
        xmin=0, xmax=npix-1, invvar=xivar, yfit=yfit, /silent

        rowind = where(total(xivar lt .001,1) eq 0, ngood)
        if ngood lt 3 then stop

        coeff = dblarr(ncoeff, nrow)
        coeff[0, *] = wset.coeff[0, *]

        poly_iter, rowind, wset.coeff[1, rowind], 1, 3, coeff=coeff1
        coeff[1, *] = poly(findgen(nrow), coeff1)

if stdev(wset.coeff[2, rowind]) eq 0 then stop
        poly_iter, rowind, wset.coeff[2, rowind], 1, 3, coeff=coeff2
        coeff[2, *] = poly(findgen(nrow), coeff2)

        coeff[3, *] = term3
        coeff[4, *] = term4
;        coeff[5, *] = term5
        wset = zset ; get dimensions right again
        wset.coeff = coeff

        traceset2xy, wset, transpose(x2), wave_est
        lamdif = (ytmp - wave_est)*transpose(xerrmask)
;      lamdif = (10.d^ytmp)-10.d^yfit
        lam2 = total(lamdif, 2)/nrow

      endrep until total(abs(lam2) GT maxerr) EQ 0 

     if keyword_set(plot) then begin 
       for i=0, nlamp-1 do atvplot, x2[*, i], ycen[*, i], color=4
       xpred = traceset2pix(wset, ytmp[*, 0])
       for i=0, nlamp-1 do atvplot, xpred[i, *], ycen[*, i], color=5
;       for i=0, nlamp-1 do print, i, lamps[i].lambda, $
;         arcsat[x2[ymid, i], ymid], slitarc[x2[ymid, i], ymid] 
     endif

     pixnorm = 0 ; to reset array
     traceset2xy, wset, pixnorm, wave
     chipedge = minmax(wave)
     traceset2xy, wset, transpose(x2), wave_est
     dlamp = fltarr(nlamp)
     dlampsig = fltarr(nlamp)
     for i=0, nlamp-1 do dlamp[i]=median(ytmp[i, *]-wave_est[i, *])
     for i=0, nlamp-1 do dlampsig[i]=stdev(ytmp[i, *]-wave_est[i, *])



     if keyword_set(plot) then begin 
;        plot, transpose(ytmp-wave_est), ps=3,yr=[-.1,.1], $
;          ytit='Delta [Ang]'
;        for i=0, nlamp-1 do oplot, [0, nrow-1]+i*nrow, [1, 1]*dlamp[i]

        if mean(ytmp) lt 7800 then xrange = [6200, 8200] $
          else xrange = [7400, 9400]


;xrange = [5500, 8500]
        plot, ytmp[*, ymid], dlamp, ps=1, yr=[-.1,.1]*5, $
          title='Arc wavelength residual', xtit='Ang', ytit='Delta [Ang]', $
          xr=xrange, /xst
        errplot, ytmp[*, ymid], dlamp-dlampsig, dlamp+dlampsig
        oplot, [0, 1e4], [1, 1]*.01, line=1
        oplot, [0, 1e4], -[1, 1]*.01, line=1
        oplot, chipedge[0]*[1, 1], [-1, 1], line=2
        oplot, chipedge[1]*[1, 1], [-1, 1], line=2

     endif
  endelse

  dline = fltarr(n_elements(lamps_in))
  dlsig = fltarr(n_elements(lamps_in))
  for i=0, nlamp-1 do begin 
     dline[where(lamps_in.lambda eq lamps[i].lambda)]=dlamp[i]
     dlsig[where(lamps_in.lambda eq lamps[i].lambda)]=dlampsig[i]
  endfor 

  if nlamp gt 1 then scatter = stdev(dlamp) else scatter = 9999.9999

  if NOT keyword_set(coeff) then coeff = polyx


  degree=n_elements(coeff[*, 0])

  if degree eq 6 then print, 'Final coefficients:', coeff[*, 0], $
	format='(A21,2f9.3,4f9.5)' else $
	print, 'Final coefficients:', coeff[*, 0]

  djsigma=djsig(dlamp)*sqrt(((nlamp-1)>1.)/((nlamp-degree) > 1.))
  sigma=scatter*sqrt(((nlamp-1)>1.)/((nlamp-degree) > 1.))
 
; correct scatter for DOF!
  print,  nlamp,  ' lines used.  Stdev/djsig (dline+DOF): ',  $
    sigma, djsigma, $
    ' Angstrom', format='(I3,A,F9.5,F9.5,A)'

 print
  print,vigstring
  print



  return, sigma
end










