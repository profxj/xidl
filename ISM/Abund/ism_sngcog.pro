;+ 
; NAME:
; ism_sngcog
;  V1.1
;
; PURPOSE:
;  Given a set of wavelengths, ew and sigew, calculate N,b for
;  a single line using standard COG analysis.
; CALLING SEQUENCE:
;   
;  ism_sngcog, cog, NLMT=Nlmt, BLMT=blmt, /CHICHK, PLTONLY=
;    NSTP=, BSTP=, PSFILE=, OUTFIL=, /EXACT, ZLBL=
;
; INPUTS:
;  wave -- Wavelengths (rest; Ang)
;  EW --  Rest equivalent widths (Ang)
;  sigEW -- Error in rest EW (Ang)
;  Nlmt] -- Range of column densities to explore (2 element array)
;  [blmt] -- Range of Doppler parameters to explore (2 element array)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /CHICHK -- Plot Chi^2 image
;  PLTONLY --  4-element array of N,b values and error for a plot
;  NSTP -- Number of steps to search N space [default: 100L]
;  BSTP -- Number of steps to search b space [default: 100L]
;  /EXACT -- Calculate EW exactly (Spline is generally good enough)
;  ZLBL= -- Label for Plot giving redshift of the absorber (string)
;  UPLIM -- file (like cog_fil) of lines to color upper limits on
;           curve
;  DEBLEND -- file with information of EW to remove for specific 
;             features because they are blended
;
; OPTIONAL OUTPUTS:
;   OUTFIL -- File with best fit values and error
;   PSFILE -- File for postscript plot
;
; COMMENTS:
;
; EXAMPLES:
;  fuse_cog, '/u/xavier/FUSE/data/PKS0405-12/Analysis/pks0405_abslin.fits', $
;    '/u/xavier/FUSE/data/PKS0405-12/Analysis/COG/Input/pks0405_z0918.cog', $
;    PSFIL='Figures/z0918_cog.ps', PLTONLY=[14.52, 38.2, 0.04, 1.8]
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   28-Mar-2006 Written by JXP  (based on fuse_cog)
;-
;------------------------------------------------------------------------------
function ism_sngcog_func, cog_strct, tau

  ;; Integrate (EXACT)
  ftau = spl_interp(cog_strct.tau, cog_strct.ftau, cog_strct.splint, $
                    tau)

  ; Other factors
  return, ftau
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ism_sngcog, cog, NLMT=Nlmt, BLMT=blmt, CHICHK=chichk, $
                NSTP=nstp, BSTP=bstp, PLOT=plot, STRCT=strct, $
                DCHISQ=dchisq, NARR=narr, BARR=barr, SILENT=silent, $
                CHISQ=chisq, NOSIG=nosig

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
      'ism_sngcog, wrest, ew, sigew, [Nlmt, blmt], /CHICHK, /PLOT,' + $
      'NSTP=, BSTP= (v1.1)' 
    return
  endif 

  if not keyword_set( NLMT ) then nlmt = [10., 18.]
  if not keyword_set( BLMT ) then blmt = [5., 50.]
  if not keyword_set( NSTP ) then nstp = 100L
  if not keyword_set( BSTP ) then bstp = 100L
  if not keyword_set(EXACT) then begin
      print, 'fuse_cog: Using spline interpolation for the COG.  '
      print, '      This assumes a Maxwellian profile'
      cogmax_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
      cog_strct = xmrdfits(cogmax_fil, 1, /silent)
  endif

  compile_opt strictarr
;  resolve_routine, 'x_calccog', /COMPILE_FULL_FILE, /EITHER
  
  ;; Reduce the EW 
  nlin = cog.nlin
  redew = cog.EW[0:nlin-1] / cog.wrest 
  redsigew = cog.sigEW[0:nlin-1] / cog.wrest 


  ;; CHISQ array
  chisq = dblarr(NSTP,BSTP)
  Nval = Nlmt[0] + (Nlmt[1]-Nlmt[0])*findgen(nstp)/float(nstp)
  bval = blmt[0] + (blmt[1]-blmt[0])*findgen(bstp)/float(bstp)

  Narr = Nval # replicate(1.,BSTP) 
  barr = replicate(1.,NSTP)  # bval

      
  for kk=0L, nlin-1 do begin
          
      ;; Get TAU
      tau = 1.497e-2*(cog.wrest[kk]*cog.f[kk]*1.e-8)*(10.^Narr)/(barr*1e5) 
      ftau = ism_sngcog_func(cog_strct, tau[*])
      ftau = reform(ftau, BSTP, NSTP)
      calcEW = 2.d*barr*ftau/(3.e5)

      ;; Calculate chisq
      chisq = chisq + ((calcEW-redew[kk])/(redsigew[kk]))^2 
  endfor

  if keyword_set( CHICHK ) then xatv, chisq/float((nlin-2)>1), /block
      
  ;; Find min chisq
  mnchi = min(chisq, ichi)
  ii = ichi MOD NSTP
  jj = (ichi/NSTP)
  Ngd = Nval[ii]
  bgd = bval[jj]

  dchisq = chisq - mnchi

  ;; Final structure
  strct = { $
          N: Ngd, $
          b: bgd, $
          Nsig: fltarr(2), $
          bsig: fltarr(2), $
          DCHISQ: dchisq, $
          barr: barr, $
          Narr: Narr $
          }
  

  ;; Plot
  if keyword_set(PLOT) then begin
      clr = getcolor(/load)
      contour, dchisq, Narr, barr, levels=[2.3, 6.17, 18.4], $
            c_annota=['1!9s!X', '2!9s!X', '3!9s!X'], color=clr.black, $
            background=clr.white, chars=1.5, c_charsi=1.5
  endif
  
  if keyword_set(NOSIG) then return

  ;; Calculate the probability array
  prob_arr = 1. - chisqr_pdf(chisq, 2)
  tot_prob = total(prob_arr)
  
  ;; Error in N
  sz = size(prob_arr)
  for i=0L, 2*NSTP/3 do begin
      if ii-i ge 0 and ii+i lt sz[1] then begin
          sub_prob = total( prob_arr[(ii-i)>0:ii+i,*] )
          if sub_prob GT 0.683*tot_prob then break
      endif else begin
          print,'ism_sngcog: cannot determine sigN'
          sigN = -9.99
          goto,assume_sigN
      endelse 
  endfor
  if i EQ 0 then stop
  sigN = (Nlmt[1]-Nlmt[0])*float(i)/float(nstp)
  assume_sigN:
  
  ;; Error in b
  for j=0L, 2*BSTP/3 do begin
      if jj-j ge 0 and jj+j lt sz[2] then begin
          sub_prob = total( prob_arr[*,(jj-j)>0:jj+j] )
          if sub_prob GT 0.683*tot_prob then break
      endif else begin
          print,'ism_sngcog: cannot determine sigb'
          sigb = -9.99
          goto,assume_sigb
      endelse
  endfor
  if j EQ 0 then stop
  sigb = (blmt[1]-blmt[0])*float(j)/float(bstp)
  assume_sigb:
  
  if not keyword_set(SILENT) then begin
      print, 'ism_sngcog: N = ', Ngd, '+/-', sigN
      print, 'ism_sngcog: b = ', bgd, '+/-', sigb
  endif

  return
end
