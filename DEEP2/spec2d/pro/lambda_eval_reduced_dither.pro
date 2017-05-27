;+
; NAME:
;   lambda_eval
;
; PURPOSE:
;   evaluates 2D lambda for a slitlet from restricted polynomial solution
;   plus the row by row correction correction term.
;
; CALLING SEQUENCE:
;   lambda = lambda_eval(calib)
; 
; INPUTS:
;   calib -- structure stored in calibSlit file
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;    double -- set to be sure to get a double-precision version of the
;              wavelength array
;
; OUTPUTS:
;   lambda -- 2-d lambda solution
;
; COMMENTS:
;    restricted polynomial is a very compact storage of full Lambda
;    solution in 2-d, but it is important to be able to quickly
;    evaluate full 2-d lambda solution on a given slitlet
;
; REVISION HISTORY:
;    md 1jul02, enroute to Kona
;    md 9jul02  dlam correction term added
;    jan 8 jul02 (hmmm...) made backwards compatible with old files
;----------------------------------------------------------------------
function lambda_eval_reduced_dither, slit,double=double

  if n_elements(double) eq 0 then double=0

  growfac=slit.dithersize
  shifts=slit.offsets

  nrow0 = slit.rawsize
  npix = (size(slit.ivar, /dimen))[0]
  nrow2 = nrow0/2
  nrow = (size(slit.ivar, /dimen))[1]
  roffset = mean(shifts)
  ldl = -10.
  udl =  10.
  if (double ne 0) then wave = dblarr(npix,nrow) else wave = fltarr(npix,nrow)

  tags = tag_names(slit)
  dlamexists = total(tags eq 'DLAM')
  if dlamexists then begin
    nexp=(size(shifts,/dimen))[0]
    dlam=fltarr(nrow)
    aveexp=fltarr(nrow)
    for i=0, nexp-1 do begin
      smin=shifts(i)
      smax=shifts(i)+nrow0-1
      dlam[smin:smax]+=slit.dlam[i]
      aveexp[smin:smax]+=1
    endfor
    dlam=dlam/aveexp
  endif else begin
    dlam = fltarr(nrow)
  endelse

;  plot,dlam[85:776]
;  oplot,dlam0


;  if total(tags eq 'LAMBDA0') then begin 
;     wave = slit.lambda0#(fltarr(nrow)+1.)+slit.lambda
;  endif 


  if total(tags eq 'TILTX') then begin 
     xx = findgen(npix)/(npix/2.) -1 ;need range -1,1 for orthogonal functions
     Px = polyleg(xx, slit.lambdax) ;4th order polynomial for lambda
     fx = polyleg(xx, slit.tiltx) ;2nd order polynomial for tilt of lines

     for i=0, nrow-1 do  wave[*, i] = Px*(1+ (i-nrow2 - roffset)*fx) + $
       ( (dlam[i] < udl) > ldl  )
;row by row correction from slit irregularities, determined from arc spectra
  endif 

  if total(tags eq 'COEFF') then begin 
     calib1 = {func: slit.func[0], xmin: slit.xmin[0], xmax: slit.xmax[0], $
            coeff: slit.coeff[*, *, 0]}
     traceset2xy, calib1, junk, wave
     for i=0, nrow-1 do  wave[*, i] = wave[*, i]+( (dlam[i] < udl) > ldl)
  endif 



 if double eq 0 then return, float(wave) else return,wave
end







