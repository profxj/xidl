;;;;;
;;;;
function brute_teff_lyman, ilambda, zem, EW_SPLINE=EW_SPLINE, bval=bval, $
                           LYMAN=lyman, POWERFN=powerfn, N_eval=n_eval, $
                           TEFF_LYMAN=teff_lyman, FNZ=fnz, DEBUG=debug, DOLYA=dolya

  if not keyword_set(bval) then bval = 35. ; km/s
  if not keyword_set(N_eval) then N_eval = 1000L
  if n_elements(ilambda) GT 1 then stop; return, 0
  lambda = ilambda[0]

  c = x_constants()

  if bval NE 35. then stop

  ;; Grab the Line list
  if not keyword_set(LYMAN_LINES) then begin
     restore, 'Lyman_lines.idl'
     lyman = lines
  endif

  ;; Will want to read the Voigt in from a FITS file
  ;; Voigt profile
  nvel = 20001L
  velo = -10000.d + findgen(nvel) ; km/s
  dvel = 1. ; km/s
  uval = velo / bval

  if keyword_set(DOLYA) then begin
     ;; Lya
     dwv = dvel * lines[0].wrest / (c.c/1e5) ;; Ang
     vd = bval/ (lines[0].wrest * 1.0d-13) ;; Frequency
     a = lines[0].gamma / (12.56637 * vd)
     vgt = voigt(a, uval)
     tau = 0.014971475*lines[0].f*vgt/vd ;; Normalized to N_HI = 1 cm^-2
  endif else begin
     ;; Lyb
     dwv = dvel * lines[1].wrest / (c.c/1e5) ;; Ang
     vd = bval/ (lines[1].wrest * 1.0d-13) ;; Frequency
     a = lines[1].gamma / (12.56637 * vd)
     vgt = voigt(a, uval)
     tau = 0.014971475*lines[1].f*vgt/vd ;; Normalized to N_HI = 1 cm^-2
  endelse
     

  if not keyword_set(POWERFN) then $
     init_madau95_fnz, powerfn
  FNZ=1
;     init_powerfn, powerfn

  ;; Find the lines
  gd_Lyman = where(lambda/(1+zem) LT lyman.wrest, nlyman)
  if nlyman EQ 0 then return, 0.

  ;; Begin looping
  
  lgNval = 11.5 + 10.5*dindgen(N_eval)/(N_eval-1) ;; This is base 10 [Max at 22]
  dlgN = lgNval[1]-lgNval[0]
  Nval = 10.d^lgNval
  teff_lyman = dblarr(nlyman)
  if keyword_set(DOLYA) then iqq = 0 else iqq = 1
  for qq=iqq,nlyman-1 do begin  ;; Would be great to do this in parallel... 
                               ;;; (Can pack together and should)
     ;; Redshift
     zeval = (lambda / lyman[qq].wrest) - 1
     if not keyword_set(FNZ) then $
        dxdz = abs(cosm_xz(zeval-0.1, /W05MAP,/silent)- $
                   cosm_xz(zeval+0.1,/W05MAP,/silent)) / 0.2 $
     else dxdz = 1. ;; Code is using f(N,z)

     ;; Evaluate f(N,X)
     log_fnz = eval_powerfn(POWERFN, lgNval, zeval, VEL_ARRAY=velo) * dxdz

     ;; Evaluate tau(z,N)
     tau_vN = Nval # tau  

     ;; Integrand
     intg = 10.d^(log_fnz) * (1. - exp(-1.*tau_vN))

     ;; Sum in N first
;     N_summed = total(intg*(Nval # replicate(1.,nvel)),  1) * dlgN * alog(10.)
;     teff_lyman[qq] = total(N_summed) * (dvel/(c.c/1e5)) * (1+zeval)
;     stop

     ;; Sum in z
     z_summed = total(intg, 2) * (dvel/(c.c/1e5)) * (1+zeval)
     teff_lyman[qq] = total(z_summed * Nval) * dlgN * alog(10.)
     stop

     ;; Debug
     if keyword_set(DEBUG) then begin
;        x_splot, lgNval, alog10(10.d^(log_fnX) * dxdz * dz * Nval), /bloc
;        x_splot, lgNval, total(10.d^(log_fnX) * dxdz * dz * Nval,/cumul) * dlgN * alog(10.) / teff_lyman[qq], /bloc
        printcol, lgnval, log_fnx, dz,  alog10(10.d^(log_fnX) * dxdz * dz * Nval)
        writecol, 'debug_file'+strtrim(qq,2)+'.dat', $
                  lgNval, restEW, log_fnX
        stop
     endif
     
  endfor
  
  return, total(teff_lyman)
end
