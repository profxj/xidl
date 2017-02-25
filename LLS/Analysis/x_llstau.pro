;+ 
; NAME:
; x_llstau
;  V1.2
;
; PURPOSE:
;    Given an NHI and b-value, calculate tau for the Lyman series + LL
;
; CALLING SEQUENCE:
;   
;   lls_struct, lls, [list]
;
; INPUTS:
;  [list] -- List of DLA base files.  [default:
;            '/u/xavier/LLS/Lists/tot_lls.lst']
;
; RETURNS:
;
; OUTPUTS:
;  lls  -- IDL DLA structure
;
; OPTIONAL KEYWORDS:
;  /ION - Input ionic column densities
;  /NOELM - Supress elemental and ion structures (saves memory)
;  /NORELM - Supress inputting Elemental [X/H] values
;  /NOHIS -- Suppress HI error in [X/H] values
;  /EW    -- Fill up the ion arrays with EW values instead of column
;            densities.  This reads the .EW files instead of .ion
;  /FINE  -- This fills up arrays for fine-structure states
;            (e.g. SiII*).  By default CII* is already considered 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_dlalst, dla, '/u/xavier/DLA/Lists/tot_dla.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   03-Jun-2008 Written by SB 
;-
;------------------------------------------------------------------------------
; these are 1.38 km/s pixels (1/50 SDSS pixels)
; wave = 10^(findgen(210000)*0.000002 + 2.7)

function x_llstau, wave, N, b, LIN_STRCT=lin_strct

  if N_PARAMS() LT 2 then begin
      print, 'Syntax  tau = x_llstau(wave, N, [b]) [v1.0]'
      return, 0
  endif

  if not keyword_set(b) then b = 30.  ; km/s
  bnorm = b/299792.4581  

  ;; Parse the HI lines
  if not keyword_set(LIN_STRCT) then begin
      readcol, getenv('XIDL_DIR')+'/Spec/Lines/Lists/atom.dat', $
               lbl, wave_rest, f, gamma, FORMAT='A,D,D,D', /sil
      HIlines = lindgen(31)
      wave_rest = wave_rest[HIlines]
      f = f[HIlines]
      gamma = gamma[HIlines]
      if arg_present(LIN_STRCT) then begin
          lin_strct = { $
                      wave_rest: wave_rest, $
                      f: f, $
                      gamma: gamma $
                      }
      endif
  endif else begin
      wave_rest = lin_strct.wave_rest
      f = lin_strct.f
      gamma = lin_strct.gamma
  endelse
  
  ;; Loop on Lyman series
  tau = wave*0.0
  for i=0,n_elements(wave_rest)-1 do begin
      
      vd = b/ (wave_rest[i] * 1.0e-13)
      vel_norm = abs(((wave/wave_rest[i])-1.0)/ bnorm)
      a = gamma[i]/ (4.0 * !Pi * vd)
      
      calc1 = where(vel_norm GE 19.0, complement=calc2) 
      
      vo = vel_norm*0.0	

      IF (calc1[0] NE -1) then begin 
          vel2 = vel_norm[calc1]*vel_norm[calc1]
          hh1 = 0.56419/vel2 + 0.846/(vel2*vel2)
          hh3 = -0.56 / (vel2 * vel2) 
          vo[calc1] = a * (hh1 + a * a * hh3) 
      ENDIF
      if (calc2[0] NE -1) then vo[calc2] = voigt(a,vel_norm[calc2]) 
      
      thistau = 0.014971475*(10.0^N)*f[i]*vo/vd
      tau = tau + thistau
  endfor

  zero = 912.7 + b/100.0        ; account for broadening in beginning of LL
  
  blueward = where(wave LT zero)
  if (blueward[0] EQ -1) then return, tau
  
  e = sqrt(zero/wave[blueward] - 1)
  corr = (zero / 912.7)^3       ; this should correct to within 0.1%
  
  newmodel = 10^(N-17.20) * (wave[blueward]/zero)^4 *  $
             exp(4.0*(1.0-atan(e)/e))/(1-exp(-6.283/e)) * corr
  
  tau[blueward] = newmodel
  
  return, tau
end
	
        
