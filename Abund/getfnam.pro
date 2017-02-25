;+ 
; NAME:
; getfnam
;  Version 1.1
;
; PURPOSE:
;    Given an atomic wavelength, return the fvalue and name of the
;    transition.  Wavelength needs to be within 0.03Ang of the 
;    desired transition.  By default, the file $XIDL_DIR/Spec/Lines/
;    Lists/qal.lst is read.  These are generally the same values as in
;    Morton (2003) or newer.
;
; CALLING SEQUENCE:
;   
;   getfnam, wave, fval, [nam], FIL=fil
;
; INPUTS:
;   wave       - ionic transition (Ang)
;
; RETURNS:
;
; OUTPUTS:
;   fval       - oscillator strength
;
; OPTIONAL KEYWORDS:
;  REF=   Reference flag (interpreted by $XIDL_DIR/Spec/Lines/Lists/ref_lin.dat)
;  TOLER= Value that a line must be within in the linelist to be
;  identified as the line of interest [Default: 0.003Ang]
;  /CLOSE =  Set TOLER = 0.2Ang if TOLER is not also set
;
; OPTIONAL OUTPUTS:
;   nam        - Name of transition (string)
;   NEWWV=     - Correct wavelength if input wavelength was an
;                'estimate'
;   REF=       - Reference number from my spectral line list
;
; COMMENTS:
;
; EXAMPLES:
;   getfnam, 1215.6701, fval, name
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro getfnam, wave, fval, nam, FIL=fil, REF=ref, CLOSE=close, $
             TOLER=toler, NEWWV=newwv, SILENT=silent

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'getfnam, wave, fval, [nam], FIL=, TOLER= /CLOSE, REF= [v1.1]'
    return
  endif 

  if not keyword_set(TOLER) then begin
      if keyword_set(CLOSE) then toler = 0.2 else toler=0.003
  endif
;
  nwv = n_elements(wave)

  dumc = ' '
  dumi = 0
  dumd1 = 0.0d
  dumd2 = 0.0d

  if not keyword_set(FIL) then $
	  fil = getenv('XIDL_DIR')+'/Spec/Lines/Lists/grb.lst'
  close, 9
  openr, 9, fil, ERROR = err
  if(err NE 0) then stop, 'File does not exist', err

  if nwv EQ 1 then begin  ; one line
      readf, 9, nlin
      for i=1,nlin do begin
          readf, 9, format='(f9.4,a13,f10.5,2x,i2)', dumd1, dumc, dumd2, dumi
          if abs(wave-dumd1) LT toler then begin
              nam = dumc
              fval = dumd2
              ref = dumi
              newwv = dumd1
              break
          endif
          if(i EQ nlin) then begin
              if not keyword_set(SILENT) then print, 'Line isnt in lindat', wave
              fval = -1
              nam = ' '
              return
          end
      endfor

  endif else begin
      ;; Read entire line list first
      readf, 9, nlin, format='(i4)'
      all_wv = dblarr(nlin)
      all_nam = strarr(nlin)
      all_fv = dblarr(nlin)
      all_ref = intarr(nlin)
      for i=0L,nlin-1 do begin
          readf, 9, format='(f9.4,a12,1x,f10.5,2x,i2)', dumd1, dumc, dumd2, $
            dumi
          all_wv[i] = dumd1
          all_nam[i] = dumc
          all_fv[i] = dumd2
          all_ref[i] = dumi
      endfor

      ;; Fill up other arrays
      fval = dblarr(nwv)
      nam = strarr(nwv)
      ref = intarr(nwv)
      newwv = dblarr(nwv)

      for j=0L,nwv-1 do begin
          a = where(abs(wave[j]-all_wv) LT TOLER, na)
          case na of
              0: begin
                  print, 'No line with wave = ', wave[j]
                  nam[j] = '          '
              end
              1: begin
                  nam[j] = all_nam[a]
                  fval[j] = all_fv[a]
                  ref[j] = all_ref[a]
                  newwv[j] = all_wv[a]
              end
              else: stop
          endcase
      endfor
  endelse
  close, 9

  return
end
