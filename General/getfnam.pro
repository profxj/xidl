;+ 
; NAME:
; getfnam
;
; PURPOSE:
;    Given an atomic wavelength, return the fvalue and name
;
; CALLING SEQUENCE:
;   
;   getfnam, wave, fval, [nam]
;
; INPUTS:
;   wave       - ionic transition
;
; RETURNS:
;   f          - oscillator strength
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   nam        - Name of transition
;
; COMMENTS:
;
; EXAMPLES:
;   getfnam, 1215.6701, fval, name
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro getfnam, wave, fval, nam

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'getfnam, wave, fval, [nam]'
    return
  endif 

;
  nwv = n_elements(wave)

  dumc = ' '
  dumd1 = 0.0d
  dumd2 = 0.0d

  fil = getenv('XIDL_DIR')+'/Spec/Lines/Lists/qal.lst'
  close, 9
  openr, 9, fil, ERROR = err
  if(err NE 0) then stop, 'File does not exist', err

  if nwv EQ 1 then begin  ; one line
      readf, 9, nlin
      for i=1,nlin do begin
          readf, 9, format='(f9.4,a11,2x,f10.5)', dumd1, dumc, dumd2
          if abs(wave-dumd1) LT 0.03 then begin
              nam = dumc
              fval = dumd2
              break
          endif
          if(i EQ nlin) then begin
              print, 'Line isnt in lindat', wave
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
      for i=0L,nlin-1 do begin
          readf, 9, format='(f9.4,a11,2x,f10.5)', dumd1, dumc, dumd2
          all_wv[i] = dumd1
          all_nam[i] = dumc
          all_fv[i] = dumd2
      endfor

      ;; Fill up other arrays
      fval = dblarr(nwv)
      nam = strarr(nwv)

      for j=0L,nwv-1 do begin
          a = where(abs(wave[j]-all_wv) LT 0.03, na)
          case na of
              0: begin
                  print, 'No line with wave = ', wave[j]
                  nam[j] = '          '
              end
              1: begin
                  nam[j] = all_nam[a]
                  fval[j] = all_fv[a]
              end
              else: stop
          endcase
      endfor
  endelse
  close, 9

  return
end
