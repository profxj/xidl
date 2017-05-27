function deimos_adu2e, adu, ampno, single=single, quick=quick
;+
; NAME:
;    deimos_adu2e
;
; PURPOSE:
;    converts adu to electrons based on nonlinear response of DEIMOS CCD mosaic
;
; CALLING SEQUENCE:
;    deimos_adu2e(adu,ampno)
;
; INPUTS:
;     adu-- input signal
;     ampno (1-16) which DEIMOS amplifier?
;
; OUTPUTS:
;    returns electron signal, corrected for gain and nonlinear response
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  uses system variables defined by DEIMOS_RED_MOSAIC

; REVISION HISTORY:
;   13apr02  md
;----------------------------------------------------------------------

if n_elements(quick) eq 0 then quick = 0

if ampno lt 1 or ampno gt 16 then begin
   print, 'amplifier number out of range'
   return,0
endif

if n_elements(single) eq 0 then single =  1

ii=ampno-1

if quick gt 0 then begin
   invgain= [1.23, 1.24, 1.23, 1.20, 1.20, 1.27, 1.25, 1.27, $
          1.26, 1.26, 1.22, 1.25, 1.40, 1.25, 1.25, 1.25 ]
   return, adu*invgain[ii]
endif else begin
   if single then $
      electrons= spl_interp(!deimos_adu1[*,ii],!deimos_electrons1[*,ii],$
         !deimos_gsset1[*,ii],adu) $
    else $
       electrons= spl_interp(!deimos_adu2[*,ii],!deimos_electrons2[*,ii],$
         !deimos_gsset2[*,ii],adu) 
endelse

return, electrons
end


