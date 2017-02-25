;+ 
; NAME:
; x_setline
;
; PURPOSE:
;    Given an atomic wavelength, return the fvalue and name in the
;    absorption line structure.
;
; CALLING SEQUENCE:
;   linstr = x_setline(wave, LINFIL=)
;
; INPUTS:
;   wave  - rest wavelength
;
; RETURNS:
;  linstr -  Abslinstrct of the relevant line
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  LINFIL=  -- FITS file containing a set of ABSLINSTRCT entries
;  /OVERRIDE -- Bust through the stop sign for lines not in the
;               database
; /CLOSE --  Return the closest line (not exact)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   linstr = x_setline( 1215.6701 )
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------
function x_setline, wave, LINFIL=linfil, CLOSE=close, OVERRIDE=override

 ;
  if  N_params() NE 1  then begin 
    print,'Syntax - ' + $
             'line = x_setline(wave, LINFIL=, /FINE, /CLOSE, /override) [v1.1]'
    return, -1
  endif 
  toler = 0.01
  if keyword_set(FINE) then toler = 1e-5

; Optional keywords

  if not keyword_set( LINFIL ) then $
    linfil = getenv('XIDL_DIR')+'/Spec/Lines/all_lin.fits'

  ;; Grab the data
;  all_lin = xmrdfits( linfil, 1, structyp='abslinstrct', /silent)
  all_lin = xmrdfits( linfil, 1, /silent)
  ;; KLUDGE FOR IDL 8.0 [which won't put f in a binary FITS table!]
  if not tag_exist(all_lin, 'f') and not keyword_set(NOF) then begin
     tmp = { $
           f: 0.d $
           }
     kludge = replicate(tmp, n_elements(all_lin))
     kludge.f = all_lin.fval
     new_lin = struct_addtags(all_lin, kludge)
     all_lin = new_lin
  endif

  nlin = n_elements(wave)
  if nlin GT 1 then begin
      for qq=0L,nlin-1 do begin
          if qq EQ 0 then $
            ret_val = x_setline(wave[qq], LINFIL=linfil, CLOSE=close) $
          else ret_val = [ret_val, $
                          x_setline(wave[qq], LINFIL=linfil, CLOSE=close)]
      endfor
      return, ret_val
  endif

  ;; Search for a match
  a = where(abs(all_lin.wrest-wave[0]) LT TOLER, na)
  case na of 
      0: begin
          if keyword_set(CLOSE) then begin
              mn = min(abs(all_lin.wrest-wave),imn)
              return, all_lin[imn]
          endif else begin
              print, 'x_setline:  No match in database!', wave
              if not keyword_set(OVERRIDE) then stop
              return, -1
          endelse
      end
      1: return, all_lin[a]
      else: begin
          mn = min(abs(all_lin[a].wrest-wave),imn)
          print, 'x_setline:  Multiple hits.  Returning closest: ', $
                 all_lin[a[imn]].ion
          return, all_lin[a[imn]]
      end
  endcase
          
  return, -1
end
