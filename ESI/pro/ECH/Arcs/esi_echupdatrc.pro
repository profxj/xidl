;+ 
; NAME:
; esi_echupdatrc   
;     Version 1.0
;
; PURPOSE:
;    Takes the arc, creates a fit and traces the lines.
;
; CALLING SEQUENCE:
;   
;  esi_echupdatrc, esi, slit, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echupdatrc, esi, slit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echupdatrc, slit, shift, GUESSARC=guessarc

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echupdatrc, slit, shift [v1.0]'
      return
  endif 
  
;  Optional Keywords
  c_s = esi_slitnm(slit)
  
; Open Guess Arc
  gdx_fil = getenv('XIDL_DIR')+'/ESI/CALIBS/ArcECH_'+c_s+'gdx.fits'
  all_gdx = mrdfits(gdx_fil, /silent)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Loop on TRC
  for qq=0,9L do begin
      ;; Open trcstr
      ordr = 15L - qq
      if ordr LT 10 then cordr = '0'+string(ordr, FORMAT='(i1)') $
      else cordr = string(ordr, FORMAT='(i2)')
      infil = 'Arcs/TRC/ArcECH_'+c_s+'trc'+cordr+'.fits'
      trcstr = mrdfits(infil, 1, /silent)

      ;; Grab the good ones from the current file
      gdpk = all_gdx[where(all_gdx[*,qq] GT 0,ngd),qq]
      newpk = gdpk

      ;; Loop on the new ones
      npk = n_elements(trcstr.xstrt)
      for i=0L,npk-1 do begin
          mn = min(abs(trcstr.xstrt[i] - gdpk - shift))
          if abs(mn) GT 4 and trcstr.xstrt[i] NE 0. then begin
              newpk = [newpk, trcstr.xstrt[i]-shift]
              print, 'esi_echupdarc: Adding ...', trcstr.xstrt[i]-shift, mn
          endif
      endfor
      nnew = n_elements(newpk)
      ;; Write to all_gdx
      print, 'esi_echupdatrc: Added ', strtrim(nnew-ngd,2), $
        ' points to order '+cordr
      all_gdx[0:nnew-1,qq] = newpk
  endfor

  ;; Output
  stop
  print, 'esi_echupdatrc: Overwriting ', gdx_fil
  mwrfits, all_gdx, gdx_fil, /create, /silent
  ;;
  print, 'esi_echupdatrc: All done!'
      
  return
end
