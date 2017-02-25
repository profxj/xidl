;+ 
; NAME:
; mike_setgain
;     Version 2.0
;
; PURPOSE:
;    Calculate the inverse gain using the routine mike_gain and the Milky Flats.
;   The value of the gain is then saved into the tag mike.gain.
;   This routine calls  x_calcagin.
;
; CALLING SEQUENCE:
;   
;  mike_setgain, mike, setup, [side]
;
; INPUTS:
;   setup   -  Setup ID
;   [side]  -  Blue (1) and/or Red (2) side  [default: [1,2] ]
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_calcgain
;
; REVISION HISTORY:
;   26-Jun-2004 Written by SMB
;   04-Apr-2005 Replaced mike_gain with x_calcgain
;-
;------------------------------------------------------------------------------

pro mike_setgain, mike, setup, side

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_setgain, mike, setup, [side] [v2.0]'
      return
  endif 
  
  if not keyword_set( SIDE ) then side = [1L,2L]

  ;; Loop on side
  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]
      if qq EQ 1 then print, 'mike_setgain: Setting the gain for the BLUE side' $
      else print, 'mike_setgain: Setting the gain for the RED side'

      ;; Indices
      mflt = where(mike.side EQ qq AND mike.flg_anly NE 0 AND $
                    strtrim(mike.type,2) EQ 'MFLT' AND $
                    mike.setup EQ setup, nflt)


      ;; Bias sub first image
      rawfil = mike[mflt[0]].rootpth+mike[mflt[0]].img_root
      rslt = mike_subbias_sngl(rawfil, qq, OVIMG=img1, /NOFITS, /SILENT)

      ;; Start looping
      svgain = fltarr(nflt-1)
      for jj=1L,nflt-1 do begin
          ;; 
          rawfil = mike[mflt[jj]].rootpth+mike[mflt[jj]].img_root
          rslt = mike_subbias_sngl(rawfil, qq, OVIMG=img2, /NOFITS, /SILENT)
          ;; Gain
          gain = x_calcgain(img1, img2)
          ;; Save 
          svgain[jj-1] = gain
          print, 'mike_setgain: ', jj, gain
          img1 = img2
      endfor

      ;; Grab the value
      djs_iterstat, svgain, sigrej=3., mean=gdgain
      
      ;; Update all instances
      all_side = where(mike.side EQ qq AND mike.setup EQ setup)
      print, 'mike_setgain: Applying gain of ', gdgain
      mike[all_side].gain = gdgain

  endfor

  print, 'mike_setgain: All done!'

  return
          

end

