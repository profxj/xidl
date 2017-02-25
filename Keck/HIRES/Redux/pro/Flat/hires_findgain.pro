;+ 
; NAME:
; hires_findgain
;     Version 1.1
;
; PURPOSE:
;    Calculate the inverse gain using the routine hires_gain. 
;   The value of the gain is then saved into hires.gain
;
; CALLING SEQUENCE:
;   
;  hires_findgain, hires, setup, [side]
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setup ID
;   [chip]  -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
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
;
; REVISION HISTORY:
;   26-Jun-2004 Written by SMB
;-
;------------------------------------------------------------------------------

pro hires_findgain, hires, setup, chip

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_findgain, hires, setup, [chip] [v1.1]'
      return
  endif 
  
  if not keyword_set( CHIP ) then chip = [1L,2L,3L]

  hires_rdxlog, 'flat_log', 'hires_findgain', /ophdr
  ;; Loop on side
  for ii=0L,n_elements(chip)-1 do begin
      qq = chip[ii]
      case qq of 
          1: mssg =  'hires_findgain: Finding the gain for the BLUE CCD' 
          2: mssg =  'hires_findgain: Finding the gain for the GREEN CCD' 
          3: mssg =  'hires_findgain: Finding the gain for the RED CCD' 
          else: stop
      endcase
      hires_rdxlog, 'flat_log', mssg

      ;; Indices
      mflt = where(hires.chip EQ qq AND hires.flg_anly NE 0 AND $
                    strtrim(hires.type,2) EQ 'PFLT' AND $
                    hires.setup EQ setup, nflt)

      if nflt EQ 0 then begin
          flg_flat = 1
          print, 'hires_findgain:  No Milky Flats, using Traces'
          mflt = where(hires.chip EQ qq AND hires.flg_anly NE 0 AND $
                       strtrim(hires.type,2) EQ 'TFLT' AND $
                       hires.setup EQ setup, nflt)
      endif else flg_flat = 0


      ;; Bias sub first image
      rawfil = hires[mflt[0]].rootpth+hires[mflt[0]].img_root
      rslt = hires_subbias_sngl(rawfil, qq, qq, $
                                OVIMG=img1, /NOFITS, /SILENT, $
                                FRAME=hires[mflt[0]].frame)

      ;; Start looping
      if nflt LE 1 then begin
          print, 'hires_findgain: Not enough flats to determine gain'
          continue
      endif
      svgain = fltarr(nflt-1)
      for jj=1L,nflt-1 do begin
          ;; 
          rawfil = hires[mflt[jj]].rootpth+hires[mflt[jj]].img_root
          rslt = hires_subbias_sngl(rawfil, qq, qq, OVIMG=img2, $
                                    FRAME=hires[mflt[0]].frame, $
                                    /NOFITS, /SILENT)
          ;; Gain
          gain = x_calcgain(img1, img2)
          ;; Save 
          svgain[jj-1] = gain
          print, 'hires_findgain: ', jj, gain
          img1 = img2
      endfor
      
      ;; Grab the value
      djs_iterstat, svgain, sigrej=3., mean=gdgain
      
      ;; Update all instances
      all_side = where(hires.chip EQ qq and hires.setup EQ setup)
      mssg =  'hires_findgain: Applying gain of '+string(gdgain)
      hires_rdxlog, 'flat_log', mssg
      hires[all_side].gain = gdgain
      
  endfor

  ;; Log
  mssg= 'hires_findgain: All done!'
  hires_rdxlog, 'flat_log', mssg, /clhdr

  return
          

end

