;+ 
; NAME:
; hires_allflat   
;     Version 1.1
;
; PURPOSE:
;    Identifies flats, processes them.  Creates one flat file
;      per slit width.  Takes TFLAT as the default
;
; CALLING SEQUENCE:
;  hires_allflat, hires, setup
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setups  (can be an array of values)
;   [chip]  -  Chip to process [default: 1,2,3]
;
; RETURNS:
;
; OUTPUTS:
;  /NOGAIN -- Do not bother to calculate the gain
;  /TFCHK  -- Check the trace flat
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_allflat, hires, [1L,2L,3L,4L]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Feb-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_allflat, hires, setup, chip, CLOBBER=clobber, $
                  CHK=chk, TFCHK=tfchk, _EXTRA=extra, NOGAIN=nogain

;
 if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_allflat, hires, setup, [chip], /clobber, /NOGAIN [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( CHIP ) then begin
     if strmatch(hires[0].ccd,'SINGLE') then chip = [-1L] else $
        chip = [1L,2L,3L]
  endif
;  if not keyword_set( iflat ) then iflat = 0L
 
  for jj=0L,n_elements(setup)-1 do begin
      
      ;; Find the gain
     if not keyword_set(NOGAIN) and  (~ strmatch(hires[0].ccd,'SINGLE')) then $
        hires_findgain, hires, setup[jj]

      for i=0, n_elements(chip) -1 do begin
          qq = chip[i]
 
          ;; COMBINE
          hires_mktflat, hires, setup[jj], qq, CLOBBER=clobber, _EXTRA=extra

          ;; TRACE
          hires_edgeflat, hires, setup[jj], qq, CLOBBER=clobber, $
            INTER=inter, CHK=tfchk, _EXTRA=EXTRA

          ;; Check for Pixel flat
          indx = where(hires.chip EQ qq AND $
                       hires.flg_anly NE 0L AND $
                       hires.setup EQ setup[jj], nindx)
          flg_f = hires_getpixflat(hires[indx[0]], FIL=flatfil)

          ;; Normalize
          if flg_f EQ -1 OR keyword_set(NORMTRC) then $
            hires_nrmflat, hires, setup[jj], qq, CLOBBER=clobber, _EXTRA=EXTRA
      endfor
  endfor

  print, 'hires_allflat: All done!'
  return
end
