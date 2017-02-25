;+ 
; hires_chktrcflat   
;     Version 1.1
;
; PURPOSE:
;   Give a visual assessment of the order traces.  This program works
;   in coordination with xatv.  An assessment of the 2D fit is given
;   by turning on the /FIT keyword.  
;
; CALLING SEQUENCE:
;   
;  hires_chktrcflat, hires, setup, chip, /NOSTOP, /FIT, XOFF=, INDX=
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   chip     -  Blue (1) or Red (2)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  XOFF -- Offset the fits in the x-direction.  This is very useful
;          for correcting the thermal gradient.
;  INDX  -- Determine the xoffset from the hires structure
;           (e.g. hires[indx].arc_xyoff[0])
;  /FIT  -- Plot the 2D fits instead of the individual traces.
;  /NOSTOP -- Run without giving a warning to first run xatv
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_chktrcflat, hires, 1, 1, /nostop, /fit
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_getfil
;
; REVISION HISTORY:
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_chktrcflat, hires, setup, chip, NOSTOP=nostop, FIT=fit, $
  XOFF=xoff, INDX=indx
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_chktrcflat, hires, setup, chip, /NOSTOP, /FIT, XOFF=, INDX= [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( XOFF ) then xoff = 0.
  if n_elements(INDX) GT 0 then xoff = hires[indx].arc_xyoff[0]

; Loop on chip
      ;; SIDE
   case chip of 
       -1: print, 'hires_chktrcflat: Single chip'
       1: print, 'hires_chktrcflat: Tracing BLUE flat'
       2: print, 'hires_chktrcflat: Tracing GREEN flat'
       3: print, 'hires_chktrcflat: Tracing RED flat'
       else: stop
   endcase

   ;; Warn about xatv!
   flatfil = hires_getfil('qtz_fil', setup, CHIP=chip, /name)
   if not keyword_set( NOSTOP ) then begin
       print, 'hires_chktrcflat:  You need to run this command before '+ $
       'running this program:  '
       print, '  IDL>  xatv, ''', flatfil, ''''
       print, 'If you havent done this yet, return out of this program'
       print, '  run the command and rerun this program with /NOSTOP'
       print, 'If you have done this command, just continue'
       stop
   endif

  xatverase


  if not keyword_set( FIT ) then begin 
      ;; Open trace structure
      trc_str = hires_getfil('tflat_str', setup, CHIP=chip)
      
      ;; Loop on order
      sz = size(trc_str.xcen, /dimensions)
      if n_elements(xoff) EQ 1 then xoff = replicate(xoff, sz[1])
      ;; PLOT
      for q=0L, sz[1]-1 do begin
          a = where(trc_str.xerr[*,q] LT 90., na)
          if na NE 0 then $
            xatvplot, trc_str.xcen[a,q] + xoff[q], [float(a)]
      endfor
  endif else begin
      ;; Order structure
      ordr_str = hires_getfil('ordr_str', setup, CHIP=chip)
      dumy = findgen(n_elements(ordr_str[0].lhedg))
      if n_elements(xoff) EQ 1 then $
        xoff = replicate(xoff, n_elements(ordr_str))
      ;; Plot
      for q=0L, n_elements(ordr_str)-1 do begin
          xatvplot, ordr_str[q].lhedg + xoff[q], dumy
          xatvplot, ordr_str[q].rhedg + xoff[q], dumy, color=2
      endfor
  endelse

  return
end
