;+ 
; NAME:
; mike_chktrcflat   
;     Version 2.0
;
; PURPOSE:
;   Give a visual assessment of the order traces.  This program works
;   in coordination with xatv.  An assessment of the 2D fit is given
;   by turning on the /FIT keyword.  
;
; CALLING SEQUENCE:
;   
;  mike_chktrcflat, mike, setup, side, /NOSTOP, /FIT, XOFF=, INDX=
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   side     -  Blue (1) or Red (2)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  XOFF -- Offset the fits in the x-direction.  This is very useful
;          for correcting the thermal gradient.
;  INDX  -- Determine the xoffset from the mike structure
;           (e.g. mike[indx].arc_xyoff[0])
;  /FIT  -- Plot the 2D fits instead of the individual traces
;           (recommended)
;  /NOSTOP -- Run without giving a warning to first run xatv
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_chktrcflat, mike, 1, 1, /nostop, /fit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_getfil
;
; REVISION HISTORY:
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_chktrcflat, mike, setup, side, NOSTOP=nostop, FIT=fit, XOFF=xoff, $
                     INDX=indx
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_chktrcflat, mike, setup, side, /NOSTOP, /FIT, XOFF=, INDX= [v2.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( XOFF ) then xoff = 0.
  if keyword_set( INDX ) then xoff = mike[indx].arc_xyoff[0]

; Loop on side
      ;; SIDE
   case side of 
       1: print, 'mike_chktrcflat: Tracing BLUE flat'
       2: print, 'mike_chktrcflat: Tracing RED flat'
       else: stop
   endcase

   ;; Warn about xatv!
   flatfil = mike_getfil('tflat_fil', setup, SIDE=side, /name)
   if not keyword_set( NOSTOP ) then begin
       print, 'mike_chktrcflat:  You need to run this command before '+ $
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
      trc_str = mike_getfil('tflat_str', setup, SIDE=side)
      
      ;; Loop on order
      sz = size(trc_str.xcen, /dimensions)
      ;; PLOT
      for q=0L, sz[1]-1 do begin
          a = where(trc_str.xerr[*,q] LT 90., na)
          if na NE 0 then $
            xatvplot, trc_str.xcen[a,q] + xoff, [float(a)]
      endfor
  endif else begin
      ;; Order structure
      ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
      dumy = findgen(n_elements(ordr_str[0].lhedg))
      ;; Plot
      for q=0L, n_elements(ordr_str)-1 do begin
          xatvplot, ordr_str[q].lhedg + xoff, dumy
          xatvplot, ordr_str[q].rhedg + xoff, dumy
      endfor
  endelse

  return
end
