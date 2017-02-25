;+ 
; uves_chktflat   
;     Version 1.1
;
; PURPOSE:
;   Give a visual assessment of the order traces.  This program works
;   in coordination with xatv.  An assessment of the 2D fit is given
;   by turning on the /FIT keyword.  
;
; CALLING SEQUENCE:
;   
;  mike_chktrcflat, mike, setup, chip, /NOSTOP, /FIT, XOFF=, INDX=
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
;  INDX  -- Determine the xoffset from the uves structure
;           (e.g. uves[indx].arc_xyoff[0])
;  /FIT  -- Plot the 2D fits instead of the individual traces.
;  /NOSTOP -- Run without giving a warning to first run xatv
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_chktrcflat, uves, 1, 1, /nostop, /fit
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_getfil
;
; REVISION HISTORY:
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uves_chktflat, uves, setup, side, NOSTOP=nostop, FIT=fit, $
  XOFF=xoff, INDX=indx
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'uves_chktflat, uves, setup, side, /NOSTOP, /FIT, XOFF=, INDX= [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( XOFF ) then xoff = 0.
  if keyword_set( INDX ) then xoff = uves[indx].arc_xyoff[0]

; Loop on chip
      ;; SIDE
   case side of 
       1: print, 'uves_chktrcflat: Tracing BLUE flat'
       2: print, 'uves_chktrcflat: Tracing Red flat'
       else: stop
   endcase

   ;; Warn about xatv!
   ix = where(uves.setup EQ setup and uves.side EQ side)
   wcen = uves[ix[0]].xdangl
   flatfil = uves_getfil('qtz_fil', setup, WCEN=wcen, /name)
   if not keyword_set( NOSTOP ) then begin
       print, 'uves_chktrcflat:  You need to run this command before '+ $
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
      trc_str = uves_getfil('tflat_str', setup, WCEN=wcen)
      
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
      ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen)
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
