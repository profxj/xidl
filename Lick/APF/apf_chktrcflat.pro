;+ 
; apf_chktrcflat   
;     Version 1.1
;
; PURPOSE:
;   Give a visual assessment of the order traces.  This program works
;   in coordination with xatv.  An assessment of the 2D fit is given
;   by turning on the /FIT keyword.  
;
; CALLING SEQUENCE:
;   
;  apf_chktrcflat, apf, setup, chip, /NOSTOP, /FIT, XOFF=, INDX=
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
;  INDX  -- Determine the xoffset from the apf structure
;           (e.g. apf[indx].arc_xyoff[0])
;  /FIT  -- Plot the 2D fits instead of the individual traces.
;  /NOSTOP -- Run without giving a warning to first run xatv
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   apf_chktrcflat, apf, 1, 1, /nostop, /fit
;
; PROCEDURES/FUNCTIONS CALLED:
;  apf_getfil
;
; REVISION HISTORY:
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro apf_chktrcflat, apf, setup, NOSTOP=nostop, FIT=fit, $
  XOFF=xoff, INDX=indx
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_chktrcflat, apf, setup, /NOSTOP, /FIT, XOFF=, INDX= [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( XOFF ) then xoff = 0.
  if n_elements(INDX) GT 0 then xoff = apf[indx].arc_xyoff[0]

   ;; Warn about xatv!
   flatfil = apf_getfil('qtz_fil', setup, /name)
   if not keyword_set( NOSTOP ) then begin
       print, 'apf_chktrcflat:  You need to run this command before '+ $
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
      trc_str = apf_getfil('tflat_str', setup)
      
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
      ordr_str = apf_getfil('ordr_str', setup)
      dumy = findgen(n_elements(ordr_str[0].lhedg))
      if n_elements(xoff) EQ 1 then $
        xoff = replicate(xoff, n_elements(ordr_str))
      ;; Plot
      for q=0L, n_elements(ordr_str)-1 do begin
         if keyword_set(ANLY) and ordr_str[q].flg_anly EQ 0 then continue
          xatvplot, ordr_str[q].lhedg + xoff[q], dumy
          xatvplot, ordr_str[q].rhedg + xoff[q], dumy, color=2
      endfor
  endelse

  return
end
