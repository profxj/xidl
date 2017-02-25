;+ 
; NAME:
; dla_xyvsab
;   Version 1.1
;
; PURPOSE:
;    Plots [X/Y] vs [A/B].  Currently, the fig only shows Echelle obs
;
; CALLING SEQUENCE:
;   
;   dla_xyvsab, [X,Y], [A,B], XR=, YR=, DLA=
;
; INPUTS:
;   X -- Atomic number of element X
;   Y -- Atomic number of element Y
;   A -- Atomic number of element A
;   B -- Atomic number of element B
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  XR= -- x-range
;  YR= -- y-range
;  DLA= -- DLA structure
;  /NOXERR -- Suppress x error bars
;  /NOLIM  -- Do not show limits
;  /NOCLOSE -- Do not close the plot
;  MXERR -- Maximum error to plot 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_xyvsab, [14,26], [14,1]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dla_xyvsab, xyatm, abatm, XR=xr, YR=yr, DLA=dla, PSFILE=psfile, NOLIM=nolim, $
                NOXERR=noxerr, NOCLOSE=noclose, MXERR=mxerr

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'dla_xyvsab, [x,y], [a,b], XR=, YR=, DLA=, /NOLIM, /NOXERR, ' + $
      '/NOCLOSE, MXERR= [v1.1]'
    return
  endif 

  ;; Optional Keywords
  if not keyword_set( MXERR ) then mxerr = 99.9

  ;; FILE
  if not keyword_set( DLA ) then parse_dlalst, dla, /NOHIS
  dla = dla[where(dla.NHI GE 20.3)]
  
  if keyword_set( PSFILE ) then begin
      x_psopen, psfile, /maxs
      xmrg = [10.5,1]
      ymrg = [5,1]
      csize = 2.0
  endif else begin
      x_tiffset
      xmrg = [9,1]
      ymrg = [5,1]
      csize = 2.5
  endelse
      
  clr = getcolor(/load)

  ;; Labels
  getabnd, xlbl, xyatm[0], xabnd, /flag
  getabnd, ylbl, xyatm[1], yabnd, /flag
  getabnd, albl, abatm[0], aabnd, /flag
  getabnd, blbl, abatm[1], babnd, /flag
  
  xi = xyatm[0]
  yi = xyatm[1]
  ai = abatm[0]
  bi = abatm[1]

  ytit = '['+xlbl+'/'+ylbl+']'
  xtit = '['+albl+'/'+blbl+']'

  ;; Full set
  gd = where( (dla.XH[xi].flgclm * dla.XH[yi].flgclm * $
             dla.XH[ai].flgclm * dla.XH[bi].flgclm NE 0) AND $
              dla.XH[xi].flginst MOD 4 LT 2 AND $
              dla.XH[yi].flginst MOD 4 LT 2 AND $
              dla.XH[ai].flginst MOD 4 LT 2 AND $
              dla.XH[bi].flginst MOD 4 LT 2, ngd)

  if ngd EQ 0 then return
  xyval = dla[gd].XH[xi].clm - dla[gd].XH[yi].clm
  abval = dla[gd].XH[ai].clm - dla[gd].XH[bi].clm

  plot, abval, xyval, color=clr.black, background=clr.white, charsize=csize,$
    xmargin=xmrg, ymargin=ymrg, xtitle=xtit, ytitle=ytit, $
    /nodata, xstyle=1, ystyle=1, xrange=xr, yrange=yr

  ;; Loop  
  for qq=0L,ngd-1 do begin
      jj=gd[qq]
      ;; Plot XY
      if dla[jj].XH[xi].flgclm EQ 1 AND dla[jj].XH[yi].flgclm EQ 1 then begin
          sigxy = sqrt(dla[jj].XH[xi].sigclm^2 + dla[jj].XH[yi].sigclm^2)
          if sigxy LT MXERR then $
            oploterror, [abval[qq]], [xyval[qq]], sigxy, errcolor=clr.black, $
            color=clr.black
      endif else begin
          if not keyword_set(NOLIM) then begin
              if (dla[jj].XH[xi].flgclm EQ 2 AND dla[jj].XH[yi].flgclm NE 2) $
                OR (dla[jj].XH[xi].flgclm NE 3 AND dla[jj].XH[yi].flgclm EQ 3) $
                then begin
                  ;; Lower limit
                  plotsym, 2, 3.5, thick=2.
                  oplot, [abval[qq]], [xyval[qq]], psym=8, color=clr.red
              endif 
              if (dla[jj].XH[xi].flgclm EQ 3 AND dla[jj].XH[yi].flgclm NE 3) $
                OR (dla[jj].XH[xi].flgclm NE 2 AND dla[jj].XH[yi].flgclm EQ 2) $
                then begin
                  ;; Upper limit
                  plotsym, 1, 3.5, thick=2.
                  oplot, [abval[qq]], [xyval[qq]], psym=8, color=clr.blue
              endif
          endif
      endelse

      ;; Plot AB
      if dla[jj].XH[ai].flgclm EQ 1 AND dla[jj].XH[bi].flgclm EQ 1 then begin
          if not keyword_set( NOXERR ) then begin
              sigxy = sqrt(dla[jj].XH[ai].sigclm^2 + dla[jj].XH[bi].sigclm^2)
              oploterror, [abval[qq]], [xyval[qq]], sigxy, [0.], $
                errcolor=clr.black, color=clr.black
          endif else oplot, [abval[qq]], [xyval[qq]], color=clr.black
      endif else begin
          if not keyword_set(NOLIM) then begin
              if (dla[jj].XH[ai].flgclm EQ 2 AND dla[jj].XH[bi].flgclm NE 2) $
                OR (dla[jj].XH[ai].flgclm NE 3 AND dla[jj].XH[bi].flgclm EQ 3) $
                then begin
                  ;; Lower limit
                  plotsym, 7, 3.5, thick=2.
                  oplot, [abval[qq]], [xyval[qq]], psym=8, color=clr.red
              endif 
              if (dla[jj].XH[ai].flgclm EQ 3 AND dla[jj].XH[bi].flgclm NE 3) $
                OR (dla[jj].XH[ai].flgclm NE 2 AND dla[jj].XH[bi].flgclm EQ 2) $
                then begin
                  ;; Upper limit
                  plotsym, 6, 3.5, thick=2.
                  oplot, [abval[qq]], [xyval[qq]], psym=8, color=clr.blue
              endif
          endif
      endelse
  endfor

  if not keyword_set( NOCLOSE ) then begin
      if keyword_set( PSFILE ) then x_psclose else x_tiffset, /unset
  endif
  return
end
