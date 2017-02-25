;+ 
; NAME:
; dla_xalphvsa
;   Version 1.1
;
; PURPOSE:
;    Plots [X/alpha] vs [alpha/H].  Currently, the fig only shows Echelle obs
;
; CALLING SEQUENCE:
;   
;   dla_xyvsab, X, XR=, YR=, DLA=
;
; INPUTS:
;   X -- Atomic number of element X
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
;   dla_xalphavsa, 7
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_xalphvsa, xatm, XR=xr, YR=yr, DLA=dla, PSFILE=psfile, NOLIM=nolim, $
                NOXERR=noxerr, NOCLOSE=noclose, MXERR=mxerr

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'dla_xalphvsa, xatm, XR=, YR=, DLA=, /NOLIM, /NOXERR, ' + $
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
  getabnd, xlbl, xatm, xabnd, /flag
  
  xi = xatm

  ytit = '['+xlbl+'/!9a!X]'
  xtit = '[!9a!X/H]'

  ;; Full set
  gd = where( dla.XH[xi].flgclm NE 0 AND dla.flgalpha EQ 1 AND $
              (dla.XH[xi].flginst MOD 4 LT 2), ngd)

  if ngd EQ 0 then return
  xyval = dla[gd].XH[xi].clm - dla[gd].alpha
  abval = dla[gd].alpha

  plot, abval, xyval, color=clr.black, background=clr.white, charsize=csize,$
    xmargin=xmrg, ymargin=ymrg, xtitle=xtit, ytitle=ytit, $
    /nodata, xstyle=1, ystyle=1, xrange=xr, yrange=yr

  ;; Loop  
  for qq=0L,ngd-1 do begin
      jj=gd[qq]
      ;; Plot XY
      case dla[jj].XH[xi].flgclm of
          1: begin
              sigxy = sqrt(dla[jj].XH[xi].sigclm^2 + dla[jj].sigalpha^2)
              if sigxy LT MXERR then $
                oploterror, [abval[qq]], [xyval[qq]], sigxy, errcolor=clr.black, $
                color=clr.black
          end
          2: begin
              if not keyword_set(NOLIM) then begin
                  ;; Lower limit
                  plotsym, 2, 3.5, thick=2.
                  oplot, [abval[qq]], [xyval[qq]], psym=8, color=clr.red
              endif 
          end
          3: begin
              if not keyword_set(NOLIM) then begin
                  ;; Upper limit
                  plotsym, 1, 3.5, thick=2.
                  oplot, [abval[qq]], [xyval[qq]], psym=8, color=clr.blue
              endif
          end
          else: stop
      endcase
      ;; Plot AB
      if not keyword_set( NOXERR ) then begin
          sigab = dla[jj].sigalpha
          oploterror, [abval[qq]], [xyval[qq]], sigab, [0.], $
            errcolor=clr.black, color=clr.black
      endif else oplot, [abval[qq]], [xyval[qq]], color=clr.black
  endfor

  if not keyword_set( NOCLOSE ) then begin
      if keyword_set( PSFILE ) then x_psclose else x_tiffset, /unset
  endif
  return
end
