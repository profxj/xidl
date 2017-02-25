;+ 
; NAME:
; dla_eplot   
;   Version 1.1
;
; PURPOSE:
;    Given an abundance file and NHI value, create a abundance
;   number plot (logarithmic with H = 12)
;
; CALLING SEQUENCE:
;   
;   dla_eplot, fil, NHI, ZMTL=, PSFILE=, XR=, YR=
;
; INPUTS:
;   fil -- Abundance file, [format: Zatm, Ncolm, Nsig, Flag, Instr]
;   NHI -- NHI value of the DLA
;
; RETURNS:
;
; OUTPUTS:
;   If PSFILE is set, will create a ps file instead of 
;       plotting to the screen
;
; OPTIONAL KEYWORDS:
;  ZMTL -- Metallicity of the gas (for overplotting Solar pattern)
;  XR   -- X range of the plot  [default: min and max of Zatm values]
;  YR   -- Y range of the plot  [default: [4., 12] ]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_eplot, fil, 20.8, ZMTL=-1.2
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  readcol
;
; REVISION HISTORY:
;   17-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_eplot, fil, NHI, ZMTL=zmtl, PSFILE=psfile, XR=xr, YR=yr

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'dla_eplot, fil, NHI, ZMTL=, PSFILE=, XR=, YR= [v1.1]'
    return
  endif 

; Optional Keywords
;  if not keyword_set( NHI ) then NHI = 21.
  if not keyword_set( YR ) then yr = [4., 12.]

; FILE
  readcol, fil, znum, nc, nsig, flg, instr, format='i,f,f,i,i'
  
; DLA
  nelm = n_elements(znum) 

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  clr = getcolor(/load)

  if not keyword_set( XR ) then begin
      xr = fltarr(2)
      xr[0] = min(znum)
      xr[1] = max(znum)
  endif


  plot, xr, yr, color=clr.black, background=clr.white, charsize=2.5,$
    xmargin=[6,1], ymargin=[3,1], xtitle='Atomic Number', ytitle='!9e!X(X)', $
    /nodata, xstyle=1, ystyle=1

  ;; Loop  
  solabd = fltarr(nelm)
  for qq=0L,nelm-1 do begin
      ;; Grab abund
      getabnd, nm, znum[qq], abnd, flag=1
      solabd[qq] = abnd
      ;; Plot
      val = nc[qq] - NHI + 12. 
      case flg[qq] of
          1: oploterror, [znum[qq]], [val], nsig[qq], errcolor=clr.black, $
            color=clr.black
          2: begin ;; Lower limit
              plotsym, 2, 3.5, thick=2.
              oplot, [znum[qq]], [val], psym=8, color=clr.red
          end
          3: begin ;; Upper limit
              plotsym, 1, 3.5, thick=2.
              oplot, [znum[qq]], [val], psym=8, color=clr.blue
          end
          else: stop
      endcase
      xyouts, znum[qq]+0.5, val, nm, color=clr.black, charsize=1.5
  endfor

  ;; Solar
  if keyword_set( ZMTL ) then $
    oplot, znum, solabd+ZMTL, color=clr.green, linestyle=1
      
  if keyword_set( PSFILE ) then x_psclose

end
