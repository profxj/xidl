;+ 
; NAME:
; cldy_uplot
;   Version 1.0
;
; PURPOSE:
;    Creates a Uplot for a given NHI, FeH, nH
;
; CALLING SEQUENCE:
;   
; cldy_Uplot, grid, NHI, FeH, nH, ions
;
; INPUTS:
;   grid  - CLOUDY grid
;   NHI - Can be an array of values
;   FeH
;   nH
;   ions  - Array of [Z,ion] vectors
;
; RETURNS:
;   
;
; OUTPUTS:
;   Creates a Plot
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_Uplot, grid, 19.0d, -1.0d, -1.0d, [[14,2], [14,3]]
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_uplot, grid, NHI, FeH, nH, ions, YMNX=ymnx, PSFIL=psfil

;
  if  N_params() LT 5  then begin 
      print, 'Syntax - ' +$
        'cldy_Uplot, grid, NHI, FeH, nH, ions, YMNX=, PSFIL= (v1.0)'
      return
  endif 

; Optional keywords

  if not keyword_set( YMNX ) then ymnx = [8.0, 18.0]


;  Loop on Plot
  n_NHI = n_elements(NHI)
;  loadct, 2
  
  for ww=0, n_NHI-1 do begin
;  Find Grid Subset

  ii = where(grid.NHI EQ NHI[ww] AND grid.FeH EQ FeH AND grid.nH EQ nH AND $
             grid.flg NE 0, count)

  if count LE 1 then return

  xmin = min(grid[ii].U, MAX=xmax)
  xplt = grid[ii].U

; Loop on the ions

  sz = size(ions)
  nions = sz[2]
  dumx = fltarr(2)
  dumy = fltarr(2)
  dumx[0] = -99.
  dumy[0] = -99.

; PS FILE
  if keyword_set(PSFIL) then begin
      device, decompose=0
      !p.thick = 5
      !p.charthick = 4
      ps_open, file=psfil, font=1, /color, bpp=8
      !y.margin = [5,2]
  endif

; Set colors
  !P.MULTI = [0, 1, n_NHI]
  xcolors = x_setclrs()
  clr = getcolor(/load)

  for q = 0, nions-1 do begin
;   Abund

      getabnd, elm, ions[0,q], abnd, flag=1
      yplt = NHI[ww] - grid[ii].X[1,1] - 12. + abnd + $
        grid[ii].X[ions[0,q],ions[1,q]]
      ; Metallicity
      if ions[0,q] NE 1 then yplt = yplt + FeH
      ;
      if q EQ 2 then psym = 1 else psym=q+1
      psym = psym<7

      if q EQ 0 then plot, xplt, yplt, psym=psym, xrange=[xmin,xmax], $
        yrange=[ymnx[0],ymnx[1]], color=clr.black, $
        background=clr.white, ystyle=1, xstyle=1, $
	xtitle='log U', ytitle='N(X)' $
      else oplot, xplt, yplt, psym=psym, color=xcolors[q]

;      Label
      case ions[1,q] of
          1: nm = elm+'I'
          2: nm = elm+'II'
          3: nm = elm+'III'
          4: nm = elm+'IV'
          5: nm = elm+'V'
          6: nm = elm+'VI'
          7: nm = elm+'VII'
          else: stop
      endcase

      dumx[1] = xmax-0.5
      dumy[1] = ymnx[0]+3-q*0.3
      oplot, dumx, dumy, psym=psym, color=xcolors[q]
      xyouts, xmax-0.4, ymnx[0]+2.9-q*0.3, nm, charsize=1.5, $
        color=xcolors[q]
      xyouts, xmin+0.2, ymnx[1]-0.3, 'N(HI)='+string(NHI[ww]), $
        charsize=1.5, color=clr.black

  endfor
endfor

  if keyword_set(PSFIL) then begin
      ps_close, /noprint, /noid
      device, decomposed=1
      !p.thick = 1
      !p.charthick = 1
      !y.margin = [4,2]
  endif
           
  !P.MULTI = [0, 1, 1]
return
end

