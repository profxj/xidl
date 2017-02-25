;+ 
; NAME:
; cldy_uplot
;   Version 1.0
;
; PURPOSE:
;    Creates a 'Uplot' for given NHI, FeH, nH values after inputting
;   a Cloudy grid and a set of ions.
;
; CALLING SEQUENCE:
;   
; cldy_uplot, grid, NHI, FeH, nH, ions
;
; INPUTS:
;   grid - Cloudy grid
;   NHI  - HI column densities [Can be an array of values]
;   FeH  - Metallicity of the gas [must match grid value]
;   nH   - Volume density of the gas [must match grid value]
;   ions  - Array of [Z,ion] vectors
;   YMNX= - Y values of the plot
;   INFIL= - Ps file to create (instead of plotting to screen)
;
; RETURNS:
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
; getcolor
; x_setclrs
; getabnd
;
; REVISION HISTORY:
;   15-Nov-2001 Written by JXP
;   27-Apr-2006 Modified to look more like published plots, KLC
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_uplot, grid, NHI, FeH, nH, ions, YMNX=ymnx, PSFIL=psfil

;
  if  N_params() LT 5  then begin 
      print, 'Syntax - ' +$
        'cldy_Uplot, grid, NHI, FeH, nH, ions, YMNX=, PSFIL= (v1.1)'
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
  if count LE 1 then begin
      print, 'cldy_uplot: No matches!'
      return
  endif

  xmin = min(grid[ii].U, MAX=xmax)
  xplt = grid[ii].U
  xsrt = sort(xplt)

; Loop on the ions

  sz = size(ions)
  nions = sz[2]
  dumx = fltarr(2)
  dumy = fltarr(2)
  dumx[0] = -99.
  dumy[0] = -99.

; PS FILE
  if keyword_set(PSFIL) then begin
      x_psopen, psfil, /maxs
      !y.margin = [5,2]
  endif

; Set colors
  !P.MULTI = [0, 1, n_NHI]
  xcolors = x_setclrs()
  clr = getcolor(/load)

  nclr = 0                      ;used to make ions of same species same color

  for q = 0, nions-1 do begin
;   Abund

      getabnd, elm, ions[0,q], abnd, flag=1
      yplt = NHI[ww] - grid[ii].X[1,1] - 12. + abnd + $
        grid[ii].X[ions[0,q],ions[1,q]]
      ; Metallicity
      if ions[0,q] NE 1 then yplt = yplt + FeH

      ; Skip using period (psym=2)
      if q EQ 2 then psym = 1 else begin
          psym=q+1
          if q ge 7 then begin 
              psym = (q mod 6)  ;rotate back through psym
              if psym eq 3 then psym = 1
          endif 
      endelse 
      ;psym = psym<7

      if q EQ 0 then plot, xplt[xsrt], yplt[xsrt], $
        psym=-psym, xrange=[xmin,xmax], $
        yrange=[ymnx[0],ymnx[1]], color=xcolors[nclr], $ 
        background=clr.white, ystyle=1, xstyle=1, $
	xtitle='logU', ytitle='logN(X)', charsize=1.5 $
      else begin
          if ions[0,q-1] ne ions[0,q] then nclr = nclr+1 ;change color
          oplot, xplt[xsrt], yplt[xsrt], psym=-psym, color=xcolors[nclr]
      endelse 
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

      ;Print title and legend, upper right-hand corner
      dumx[1] = xmin+0.2        
      dy = ymnx[1]-ymnx[0]      ;use to scale vertical placement
      dx = xmax-xmin 

      dumy[1] = ymnx[1]-0.1*dy-q*0.03*dy 
      oplot, dumx, dumy, psym=psym, color=xcolors[nclr]
      xyouts, xmin+0.07*dx, ymnx[1]-0.11*dy-q*0.03*dy, nm, charsize=1.5, $
        color=xcolors[nclr]

      xyouts, xmin+0.03*dx, ymnx[1]-0.05*dy, $
        'logN(HI)='+string(NHI[ww],format='(f5.2)') + $
        '   [M/H]='+string(FeH[ww],format='(f5.2)'),$
        charsize=1.5, color=clr.black

  endfor
endfor

  if keyword_set(PSFIL) then begin
      x_psclose
      !y.margin = [4,2]
  endif
           
  !P.MULTI = [0, 1, 1]
return
end

