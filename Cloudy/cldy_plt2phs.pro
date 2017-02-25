;+ 
; NAME:
; cldy_plt2phs
;   Version 1.1
;
; PURPOSE:
;  Plots a 2 phase solution (one purely neutral) given a Cloudy
;  grid and related values.  The user also inputs a pair of ions and
;  a ratio between them.
;    
;
; CALLING SEQUENCE:
;   
; cldy_plt2phs, grid, NHI, FeH, nH, obsi, val, ions
;
; INPUTS:
;   grid  - CLOUDY grid
;   NHI   - N(HI) value
;   FeH   - Metallicity of the gas
;   nH    - Hydrogen volume density
;   obsi  - Observed pair of ions
;   val   - Ratio of observed ions
;   ions  - Array of [Z,ion] vectors to plot
;   [YMNX=] - Y limits of the plot
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
;  Developed for the analysis in Prochaska et al. 2002 (Q1755)
;
; EXAMPLES:
;   cldy_plt2phs, grid, 19.0d, -1.0d, -1.0d, [ [26,3], [26,2] ], -0.5,
;   [[14,2], [14,3]]
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_plt2phs, grid, NHI, FeH, nH, obsi, val, ions, YMNX=ymnx

;
  if  N_params() LT 5  then begin 
      print, 'Syntax - ' +$
        'cldy_plt2phs, grid, NHI, FeH, nH, obsi, val, ions, YMNX= (v1.1)'
      print, 'Assumes obsi has [high,low] order'
      return
  endif 

; Optional keywords

  if not keyword_set( YMNX ) then ymnx = [-0.5, 0.5]


;  Loop on Plot
  !P.MULTI = [0, 1, 2]
  !P.CHARSIZE = 1.5
  clr = getcolor(/load)
  xclr = x_setclrs(/white)
  
;  Find Grid Subset

  ii = where(grid.NHI EQ NHI AND grid.FeH EQ FeH AND grid.nH EQ nH $
             AND grid.flg EQ 1 $
             AND grid.X[obsi[0,0],obsi[1,0]] - grid.X[obsi[0,1],obsi[1,1]] GT val,$
             count)             ; constraint on obsi
  if count LE 1 then return

; Calculate fN()

  fN = alog10( 1. - 10^(val - $
                        (grid[ii].X[obsi[0,0],obsi[1,0]] - $
                         grid[ii].X[obsi[0,1],obsi[1,1]])))
  xmin = min(fN) - 0.1
  xmax = 0.03

; Loop on the ions

  sz = size(ions)
  nions = sz[2]
  if nions LT 2 then return

  for ww=0,1 do begin
  for q = 0, nions-1 do begin
;   Abund

      getabnd, elm, ions[0,q], abnd, flag=1
      XY = alog10( 10^fN + (1-10^fN)*10^(grid[ii].X[ions[0,q],ions[1,q]]-$
                                          grid[ii].X[obsi[0,1],obsi[1,1]])) ; Y
      XHI = alog10( 10^fN + (1-10^fN)*10^(grid[ii].X[ions[0,q],ions[1,q]]-$
                                          grid[ii].X[1,1])) ; HI
                                          

      if q EQ 2 then psym = 1 else psym = q+1
;   XY
      case ww of
          0: begin
              if q EQ 0 then plot, fN, XY, psym=psym, xrange=[xmin,xmax], $
                yrange=[ymnx[0],ymnx[1]], color=clr.white, $
                background=clr.black, $
                xmargin=[12., 2.], ymargin=[5.,2.], $
                xtitle='fN', ytitle='[X/Fe+]' $
              else oplot, fN, XY, psym=psym, color=xclr[q+1]
          end
          1: begin
              if q EQ 0 then begin
                  plot, fN, XHI, psym=psym, xrange=[xmin,xmax], $
                    yrange=[ymnx[0],ymnx[1]], color=clr.white, $
                    xmargin=[12., 2.], ymargin=[5.,2.], $
                    background=clr.black, $
                    xtitle='fN', ytitle='[X/HI]' 
;                  XHI = alog10( 10^fN + $
;                                (1-10^fN)*10^(grid[ii].X[obsi[0,1],obsi[1,1]]-$
;                                              grid[ii].X[1,1])) ; HI
;                  oplot, fN, XHI, psym=psym, color=nions*30
              end else oplot, fN, XHI, psym=psym, color=xclr[q+1]
          end
      endcase

;      Label
      case ions[1,q] of
          1: nm = elm+'I'
          2: nm = elm+'II'
          3: nm = elm+'III'
          4: nm = elm+'IV'
          5: nm = elm+'V'
          6: nm = elm+'VI'
      endcase

;       IONS
      oplot, [xmin+0.02], [ymnx[1]-(q+1)*0.05], psym=psym, color=xclr[q+1]
      xyouts, xmin+0.03, ymnx[1]-(q+1)*0.05-0.03, nm, charsize=1.3, $
        color=xclr[q+1]
;       NHI      
      xyouts, xmin+0.05, ymnx[0]+0.05, 'N(HI)='+strmid(strtrim(NHI,2),0,5), $
        charsize=1.5, color=clr.white
     
      if q EQ 0 then oplot, [xmin, xmax], [0., 0.], color=clr.white

      

  endfor
endfor

           
  !P.MULTI = [0, 1, 1]
return
end

