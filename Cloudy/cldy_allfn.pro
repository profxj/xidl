;+ 
; NAME:
; cldy_allfn
;   Version 1.1
;
; PURPOSE:
;    Plots [X/Fe+], [X/H0] vs. Elem for a range of the grid
;
; CALLING SEQUENCE:
;   
; cldy_allfn, grid, obsi, val, ions
;
; INPUTS:
;   grid  -- CLOUDY grid
;   obsi  -- Observed pair of ions
;   val   -- Ratio of observed ions
;   ions  -- Array of [Z,ion] vectors to plot
;   fN_in -- Fraction of the ratio contributed by entirely ionzed gas
;
; RETURNS:
;   
;
; OUTPUTS:
;   Creates a Plot
;
; OPTIONAL OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   NHI   
;   FeH
;   nH
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_allfn, grid, [ [26,3], [26,2] ], -0.5,  [[14,2], [13,2]]
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_allfn, grid, obsi, val, ions, fN_in, $
		YMNX=ymnx, NHI=NHI, FeH=FeH, nH=nH

;
  if  N_params() LT 5  then begin 
      print, 'Syntax - ' +$
        'cldy_allfn, grid, obsi, val, ions, fN_in, NHI=, FeH=, nH= (v1.1)'
      print, 'Assumes obsi has [high,low] order'
      return
  endif 

; Optional keywords

  if not keyword_set( YMNX ) then ymnx = [-0.5, 0.5, -0.5, 1.0]
  if not keyword_set( NHI ) then NHI = [0.0, 50.0]
  if not keyword_set( FeH ) then FeH = [-10., 1.]
  if not keyword_set( nH ) then nH = [-10., 10.]


;  Set up plot
  !P.MULTI = [0, 1, 2]
  !P.CHARSIZE = 1.5
  clr = getcolor(/load)
  xclr = x_setclrs()
  
;  Find Grid Subset

  ii = where(grid.NHI GE NHI[0] AND grid.NHI LE NHI[1] AND $
             grid.FeH GE FeH[0] AND grid.FeH LE FeH[1] AND $
             grid.nH GE nH[0] AND grid.nH LE nH[1]  $
             AND grid.flg NE 0 AND $
             grid.X[obsi[0,0],obsi[1,0]] - grid.X[obsi[0,1],obsi[1,1]] GT val,$
             count)             ; constraint on obsi
  if count LE 1 then return

; Calculate fN()

  fN = alog10( 1. - 10^(val - $
                        (grid[ii].X[obsi[0,0],obsi[1,0]] - $
                         grid[ii].X[obsi[0,1],obsi[1,1]])))

; Constrain on fN()
  jj = where(abs(fN-fN_in) LT 0.05, count)
  if count EQ 0 then begin
      print, 'No fN values work!'
      return
  endif

; Loop on the ions

  sz = size(ions)
  nions = sz[2]
  if nions LT 2 then return
;      Label
  nm = strarr(nions)
  for q=0,nions-1 do begin
      getabnd, elm, ions[0,q], abnd, flag=1
      case ions[1,q] of
          1: nm[q] = elm+'I'
          2: nm[q] = elm+'II'
          3: nm[q] = elm+'III'
          4: nm[q] = elm+'IV'
          5: nm[q] = elm+'V'
          6: nm[q] = elm+'VI'
      endcase
  endfor
          

  for ww=0,1 do begin           ; Plot X/Fe and X/H
      ; Setup the plot
      !x.ticks = nions-1
      !x.tickname=nm
      case ww of 
          0: plot, [0.0, nions-1], [ymnx[0], ymnx[1]], /nodata, $
            color=clr.white, background=clr.black, $
            xtitle='Elements', ytitle='[X/Fe+]', $
            xminor=0, xmargin=[12., 2.], $
            ymargin=[5.,2.]
          1: plot, [0, nions-1], [ymnx[2], ymnx[3]], /nodata, $
            color=clr.white, background=clr.black, $
            xtitle='Elements', ytitle='[X/HI]', $
            xminor=0, xmargin=[12., 2.], $
            ymargin=[5.,2.]
      endcase
      
      
  for q = 0, nions-1 do begin
;   Abund

      XY = alog10( 10^fN[jj] + (1-10^fN[jj])*$
                   10^(grid[ii[jj]].X[ions[0,q],ions[1,q]]-$
                       grid[ii[jj]].X[obsi[0,1],obsi[1,1]])) ; Y
      XHI = alog10( 10^fN[jj] + (1-10^fN[jj])* $
                    10^(grid[ii[jj]].X[ions[0,q],ions[1,q]]-$
                        grid[ii[jj]].X[1,1])) ; HI
                                          
      xplt = fltarr(n_elements(XY))+float(q)

;   XY
      case ww of
          0: oplot, xplt, XY, psym=q+1, color=xclr(q+1)
          1: oplot, xplt, XHI, psym=q+1, color=xclr(q+1)
      endcase
          
;       IONS
;      xyouts, q, ymnx[0]-0.1, nm, ALIGNMENT=0.5, charsize=1.6, color=xclr(q+1)

      endfor
  endfor

           
  !P.MULTI = 0
  !x.ticks = 0
  !x.tickname = ''
return
end

