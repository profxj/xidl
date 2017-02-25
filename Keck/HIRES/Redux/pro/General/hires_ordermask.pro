;+ 
; NAME:
; hires_ordermask
;     Version 0.1
;
; PURPOSE:
;    Returns a integer image which identifies each pixel with a given order
;      or inter-order gap.
;
; CALLING SEQUENCE:
;   
;  maskimage = hires_ordermask(ncol, nrow, ordr_str, trim=trim)
;
;    (note: it's a little silly to pass ncol and nrow, but had no better
;           idea at the time.  I should have just passed hires structure,
;           but it doesn't really need it to run.)
;
;
; INPUTS:
;   ncol     - number of columns in image
;   nrow     - number of rows
;   ordr_str - Structure array describing order format
;
; RETURNS:
;
; OUTPUTS:
;  An [NCOL,NROW] image with each pixel exclusively assigned to an order or gap
;
; OPTIONAL KEYWORDS:
;   trim     - The buffer included in each order (default 0. pixels).
;              Increasing trim increases the width of each order!
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Crude routine to make a map of order placement
;     Orders fall at positive Maps  [1...nord]
;     Gaps  are marked with [-1 ... -NGAP]
;     Pixel leftward of 1st order or rightward of last order in 
;           ordr_str are set to 0
;
; EXAMPLES:
;  maskimage = hires_ordermask(ncol, nrow, ordr_str, trim=trim)
;
; PROCEDURES/FUNCTIONS CALLED:
;  gap_index  (internal routine to find those gaps!)
;
; REVISION HISTORY:
;   15-Jul-2003 Checked in by SB
;     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  gap_index jumps through hoops to find the list of image indices which
;       lie between lhedg <= column <= rhedg 
;-
;------------------------------------------------------------------------------

function gap_index, lhedg, rhedg, ncol

     spot = long(total(long(rhedg - lhedg+1),/cumulative))
     nrow = n_elements(spot)
     ;; JXP -- bug fix?
;     if spot[nrow-1] EQ 0 then return, -1L
     if spot[nrow-1] LE 0 then return, -1L

     place = lindgen(spot[nrow-1]) 
     place2 = place*0L
     place2[spot[0:nrow-2]] = 1
     pt2 = long(total(place2,/cumulative))*ncol
      
     ptcol = pt2*0
     ptcol[0:spot[0]-1] = place[0:spot[0]-1] + lhedg[0]
     for i=1, nrow-1 do ptcol[spot[i-1]:spot[i]-1] = $
          place[spot[i-1]:spot[i]-1] - spot[i-1] + lhedg[i]

    
     ptfinal = long(ptcol + pt2)

      ptgood = where(ptfinal GE 0 AND ptfinal LT ncol*nrow $ 
                  AND ptcol GE 0 AND ptcol LT ncol)

      if ptgood[0] EQ -1 then return, -1L $
      else return, ptfinal[ptgood]

end

     
function hires_ordermask, ncol, nrow, ordr_str, trim=trim, noextra=noextra

    if (size(trim,/tname) EQ 'UNDEFINED') then trim=0.0
    if NOT keyword_set(noextra) then noextra=0


    mask = intarr(ncol, nrow)
    fill = replicate(1, ncol)

    nord = n_elements(ordr_str) 
;
;     Include the 0.5 to round properly
;
;

    orders = ordr_str.order

    if noextra EQ 0 then begin

;;;; just to ensure non-negativity
      if max(orders) LT nord+4 then orders=orders+1000L
      maxorder = max(orders)
      minorder = min(orders)
      midorder = (minorder+maxorder)/2.
      orderrange = (maxorder-minorder)

;
;     add a buffer of 3 orders on each end
;
      ordernorm = 2.0d * (1.0d*orders - midorder)/(orderrange + 6.)
      hh = fchebyshev(ordernorm, 5)

      neworders = orders[0] EQ maxorder ? $
          reverse(minorder+findgen(nord+6)-3) : minorder+findgen(nord+6)-3

      orderspot = lindgen(nord) + 3

      hfinal = fchebyshev(2.0d*(neworders - midorder)/ (orderrange + 6), 5)

      alpha = 1.0d * transpose(hh) # hh
      alpha_work = invert(alpha) # transpose(hfinal)


      lhfit = (hh ## ordr_str.lhedg) # alpha_work
      rhfit = (hh ## ordr_str.rhedg) # alpha_work

      lhedg = fix(lhfit - trim + 0.5) 
      rhedg = fix(rhfit + trim + 0.5) 

      ; put original orders back in..., no change to input orders

      lhedg[*,orderspot] =  long(ordr_str.lhedg - trim + 0.5) 
      rhedg[*,orderspot] =  long(ordr_str.rhedg + trim + 0.5) 
      

    endif else begin

      ; just keep orders as reported in ordr_str

      neworders=orders
      lhedg = long(ordr_str.lhedg - trim + 0.5) 
      rhedg = long(ordr_str.rhedg + trim + 0.5) 
    endelse

    final_nord = n_elements(neworders)

    for i=0, final_nord-1 do begin


      lower2 = lhedg[*,i]
      upper2 = rhedg[*,i]
      if i LT final_nord-1 then begin
        upper = (lhedg[*,i+1]-1) > (rhedg[*,i]+1)
        gapo = gap_index(rhedg[*,i]+1,upper,ncol)
        if gapo[0] NE -1 then mask[gapo] = -1.0*neworders[i]

        lower = (rhedg[*,i]) < lhedg[*,i+1] 
        gapl = gap_index(lower, rhedg[*,i],ncol)
        if gapl[0] NE -1 then mask[gapl] = 9999
 
        upper2 = upper2 < (lhedg[*,i+1]-1) 
       
      endif

      if i GT 0 then lower2 = lower2 > (rhedg[*,i-1] + 1)

;
;       which pixels fall within 1/2 pixel of  lh-trim<= x <= rh+trim
      gapi = gap_index(lower2,upper2,ncol)
      if gapi[0] NE -1 then mask[gapi] = neworders[i]

    endfor

    return, mask
end

       
