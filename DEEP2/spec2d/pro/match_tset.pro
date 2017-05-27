;+
; NAME:
;   match_tset
;
; PURPOSE:
;   Match "left" and "right" edge tsets so that they delimit a slit
;
; CALLING SEQUENCE:
;   match_tset, tset1, tset2
;
; INPUTS:
;   tset1      - "left" traceset.  
;   tset2      - "right" traceset. 
;
; OUTPUTS:
;   tset1      - trace set structure for "left" traces (extras discarded)
;          tset1 = $
;            { func    :    'legendre'  , $
;              xmin    :    0  , $
;              xmax    :    4095   , $
;              coeff   :    array[ncoeff, ntrace] $
;            }
;   tset2      - trace set structure for "right" traces (extras discarded)
;
; COMMENTS: 
;   tset1 should contain the traceset of all edges of increasing
;   brightness (moving left to right) and so are called "left" edges.
;   tset2 contains decreasing brightness ("right") edges.  This
;   procedure sorts the positions of the traces (both left and right)
;   and then makes sure that each "left" trace is followed immediately
;   by exactly 1 "right" trace.  Any excess traces are discarded.
;
;   If slitmask information were input, this routine could be expanded
;   to determine WHICH left & right traces go with each slit, and
;   suggest traces that might be missing or redundant. 
;
; BUGS:
;
; EXAMPLES:
;
; REVISION HISTORY:
;   01-Dec-2000  Written by D. Finkbeiner, Berkeley
;-
;------------------------------------------------------------------------------
pro match_tset, tset1, tset2

; take a slice through the middle of the trace set  (doug version)
;  ny = (tset1.xmax-tset1.xmin+1)/2
;  traceset2xy, tset1, ny[0], x1    ; ny[0] so ny doesn't get clobbered
;  traceset2xy, tset2, ny[0], x2

; md version
  dims =  size(tset1.coeff, /dim)
  ny = fltarr(1, dims[1])
  ny[*] = (tset1.xmax-tset1.xmin+1)/2
  traceset2xy, tset1, ny, x1    
  traceset2xy, tset2, ny, x2

  n1 = n_elements(x1)
  n2 = n_elements(x2)

  id = [fltarr(n1)+1, fltarr(n2)+2]
  x = [reform(x1, n1), reform(x2, n2)]
  ind = sort(x)

; where is a 1 followed immediately by a 2?  (Beware wrap-around)
  match = (id[ind] EQ 1) AND (shift(id[ind], -1) EQ 2)
  match[n1+n2-1] = 0B ; take care of wrap
  mind = match OR shift(match, 1)
  m = mind-mind
  m[ind] = mind  ; reorder
  m1 = m[0:n1-1]
  m2 = m[n1:*]

  w1 = where(m1, ct1)
  w2 = where(m2, ct2)

  nbad = n1+n2-ct1-ct2
  print, 'Number rejected by MATCH: ', nbad
  if (ct1 EQ 0) OR (ct2 EQ 0) then message, 'No good traces!!' $
  else begin 
     tset1 = $
       { func    :    tset1.func   , $
         xmin    :    tset1.xmin   , $
         xmax    :    tset1.xmax   , $
         coeff   :    tset1.coeff[*, w1] $
       }
     tset2 = $
       { func    :    tset2.func   , $
         xmin    :    tset2.xmin   , $
         xmax    :    tset2.xmax   , $
         coeff   :    tset2.coeff[*, w2] $
       }
  endelse 

  return
end
