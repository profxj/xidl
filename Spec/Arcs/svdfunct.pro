; $Id: //depot/idl/IDL_71/idldir/lib/svdfunct.pro#1 $
;
; Distributed by ITT Visual Information Solutions.
;
;       Default function for SVDFIT
;
;       Accepts scalar X and M, returns
;       the basis functions for a polynomial series.
;
; History:
;    08 Jun 2011  faster version, called by x_findnpeaks & co., KLC
;
function svdfunct,X,M
  
  compile_opt idl2
  ;; ;;;;;;;;;
  ;; Original code
  ;; x_findnpeaks calls zfitmax (IDLSPEC2D_DIR) which calls
  ;; svdfit, which calls this function. It's MUCH faster to
  ;; call it without size() and without looping.
  ;; Just need to make sure $XIDL_DIR is before $IDL_DIR/lib in
  ;; $IDL_PATH 
;        XX=X[0]                 ; ensure scalar XX
;	sz=reverse(size(XX))    ; use size to get the type
;        IF sz[n_elements(sz)-2] EQ 5 THEN $
;		basis=DBLARR(M) else basis=FLTARR(M)
;;
;;       Calculate and return the basis functions
;;
;        basis[0]=1.0
;        FOR i=1,M-1 DO basis[i]=basis[i-1]*XX
  ;; ;;;;;;;;;

        basis=x[0]^indgen(m)  ; IDL will convert to type of x[0]
	return,basis
end
