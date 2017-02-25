;+ 
; NAME:
; cldy_2phas   
;   Version 1.0
;
; PURPOSE:
;    Finds the best 2-phase solution for a series of ratio
;    constraints.  This is not a well devloped or tested routine.
;
; CALLING SEQUENCE:
;   
; cldy_2phas, grid, ratio, sig, ion1, ion2, restrct, NHILMT=, FEHLMT=, /UONLY
;
; INPUTS:
;   grid  - CLOUDY grid
;   ratio -  Observed ionic ratios
;   sig  -  Error on the ratios
;   ion1 - Z, i for ion1  (can be an array of Z,i pairs)
;   ion2 - Z, i for ion2  (can be an array of Z,i pairs)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   NHILMT - 2-element array giving NHI min,max
;   FEHLMT - 2-element array giving FeH min,max
;   UONLY  - Value of nH to use
;
; OPTIONAL OUTPUTS:
;   retstrct - Returns a structure of the output
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_2phas, grid, [0.3], [0.1], [14,4], [14,2], soltn
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_2phas, grid, ratio, sig, ion1, ion2, retstrct, $
                NHILMT=nhilmt, FEHLMT=fehlmt, UONLY=uonly

;
  if  N_params() LT 5  then begin 
      print, 'Syntax - ' +$
        'cldy_2phas, grid, ratio, sig, ion1, ion2, retstrct, ' 
      print, ' NHILMT=, FEHLMT=, UONLY= (v1.0)'
      return
  endif 

  flg_subset = 0

; Optional keywords

  if keyword_set( NHILMT ) then flg_subset = flg_subset + 1
  if keyword_set( FEHLMT ) then flg_subset = flg_subset + 2
  if keyword_set( UONLY ) then flg_subset = flg_subset + 4

  npnt = n_elements(ratio)
  Aval = dblarr(npnt)
  sigA = dblarr(npnt)
  val = dblarr(npnt)


; Weight is inverse variance 

   chk = where(sig EQ 0.0, cnt)
   if cnt GT 0 then begin
       print, 'Cant have zero error!'
       return
   endif


;  Find Grid Subset

   case flg_subset of
       0: sub = where(grid.flg NE 0)
       1: sub = where(grid.flg NE 0 AND grid.NHI GE nhilmt[0] AND $
                     grid.NHI LT nhilmt[1])
       2: sub = where(grid.flg NE 0 AND grid.FeH GE fehlmt[0] AND $
                     grid.FeH LT fehlmt[1])
       3: sub = where(grid.flg NE 0 AND $
                      grid.FeH GE fehlmt[0] AND grid.FeH LT fehlmt[1] AND $
                      grid.NHI GE NHIlmt[0] AND grid.NHI LT NHIlmt[1])
       4: sub = where(grid.flg NE 0 AND grid.nH EQ double(Uonly))
       5: sub = where(grid.flg NE 0 AND grid.nH EQ double(Uonly) AND $
                      grid.NHI GE NHIlmt[0] AND grid.NHI LT NHIlmt[1])
       6: sub = where(grid.flg NE 0 AND grid.nH EQ double(Uonly) AND $
                      grid.FeH GE fehlmt[0] AND grid.FeH LT fehlmt[1] )
       7: sub = where(grid.flg NE 0 AND grid.nH EQ double(Uonly) AND $
                      grid.FeH GE fehlmt[0] AND grid.FeH LT fehlmt[1] AND $
                      grid.NHI GE NHIlmt[0] AND grid.NHI LT NHIlmt[1])
       else: return
   endcase
           
   nsub = n_elements(sub)
   chisq = fltarr(nsub,nsub) - 1.
   ans = fltarr(nsub,nsub) - 1.
   sigans = fltarr(nsub,nsub) - 1.

;  Calculate the Ratios
   lgrid = fltarr(nsub, npnt)
   for i=0,nsub-1 do begin
       for k=0,npnt-1 do $
         lgrid[i,k] = grid[sub[i]].X[ion1[0,k],ion1[1,k]]-$
                   grid[sub[i]].X[ion2[0,k],ion2[1,k]]
   endfor
                   

;  Calculate weighted zeropoint

   min = [-1L,-1L]
   minchi = 9 * 10^9

   for i=0,nsub-1 do begin
       for j=0,nsub-1 do begin
;  Don't calculate for i=j
           if(i EQ j) then continue


           for k=0,npnt-1 do begin
               denom = lgrid[i,k] - lgrid[j,k]
               if(denom EQ 0.) then continue
;          Calculate A, sigA
               Aval[k] = (ratio[k]-lgrid[j,k])/denom
               sigA[k] = sig[k]/abs(denom)
           endfor

           weight = 1./sigA^2
           ans[i,j] = total(Aval*weight)/total(weight)
;   Require 0 > ans > 1
           if (ans[i,j] LT 0 OR ans[i,j] GT 1) then continue
           sigans[i,j] = sqrt(1./total(weight))
;  Chisq
           dumchi = 0.d
           for k=0,npnt-1 do begin
               dumchi = dumchi + ((ratio[k] - ans[i,j]*lgrid[i,k] - $
                                   (1.-ans[i,j])*lgrid[j,k])/sig[k])^2
           endfor
           chisq[i,j] = dumchi

;      Minimum
           if chisq[i,j] LT minchi then begin
               min = [i,j]
               minchi = chisq[i,j]
               for k=0,npnt-1 do $
                 val[k] = ans[i,j]*lgrid[i,k] + (1.-ans[i,j])*lgrid[j,k]
           endif
           
       endfor
   endfor

; Make the Structure

   retstrct = { sub: sub, $
                ans: ans, $
                sigans: sigans, $
                chisq: chisq, $
                min: min, $
                val: val }
           

   delvarx, lgrid

   return
end

