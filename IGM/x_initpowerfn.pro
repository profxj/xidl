;+ 
; NAME: 
; x_initpowerfn   
;    Version 1.1
;
; PURPOSE:
;    Initialize a power-law f(N,X)
;
; CALLING SEQUENCE:
;   
;  x_initpowerfn, powerfn
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;  powerfn -- Structure containing a power-law f(N,X)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES/FUNCTIONS CALLED:
;  
;
; REVISION HISTORY:
;   May-2011 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Initalize the powerfn structure using empirical measures
pro x_initpowerfn, powerfn_strct

  ;; Create
  powerfn_strct = replicate({powerfnstrct}, 10)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;
  ;; z ~ 2
  idx = 2
  powerfn_strct[idx].zmnx = [0.5, 3.0]
  powerfn_strct[idx].npivot = 7  ;; # of pivots without padding
  powerfn_strct[idx].pivots[0:8] = [0, 12., 14.5, 16.9, 17.5, 19., 20.3, 21.5, 99]

  ;; Cutoff in Lya forest 
  powerfn_strct[idx].fn_pivot[0] = -99.
  powerfn_strct[idx].beta[0] = 0.

  ;; Low z Lya forest (from Ribaudo for now)
  powerfn_strct[idx].fn_pivot[1] = -9.58
  powerfn_strct[idx].beta[1] = -1.50

  ;; LLS
  powerfn_strct[idx].fn_pivot[3] = -17.66
  powerfn_strct[idx].beta[3] = -1.8

  ;; Interpolate across between Lya forest and LLS
  powerfn_strct[idx].fn_pivot[2] = powerfn_strct[idx].fn_pivot[1] + $
                                   (powerfn_strct[idx].pivots[2]-powerfn_strct[idx].pivots[1])*$
                                   powerfn_strct[idx].beta[1]
  powerfn_strct[idx].beta[2] = (powerfn_strct[idx].fn_pivot[3]- powerfn_strct[idx].fn_pivot[2]) / $
                               (powerfn_strct[idx].pivots[3]- powerfn_strct[idx].pivots[2])
                               
  ;; LLS
  powerfn_strct[idx].fn_pivot[4] = -18.9
  powerfn_strct[idx].beta[4] = -0.9  ;; Set this from l(z) of tau>2

  ;; SLLS
  powerfn_strct[idx].fn_pivot[5] = powerfn_strct[idx].fn_pivot[4] + $
                                   (powerfn_strct[idx].pivots[5] - $
                                    powerfn_strct[idx].pivots[4])* $
                                   powerfn_strct[idx].beta[4] 
  powerfn_strct[idx].beta[5] = -1.2  ;; Set this from DLA

  ;; DLA (set this from PW09)
  powerfn_strct[idx].fn_pivot[6] = powerfn_strct[idx].fn_pivot[5] + $
                                   (powerfn_strct[idx].pivots[6] - $
                                    powerfn_strct[idx].pivots[5])* $
                                   powerfn_strct[idx].beta[5] 
  powerfn_strct[idx].beta[6] = -1.7  ;; Set this from DLA

  ;; DLA2 (set this from PW09)
  powerfn_strct[idx].fn_pivot[7] = powerfn_strct[idx].fn_pivot[6] + $
                                   (powerfn_strct[idx].pivots[7] - $
                                    powerfn_strct[idx].pivots[6])* $
                                   powerfn_strct[idx].beta[6] 
  powerfn_strct[idx].beta[7] = -6.5  ;; Set this from DLA
  

  ;; Redshift evolution of XXX (takes into account dX/dz)
  powerfn_strct[idx].zpivot = 2.
  powerfn_strct[idx].gamma = 0.5   ;; Shallow evolution in dX
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  return
end
