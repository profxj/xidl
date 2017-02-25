;+ 
; NAME: 
; x_initm95_fnz
;    Version 1.1
;
; PURPOSE:
; Initalize the powerfn structure using empirical measures
; This is for f(N,z) as modeled in Madau 1995
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
pro x_initm95_fnz, powerfn_strct

  ;; Create
  powerfn_strct = replicate({powerfnstrct}, 1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;
  ;; z ~ 2
  idx = 0
  powerfn_strct[idx].zmnx = [0., 5.]
  powerfn_strct[idx].npivot = 3  ;; # of pivots without padding
  powerfn_strct[idx].pivots[0:4] = [0, 11.3, 17.2, 20., 99]

  ;; Cutoff in Lya forest 
  powerfn_strct[idx].fn_pivot[0] = -99.
  powerfn_strct[idx].beta[0] = 0.

  ;; Lya forest 
  powerfn_strct[idx].fn_pivot[1] = -9.57 ;; For 11.3
;  powerfn_strct[idx].fn_pivot[1] = -11.07 ;; For 12.3
  powerfn_strct[idx].beta[1] = -1.5

  ;; LLS + DLA
  powerfn_strct[idx].fn_pivot[2] = -17.52
  powerfn_strct[idx].beta[2] = -1.5

  ;; LLS + DLA
  powerfn_strct[idx].fn_pivot[3] = -99.
  powerfn_strct[idx].beta[3] = -4.5

  ;; Redshift evolution  (does not include dX/dz)
  powerfn_strct[idx].zpivot[1:2] = 0.
  powerfn_strct[idx].gamma[1:2] = [2.46, 0.68]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  return
end
