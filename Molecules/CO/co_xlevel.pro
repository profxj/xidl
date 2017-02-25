;+ 
; NAME:
; co_xlevel
;  (V1.0)
;
; PURPOSE:
;    Calculate the energy levels of the CO molecule ground state (X)
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Sep-2008 Written by JXP with guidance from Y Sheffer
;-
;------------------------------------------------------------------------------
function co_xlevel, vpp, jpp

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'engy = co_xlevel(vpp, jpp) [v1.0]'
    return, -1
  endif 
  ;; All values are from George et al. 1994, J. Mol. Spec. 165, 500

  ;; Masses
  M_C12 = 12.0d
  M_O16 = 15.99491463d
  me = 5.48579903d-4
;  X_00 = 1081.585681888d  ;; Energy above potential minimum (cm^-1)

  ;; U values  :: Units are MHz x AMU^(k/2+l)
  U = dblarr(10,7)   ;; k,l
  U[1,0] =  0.1703231494924D9
  U[2,0] = -0.27313142748D7
  U[3,0] =  0.56137102D4
  U[4,0] =  0.9507155D2
  U[5,0] =  0.1217222D1
  U[6,0] = -0.565173D-1
  U[7,0] = -0.158669D-1
  U[8,0] =  0.601350D-3
  U[9,0] = -0.107034D-4
  
  U[0,1] =  0.397029019829D6
  U[1,1] = -0.9422092076D4
  U[2,1] =  0.8790883D0
  U[3,1] = -0.2603192D-1
  U[4,1] =  0.3374105D-1  
  U[5,1] = -0.2752925D-2  
  U[6,1] =  0.6317853D-4  
  U[7,1] = -0.3686248D-5  

  U[0,2] = -0.86293759048D1  
  U[1,2] =  0.37595887D-2    
  U[2,2] = -0.17588271D-2    
  U[3,2] =  0.4487467D-4     
  U[4,2] = -0.6912644D-5     

  U[0,3] =  0.5687165421D-4  
  U[1,3] = -0.3627044D-5     
  U[2,3] = -0.6903992D-7     

  U[0,4] = -0.23964487D-8
  U[1,4] = -0.126055D-9   
  U[2,4] = -0.22958D-11   

  U[0,5] = -0.208015D-13 
  U[1,5] = -0.704332D-14 

  U[0,6] = -0.473222D-17 

  ;; Delta values  
  DtC = dblarr(10,7)   ;; k,l
  DtC[1,0] =  0.694412285D0
  DtC[2,0] =  0.3159935D0  
  DtC[3,0] = -0.122868D2   

  DtC[0,1] = -0.20563695D1 
  DtC[1,1] = -0.1666043D1  

  DtC[0,2] = -0.739681645D1

  DtO = dblarr(10,7)   ;; k,l
  DtO[1,0] = -0.167780221D0
  DtO[2,0] = -0.109806D1
  DtO[3,0] = -0.3534D1

  DtO[0,1] = -0.20976665D1
  DtO[1,1] = -0.1706614D1

  DtO[0,2] =  0.135D1

  ;; Basis
  vv = (vpp + 0.5)
  jj = jpp*(jpp + 1.)

  ;; Sum it up
  inv_mu = 1./M_C12  + 1./M_O16

  engy = 0.
  for k=0,9 do begin
      for l=0,6 do begin

          ;; Ykl
          Ykl = inv_mu^(k/2. + l) * (1 + me * (DtC[k,l]/M_C12 + $
                                               DtO[k,l]/M_O16) ) * U[k,l]
          ;; Energy (MHz)
          engy = engy + Ykl * vv^k * jj^l
;          print, k, l, Ykl*vv^k*jj^l, engy

      endfor
  endfor

  ;; Convert to cm^-1 (relative to potential minimum not vpp=0, jpp=0)
  engy = engy * 1e6 / 2.9979246d+10 

  ;;
  return, engy
end
