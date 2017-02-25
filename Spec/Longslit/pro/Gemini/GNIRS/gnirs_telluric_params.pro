; NAME:
;   gnirs_telluric_params
;
; PURPOSE:
;   Compute physical parameters for a given stellar type
;
; CALLING SEQUENCE:
;   gnirs_telluric_params, type, [ logR = , logT = , logg = , M_V = M_V]
;
; INPUTS:
;   type         - Stellar type. A Dwarf (V) is assumed. 
;
; OPTIONAL OUTPUTS:
;   logR         - Log of R/R_sol
;   logT         - Log of T_eff in K
;   logg         - Log of g in cm/s^2
;   M_V          - V-band absolute magnitude
;
; COMMENTS:
;   This routine uses the table provided in Schmidt-Kaler (1982) pg 456 
;   to determine T_eff and logg. For stellar types not listed in the table
;   i.e. G0, the routine will interpolate between the two nearest types. 
;  
; EXAMPLES:
; gnirs_telluric_params, 'G0V', logR = logR, logT = logT, logg = logg, M_V = MV
; 
; BUGS:
;    
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   29-May-2006 Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------
PRO GNIRS_TELLURIC_PARAMS, type, logR = logR, T = T, logg = logg $
                           , M_V = M_V

logg_sol = alog10(6.67259d-8) + alog10(1.989d33) - 2.0D*alog10(6.96d10)

type_let = strmid(type, 0, 1)
type_num =  long(strmid(type, 1, 2))
file = getenv('XIDL_DIR') + $
       '/Spec/Longslit/calib/standards/kurucz93/schmidt-kaler_table.txt'
readcol, file, stype1, logT1, T1, BminV1, M_V1, BC1, M_bol1, L1 $
         , format = 'A,F,F,F,F,F,F,F'
stype_let = strmid(stype1, 0, 1)
stype_num = long(strmid(stype1, 1, 2))

ind = WHERE(stype_let EQ type_let AND stype_num EQ type_num, nmatch)
IF nmatch EQ 1 THEN BEGIN
    M_bol = M_bol1[ind]
    M_V = M_V1[ind]
    logT = logT1[ind]
    T = T1[ind]
    L = L1[ind]
ENDIF ELSE IF nmatch EQ 0 THEN BEGIN
    typeind = WHERE(stype_let EQ type_let, ntype)
    IF ntype EQ 0 THEN message, 'Problem with telluric params'
    M_bol_vec = M_bol1[typeind]
    M_V_vec   = M_V1[typeind]
    logT_vec = logT1[typeind]
    T_vec    = T1[typeind]
    L_vec    = L1[typeind]
    num_vec = stype_num[typeind]
    M_bol = interpol(M_bol_vec, num_vec, type_num)
    M_V   = interpol(M_V_vec, num_vec, type_num)
    logT  = interpol(logT_vec, num_vec, type_num)
    T = interpol(T_vec, num_vec, type_num)
    L = interpol(L_vec, num_vec, type_num)
ENDIF

; relation between radius, temp, and bolometric luminosity
logR = 0.2*(42.26D - M_bol - 10.0D*logT)
; mass-bolometric luminosity relation 
; from schimdt-kaler p28 valid for M_bol < 7.5
logM = 0.46D - 0.10*M_bol
logg = logM-2.0*logR + logg_sol

RETURN
END
