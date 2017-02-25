;+ 
; NAME:
; sdss_filter_nu
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
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
;   2005  Written by Joe Hennawi UCB 
;-
;------------------------------------------------------------------------------
FUNCTION sdss_filter_nu, lognu, filter

COMMON FILTER, lognu_u, F_u, lognu_g, F_g, lognu_r, F_r, lognu_i, F_i $
  , lognu_z, F_z, lognu_B, F_B, lognu_B_J, F_B_J $
  , lognu_J, F_J, lognu_H, F_H, lognu_K, F_K $
  , SPL_U, SPL_G, SPL_R, SPL_I, SPL_Z, SPL_B, SPL_B_J $
  , SPL_J, SPL_H, SPL_K

IF NOT KEYWORD_SET(SPL_U) THEN BEGIN
; evaluate splines for all filters    
    c = 2.9979246d10
    angstrom = 1.0d-8
; files for filter curves
    sfront = getenv('XIDL_DIR')+'/SDSS/Photo/'
    mfront = getenv('XIDL_DIR')+'/SDSS/m912/'
    u_file = sfront+'sdss_u0.res'
    g_file = sfront+'sdss_g0.res'
    r_file = sfront+'sdss_r0.res'
    i_file = sfront+'sdss_i0.res'
    z_file = sfront+'sdss_z0.res'
    B_file = mfront+'john_b.dat'
    BJ_file = mfront+'bj.dat'
    J_file = mfront+'niri_J.dat'
    H_file = mfront+'niri_H.dat'
    K_file = mfront+'niri_K.dat'
; read in sdss m quantum efficiency curves
    readcol, u_file, lam_u, junk1, junk2, F_u13, COMMENT = '\' $
             , FORMAT = 'F,F,F,F,F,F'
    readcol, g_file, lam_g, junk1, junk2, F_g13, COMMENT = '\' $
             , FORMAT = 'F,F,F,F,F,F'
    readcol, r_file, lam_r, junk1, junk2, F_r13, COMMENT = '\' $
             , FORMAT = 'F,F,F,F,F,F'
    readcol, i_file, lam_i, junk1, junk2, F_i13, COMMENT = '\' $
             , FORMAT = 'F,F,F,F,F,F'
    readcol, z_file, lam_z, junk1, junk2, F_z13, COMMENT = '\' $
             , FORMAT = 'F,F,F,F,F,F'
; These are energy efficiencies, so need to be divided by lam to be
; canonical quantum efficiencies. The normalization is such that
; HD19445 has B-V of 0.46
; B filter is at 1 airmass whereas sdss is at 1.3???
    readcol, B_file, lam_B, junk1, F_B1_temp, COMMENT = '\' $
             , FORMAT = 'F,F,F'
    F_B1 = F_B1_temp/lam_B
; #UKSTU B_J - 2mm GG395 with UJ10738P transmission, UKST optics and unit 
; airmass
    readcol, BJ_file, lam_B_J, F_BJ1, comment = '#', format = 'F,F'
; These are warm NIRI IR filters take from the NIRI webpage
; http://www.gemini.edu/sciops/instruments/niri/NIRIFilterList.html
    readcol, J_file, lam_J, F_J1, COMMENT = '#' $
             , FORMAT = 'F'
    readcol, H_file, lam_H, F_H1, COMMENT = '#' $
             , FORMAT = 'F'
    readcol, K_file, lam_K, F_K1, COMMENT = '#' $
             , FORMAT = 'F'
; now compute splines
    nu_u = reverse(c/(lam_u*angstrom))
    nu_g = reverse(c/(lam_g*angstrom))
    nu_r = reverse(c/(lam_r*angstrom))
    nu_i = reverse(c/(lam_i*angstrom))
    nu_z = reverse(c/(lam_z*angstrom))
    nu_B = reverse(c/(lam_B*angstrom))
    nu_B_J = reverse(c/(lam_B_J*angstrom))
;   These were in descending order so no reverse
    nu_J = c/(lam_J*10.0D*angstrom)
    nu_H = c/(lam_H*10.0D*angstrom)
    nu_K = c/(lam_K*10.0D*angstrom)
;   Take log of frequencies
    lognu_u = alog10(nu_u)
    lognu_g = alog10(nu_g)
    lognu_r = alog10(nu_r)
    lognu_i = alog10(nu_i)
    lognu_z = alog10(nu_z)
    lognu_B = alog10(nu_B)
    lognu_B_J = alog10(nu_B_J)
    lognu_J = alog10(nu_J)
    lognu_H = alog10(nu_H)
    lognu_K = alog10(nu_K)
;   filter transmissions in aribtrary units
    F_u = (nu_u/1.0e15)^2*reverse(F_u13)
    F_g = (nu_g/1.0e15)^2*reverse(F_g13)
    F_r = (nu_r/1.0e15)^2*reverse(F_r13)
    F_i = (nu_i/1.0e15)^2*reverse(F_i13)
    F_z = (nu_z/1.0e15)^2*reverse(F_z13)
    F_B = (nu_B/1.0e15)^2*reverse(F_B1)
    F_B_J = (nu_B_J/1.0e15)^2*reverse(F_BJ1)
    F_J = (nu_J/1.0e15)^2*(F_J1)
    F_H = (nu_H/1.0e15)^2*(F_H1)
    F_K = (nu_K/1.0e15)^2*(F_K1)
;   Compute splines    
    SPL_u = SPL_INIT(lognu_u, F_u, /DOUBLE)
    SPL_g = SPL_INIT(lognu_g, F_g, /DOUBLE)
    SPL_r = SPL_INIT(lognu_r, F_r, /DOUBLE)
    SPL_i = SPL_INIT(lognu_i, F_i, /DOUBLE)
    SPL_z = SPL_INIT(lognu_z, F_z, /DOUBLE)
    SPL_B = SPL_INIT(lognu_B, F_B, /DOUBLE)
    SPL_B_J = SPL_INIT(lognu_B_J, F_B_J, /DOUBLE)
    SPL_J = SPL_INIT(lognu_J, F_J, /DOUBLE)
    SPL_H = SPL_INIT(lognu_H, F_H, /DOUBLE)
    SPL_K = SPL_INIT(lognu_K, F_K, /DOUBLE)
ENDIF


CASE FILTER OF 
    'u': BEGIN
        SPL = SPL_u
        lognu_tab = lognu_u
        F_tab = F_u
    END
    'g': BEGIN
        SPL = SPL_g
        lognu_tab = lognu_g
        F_tab = F_g
    END
    'r': BEGIN
        SPL = SPL_r
        lognu_tab = lognu_r
        F_tab = F_r
    END
    'i': BEGIN
        SPL = SPL_i
        lognu_tab = lognu_i
        F_tab = F_i
    END
    'z': BEGIN
        SPL = SPL_z
        lognu_tab = lognu_z
        F_tab = F_z
    END
    'B': BEGIN
        SPL = SPL_B
        lognu_tab = lognu_B
        F_tab = F_B
    END
    'B_J': BEGIN
        SPL = SPL_B_J
        lognu_tab = lognu_B_J
        F_tab = F_B_J
    END
    'J': BEGIN
        SPL = SPL_J
        lognu_tab = lognu_J
        F_tab = F_J
    END
    'H': BEGIN
        SPL = SPL_H
        lognu_tab = lognu_H
        F_tab = F_H
    END
    'K': BEGIN
        SPL = SPL_K
        lognu_tab = lognu_K
        F_tab = F_K
    END
    ELSE: message, 'Unrecognized filter'
ENDCASE

min_log = min(lognu_tab, j_min)
max_log = max(lognu_tab, j_max)

ind_spline = WHERE(lognu GT min_log AND lognu LT max_log)
ind_min = WHERE(lognu LE min_log)
ind_max = WHERE(lognu GE max_log)
    
F_out = dblarr(n_elements(lognu))

IF ind_spline[0] NE -1 THEN $
  F_out[ind_spline] = SPL_INTERP(lognu_tab, F_tab, SPL, lognu[ind_spline])
IF ind_min[0] NE -1 THEN  F_out[ind_min] = 0.0D 
IF ind_max[0] NE -1 THEN  F_out[ind_max] = 0.0D
LN10 = alog(10.0D)
; normalize to unity in integral against lognu.
norm = INT_TABULATED(lognu, LN10*F_out, /DOUBLE)
F_out = F_out/norm

RETURN, F_out
END
