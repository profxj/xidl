function sysinit,mu,gmm=gmm,o3=o3,roll3=roll3,coll_angle=coll_angle,coll_phi=coll_phi,coll_zdst=coll_zdst,t2=t2,p2=p2,p3=p3,cam_angle=cam_angle,cam_phi=cam_phi,cam_foc=cam_foc,x_optaxis=x_optaxis,y_optaxis=y_optaxis,mos_rotation=mos_rotation,norder=norder
;+
; NAME:
;        SYSINIT
;
; PURPOSE:
;        Initialize a structure containing parameters of DEIMOS 
;
; CATEGORY:
;        DEIMOS optical model
;
; CALLING SEQUENCE:
;        sysstructure=sysinit(mu, [gmm= , o3=, roll3= , coll_angle=,
;          coll_phi=, coll_zdst=, t2=, p2= ,p3=, cam_angle=,
;          cam_phi=, cam_foc=, x_optaxis=, y_optaxis=, mos_rotation=, norder=])
; 
; Inputs:
;        mu - Grating tilt (deg)
;
; OPTIONAL INPUTS:
;        Non-default values of parameters are set by optional keywords.
;
; OPTIONAL KEYWORD PARAMETERS:
;        gmm   - Grating grooves/mm (default 900)
;        o3    - Grating 3rd yaw angle (deg; default -0.057)
;        roll3 - Grating roll (deg; default 0.107)
;        coll_angle - Collimator tilt error (deg; default 0.002)
;        coll_phi   - Collimator tilt phi angle  (deg; default 0.0)
;        coll_zdst  - Collimator Z-location (mm; default 2197.1)
;        t2 - 2nd Normal Theta angle (deg; default -0.45)
;        p2 - 2nd Normal Phi angle (deg; default 0.05)
;        p3 - 3rd Normal Phi angle (deg) [currently ignored]
;        cam_angle - Camera angle (nom. 2.33) (deg)
;        cam_phi   - Camera phi (nom. 90) (deg)
;        cam_foc   - Camera focal length (mm; default 382.73)
;        x_optaxis - Camera optical axis X(in ICS coordinates) (pix; def. -11)
;        y_optaxis - Camera optical axis Y(in ICS coordinates) (pix; def. -2)
;        mos_rotation - rotation of CCD mosaic (deg; default -0.148)
;        norder       - grating order (default 1)
;
; OUTPUTS:
;        System structure, as used by QTRMAP and QMODEL
;
; OPTIONAL OUTPUTS:
;        None.
; COMMON BLOCKS:
;        None.
; SIDE EFFECTS:
;        None.
; RESTRICTIONS:
;        None.
;
; EXAMPLE:
;        sys=sysinit(10.,gmm=1200)
;        qtrmap, sys
;        qmodel,sys,xmm,ymm,lambda,xics,yics,ccdnum,xpix,ypix,/cubic
;
;
; MODIFICATION HISTORY:
;        Based on routines by Drew Phillips
;        Finished testing 12/5/01 JAN
;-

if n_elements(coll_angle) eq 0 then coll_angle=0.002    ;0.023
if n_elements(coll_phi) eq 0 then coll_phi=0.
if n_elements(coll_zdst) eq 0 then coll_zdst=2197.1
if n_elements(t2) eq 0 then t2 = -0.5   ;-0.45     ;-0.30100
if n_elements(p2) eq 0 then p2 = 0.0810 ;0.05   ;0.
if n_elements(p3) eq 0 then p3 = 0.0   
if n_elements(o3) eq 0 then o3 = -0.057
if n_elements(roll3) eq 0 then roll3 = 0.107
if n_elements(cam_angle) eq 0 then cam_angle= 2.752 ;2.33
if n_elements(cam_phi) eq 0 then cam_phi=90.d0
;if n_elements(cam_foc) eq 0 then cam_foc=381.
if n_elements(cam_foc) eq 0 then cam_foc=382.00 ;382.73 ;382.75
;if n_elements(x_optaxis) eq 0 then x_optaxis=0.
;if n_elements(y_optaxis) eq 0 then y_optaxis=0.
if n_elements(x_optaxis) eq 0 then x_optaxis=-15.6  ;-11 ;-9.
if n_elements(y_optaxis) eq 0 then y_optaxis=-254.8 ;-2; 9.
;if n_elements(mos_rotation) eq 0 then mos_rotation=0.
if n_elements(mos_rotation) eq 0 then mos_rotation=0.021 ;-0.148   ;-0.221
if n_elements(norder) eq 0 then norder=1
if n_elements(gmm) eq 0 then gmm = 900


t1=coll_angle
p1=coll_phi
t3=mu

sysstruct=create_struct('CAM_FOC',double(cam_foc), $
	'ORDER',double(norder),$
	'GRLINES',1.d-3*gmm, $
	'X_OPT',double(x_optaxis), $
	'Y_OPT',double(y_optaxis), $
	'MOS_ROT',!dtor*double(mos_rotation), $
	'COL_DST',double(coll_zdst), $
	'COL_ERR',!dtor*double(t1), $
	'COL_PHI',!dtor*double(p1), $
	'TNT_ANG',!dtor*double(71.5+t2), $
	'TNT_PHI',!dtor*(90.+p2), $
	'MU',!dtor*double(t3), $
	'GR_YERR',!dtor*double(roll3), $
	'GR_ZERR',!dtor*double(o3), $
	'CAM_ANG',!dtor*double(cam_angle), $
	'CAM_PHI',!dtor*double(cam_phi), $
	'CN_XERR',dblarr(8), $
	'CN_YERR',dblarr(8), $
	'CN_RERR',dblarr(8))

return,sysstruct
end
