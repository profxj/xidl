pro deimos_red_mosaic,  ampmode
;+
; NAME:
;    deimos_red_mosaic
;
; PURPOSE:
;    generates polynomial solution to nonlinear resonse of DEIMOS CCD mosaic
;
; CALLING SEQUENCE:
;    deimos_red_mosaic, ampmode
;
; INPUTS:
;     ampmode -- FITS keyword in DEIMOS header
;
; OUTPUTS:
;    defines system variables used by DEIMOS_ADU2E 
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;  this should be called once at startup
; COMMENTS:
;  nonlinear behavior seems pretty weird: best fit with cubic spline.
; 
; REVISION HISTORY:
;   13apr02  md
;----------------------------------------------------------------------

; Amplifier     Vdd voltage gain(e/dn)   ange of deviation
; ---------     -----------     ------------------
; 1A            19.0          1.23   +0.15 to -0.41 %  
; 1B            19.0          1.24   +0.20 to -0.31 %
; 2A            19.5          1.23   +0.48 to -0.31 %  
; 2B            18.5          1.20   +0.41 to -0.21 %  
; 3A            18.5          1.20   +0.19 to -0.12 %  
; 3B            18.5          1.27   +0.37 to -0.38 %  
; 4A            19.0          1.25   +0.20 to -0.40 %
; 4B            19.0          1.27   +0.14 to -0.42 %  
; 5A            20.0          1.26   +0.15 to -0.25 % 
; 5B            19.5          1.26   +0.40 to -0.29 % 
; 6A            19.5          1.22   +0.32 to -0.18 % 
; 6B            18.5          1.25   +0.35 to -0.22 % 
; 7A            18.5          1.40   +0.10 to -0.11 % 
; 7B            19.5          1.25   +0.24 to -0.47 %
; 8A            18.5          1.25   +0.43 to -0.27 %
; 8B            19.3          1.25   +0.37 to -0.30 %


time =  [0., 10., 20.,  40.,  80., 160., 320.]
adu    = [ [0., 1744.7, 3488.0, 6976.1, 13958.2, 27933.9, 55555.0], $ ; 1A
           [0., 1706.5, 3410.7, 6821.9, 13657.3, 27344.4, 54410.7], $ ; 1B
           [0., 1611.8, 3229.6, 6468.8, 13040.7, 26356.5, 53082.0], $ ; 2A
           [0., 1669.9, 3339.5, 6687.5, 13407.2, 26884.0, 53510.7], $ ; 2B
           [0., 1344.7, 2689.5, 5378.5, 10768.3, 21579.5, 43162.6], $ ; 3A
           [0., 1671.1, 3340.5, 6684.6, 13392.7, 26845.8, 53287.0], $ ; 3B
           [0., 1707.3, 3414.9, 6833.4, 13674.0, 27371.1, 54417.0], $ ; 4A
           [0., 1751.0, 3502.3, 7005.1, 14018.1, 28041.3, 55775.1], $ ; 4B
           [0., 1689.5, 3379.8, 6760.7, 13526.7, 27075.0, 53931.9], $ ; 5A
           [0., 1812.4, 3628.6, 7265.4, 14563.5, 29200.1, 58166.4], $ ; 5B
           [0., 1727.4, 3457.3, 6917.6, 13851.0, 27777.0, 55381.1], $ ; 6A
           [0., 1713.7, 3429.2, 6868.6, 13758.4, 27573.9, 54902.2], $ ; 6B
           [0., 1587.1, 3176.8, 6355.8, 12719.2, 25442.9, 50778.1], $ ; 7A
           [0., 1678.7, 3356.8, 6712.9, 13434.5, 26908.4, 53436.5], $ ; 7B
           [0., 1668.2, 3336.4, 6678.1, 13389.8, 26878.3, 53642.0], $ ; 8A
           [0., 1604.4, 3214.2, 6437.0, 12897.3, 25841.9, 51431.5] ]

invgain= [1.23, 1.24, 1.23, 1.20, 1.20, 1.27, 1.25, 1.27, $
          1.26, 1.26, 1.22, 1.25, 1.40, 1.25, 1.25, 1.25 ]

  single = strmid(ampmode,0,6) EQ 'SINGLE' ; is it single-amp mode?

; mapping of DEIMOS ccd amplifier to FITS HDU
hdu_amp = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] 

; if A side
if (ampmode eq 'SINGLE:A') then $
       hdu_amp = [2,4,6,8,9,11,13,15]

; if B side
if (ampmode eq 'SINGLE:B') then $
       hdu_amp = [1,3,5,7,10,12,14,16]

esignal = adu*0.
sset=adu*0.
electrons=esignal

electron_fit = fltarr(2,16)

for i=0,15 do begin
   esignal[*,i] = adu[*,i] * invgain[i]
   electron_fit[*,i] = poly_fit(time[0:4],esignal[0:4,i], 1, measure_errors= $
          sqrt(esignal[0:4,i])+.0001,  yfit=outfit )
; slope is fit only over first 4 data points to be in very linear domain.

;extrapolate electrons to larger integration times
   electrons[*,i] = electron_fit[0,i] + electron_fit[1,i]*time
;input electron signal
;   if i eq 0 then plot, electrons[*,i], esignal[*,i]-electrons[*,i], yr=[-400,100] $
;     else oplot, electrons[*,i], esignal[*,i]-electrons[*,i]
 
;get spline fit of nonlinear behavior
   sset[*,i] = spl_init(adu[*,i],electrons[*,i]) 


endfor 

;copy to newly defined system variables
read_only = 1

if single then begin
	defsysv, '!deimos_adu1',adu,read_only
	defsysv, '!deimos_electrons1', electrons,read_only
	defsysv, '!deimos_gsset1',sset,read_only
	defsysv, '!deimos_hduamp1',hdu_amp,read_only
endif else begin
       defsysv, '!deimos_adu2',adu,read_only
        defsysv, '!deimos_electrons2', electrons,read_only
        defsysv, '!deimos_gsset2',sset,read_only
        defsysv, '!deimos_hduamp2',hdu_amp,read_only
endelse


return
end
