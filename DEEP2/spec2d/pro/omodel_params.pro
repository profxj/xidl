pro omodel_params, grname, slider, tltval, roll, o3, mu
;+
; NAME:
;
; OMODEL_PARAMS
;
; PURPOSE:
;
; Given grating, slider, and tilt, determine required parameters to
; apply the optical model
;
; CATEGORY:
;
; optical model
;
; CALLING SEQUENCE:
;
; omodel_params, grname, slider, tltval, roll, o3, mu
;
; INPUTS:
;      grname  -- name of grating (as integer, 0 = mirror)
;      slider  -- slider grating is installed in
;      tltval  -- G3TLTVAL/G4TLTVAL, depending on slider
;
; OUTPUTS:
;      roll    -- roll3 value [Grating roll (deg)]
;      o3      -- o3 value [Grating 3rd yaw angle (deg)]
;      mu      -- mu value (effective grating tilt)
; MODIFICATION HISTORY:
;  jan 2003mar05
;-


; mirror
if slider eq 2 and grname eq 0 then begin
    roll=0.
    o3=0.
    mu=-19.423
endif


; values taken from: http://www.ucolick.org/~phillips/deimos_ref/

if slider eq 3 then begin

  case grname of
    600:   begin
             roll=0.145d0
             o3= -0.008d0
             mu = tltval*(1-5.6d-4) -0.182
           end
    831:   begin ; no values - guessing
             roll=0.143d0
             o3= 0.00d0
             mu = tltval*(1-5.6d-4) -0.182
           end
    900:   begin 
             roll=0.141d0
             o3= 0.008d0
             mu = tltval*(1-5.6d-4) -0.134
           end
    1200:  begin 
             roll=0.145d0
             o3= 0.055d0
             mu = tltval*(1-5.6d-4) -0.181
           end
    else:  begin 
             roll=0.145d0
             o3= 0.0d0
             mu = tltval*(1-5.6d-4) -0.182
           end             
    endcase
endif 

if slider eq 4 then begin

  case grname of
    600:   begin ; no values - guesses [ from Drew, 3/10/03]
             roll=-0.065d0
             o3= 0.063d0
             mu = tltval*(1-6.9d-4) -0.298
           end
    831:   begin 
             roll=-0.034d0
             o3= 0.060d0
             mu = tltval*(1-6.9d-4) -0.196
           end
    900:   begin 
             roll=-0.064d0
             o3= 0.083d0
             mu = tltval*(1-6.9d-4) -0.277
           end
    1200:  begin 
             roll=-0.052d0
             o3= 0.122d0
             mu = tltval*(1-6.9d-4) -0.294
           end
    else:  begin 
             roll=-0.05d0
             o3= 0.08d0
             mu = tltval*(1-6.9d-4) -0.25
           end             
    endcase

endif

return

end
