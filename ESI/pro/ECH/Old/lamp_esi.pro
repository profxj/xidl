;+
; NAME:
;   lamp.pro
;
; PURPOSE:
;   Create a structure containing information about tabulated arc lines 
;        lamp = 
;           { lambda    :  0.0d0
;             loglam    :  0.0d0
;             intensity :  0.0d0
;             quality   :  0.0d0            
;            }          
;
; CALLING SEQUENCE:
;      rsh=lamp(name,quality='quality')
;
; INPUTS:
;   name    -string containing the name of the lamp file 
;            Supposed to be in  /deep1/marinoni/deep/spec2d/pro/MIA
;
; OPTIONAL INPUTS:
;   
;
; REQUIRED KEYWORDS:
;   
;
; OPTIONAL KEYWORDS: 
;   quality  -can be set to string 'GOOD'  
;
; OUTPUTS:
;   
;
; OPTIONAL OUTPUTS:
;   
;   
; COMMENTS:
;
;
; EXAMPLES:
;
;
; BUGS:
;
;
; PROCEDURES CALLED:
;   
;
; REVISION HISTORY: written  by df
;                   modified by cm 08-01-2001
;   
;-
;--------------------------------------------------------------------------- 

function lamp, name, quality=quality

  if NOT keyword_set(quality) then quality = 'GOOD'
;  lampfilename = '/deep1/marinoni/deep/spec2d/pro/MIA/'+name 
  lampfilename = name 
  readcol, lampfilename, lampwave, lampinten, lampquality, format='D,F,A', /sil
  lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0.0d0}
  lamps = replicate(lamps, N_elements(lampwave))
  lamps.lambda = lampwave
  lamps.loglam = alog10(lampwave)
  lamps.intensity = lampinten
  lamps.good = strupcase(lampquality) EQ 'GOOD' AND lampinten GT 0
  
  if quality eq 'GOOD' then begin
     itrim = where(lamps.good, ct)  
     lamp = lamps(itrim)
  endif 
  return, lamp
end
