;+ 
; NAME:
; hamspec_headid   
;     Version 1.1
;
; PURPOSE:
;    Guess the type of HIRES image based on exposure time, counts,
;    and other characteristics of the image.
;
; CALLING SEQUENCE:
;   
;  guess = hamspec_headid(head, keywd)
;
; INPUTS:
;    head -- Header
;    keywd -- Structure contating the tags describing the keywords
;
; RETURNS:
;  Image type::  'UNK', 'BAD', 'ZRO', 'ARC', 'TFLT', 'PFLT', 'STD',
;  'OBJ', 'MFLT', 'TWI'
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   guess = hamspec_headid(data, xbin, ybin, filename='file')
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------
function hamspec_headid, head, keywd, SILENT=silent

  ;;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'guess = hamspec_headid(header, keywd, /SILENT) [v1.1]'
      return, -1
  endif 

  ;; Exposure
  exp = sxpar(head,keywd.expcrd)
  deck = strtrim(sxpar(head,keywd.decker),2)

  ;; Hatch
;  hat = sxpar(head,keywd.hatch)
;  xdangl = sxpar(head,keywd.xdangl)
  

  ;; Calibs
  guess = ''
  case strcompress(sxpar(head,keywd.lamp),/remove_all) of
     'Thorium-Argon': guess = 'ARC'
     'PolarQuartz': guess = 'TFLT'
     'Polar-Quartz': guess = 'TFLT'
     'Off': 
     '0': 
     else: stop
  endcase
  if strlen(guess) GT 0 then return, guess


  ;; Dark or Bias
  if exp LE 0. then guess = 'ZRO' else $
     if strcompress(sxpar(head,keywd.obstyp)) EQ 'DARK' then guess = 'DRK'
  if strlen(guess) GT 0 then return, guess

  ;; Object or Standard
  if exp LE 60. then guess = 'STD' else guess = 'OBJ'

  return, guess

end
    
