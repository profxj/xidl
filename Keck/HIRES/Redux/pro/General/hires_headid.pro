;+ 
; NAME:
; hires_headid   
;     Version 1.1
;
; PURPOSE:
;    Guess the type of HIRES image based on exposure time, counts,
;    and other characteristics of the image.
;
; CALLING SEQUENCE:
;   
;  guess = hires_headid(head, keywd)
;
; INPUTS:
;    head -- Header
;    keywd -- Structure contating the tags describing the keywords
;
; RETURNS:
;  Image type::  'UNK', 'BAD', 'ZRO', 'ARC', 'TFLT', 'TFLT', 'STD',
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
;   guess = hires_headid(data, xbin, ybin, filename='file')
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------
function hires_headid, head, keywd, SILENT=silent

  ;;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'guess = hires_headid(header, keywd, /SILENT) [v1.1]'
      return, -1
  endif 

  ;; Exposure
  exp = sxpar(head,keywd.expcrd)
  deck = strtrim(sxpar(head,keywd.decker),2)

  ;; Hatch
  hat = sxpar(head,keywd.hatch)
  xdangl = sxpar(head,keywd.xdangl)

  if hat then begin ;; On Sky or Dome, assuming Sky for now
      ;; 
      if exp LE 5. then guess = 'STD' else guess = 'OBJ'
      case strtrim(sxpar(head,keywd.lamp)) of
          'none': 
          'ThAr2': guess = 'ARC'
          'ThAr1': guess = 'ARC'
          'quartz': flg = 2
          'quartz1': flg = 2
          'quartz2': flg = 2
          'undefined': 
          else: stop
      endcase
      if keyword_set(flg) then begin
          if deck EQ 'D5' then guess = 'TRC' else guess='TFLT'
          if xdangl LT -2 then guess = 'PFLT'  ;; Pixel flat
      endif
  endif else begin ;; Calibration of some sort
      flg = 0
      case strtrim(sxpar(head,keywd.lamp)) of
          'none': begin
              if exp EQ 1 then guess = 'ZRO' else guess = 'DRK'
          end
          'ThAr2': guess = 'ARC'
          'ThAr1': guess = 'ARC'
          'quartz': flg = 2
          'quartz1': flg = 2
          'quartz2': flg = 2
          'undefined': guess = 'UNK'
          else: stop
      endcase
      ;; 
      if not keyword_set( GUESS ) then begin
          case flg of
              2: begin          ; QTZ lamp
                  if deck EQ 'D5' then guess = 'TRC' else guess='TFLT'
                  if xdangl LT -2. then guess = 'PFLT'  ;; Pixel flat
              end
              else: stop
          endcase
      endif
  endelse

  ;; Better have a guess now
  if not keyword_set( GUESS ) then stop
          

  ;; Final checks
  case strtrim(sxpar(head,keywd.obstyp),2) of
      'Object': user = 'OBJ'
      'IntFlat': user = 'TFLT'
      'SkyFlat': user = 'PFLT'
      'DmFlat': user = 'PFLT'
      'Line': user = 'ARC'
      'Bias': user = 'ZRO'
      'Dark': user = 'DRK'
      else: user = 'UNK' ;; For the Old Chip
  endcase
  if guess NE user then begin
      if guess NE 'TRC' and guess NE 'STD' $
        and not keyword_set(SILENT) then begin
          print, 'hires_headid: Warning, my guess does not match ' + $
            'your designation!'
          print, guess, ' NE ', user
      endif
  endif
  
  return, guess

end
    
