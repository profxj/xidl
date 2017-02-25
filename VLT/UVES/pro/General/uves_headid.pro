;+ 
; NAME:
; uves_headid   
;     Version 1.1
;
; PURPOSE:
;    Guess the type of MIKE image based on exposure time, counts,
;    and other characteristics of the image.
;
; CALLING SEQUENCE:
;   
;  guess = uves_headid(data, xbinguess, ybinguess, $
;               filename=img[q], hdr=head, /silent)
;
; INPUTS:
;    data      - 2D image (generally set to 0 and read from filename)
;
; RETURNS:
;  Image type::  'UNK', 'BAD', 'ZRO', 'ARC', 'TFLT', 'TFLT', 'STD',
;  'OBJ', 'MFLT', 'TWI'
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FILENAME   - Name of the image to be considered
;   XBIN       - x binning
;   YBIN       - y binning
;   SILENT     - Suppress print statements
;   SATURATED  - Saturation level (default = 50000.)
;   
;
; OPTIONAL OUTPUTS:
;   HDR     -  Image header if filename specified
;   GUESS_EXPTIME -- Guess of exposure time based on CR hits
;
; COMMENTS:
;   Need to do some consistency checks between red and blue side,
;         i.e. they should be the same type
;
; EXAMPLES:
;   guess = uves_headid(data, xbin, ybin, filename='file')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------
function uves_headid, head, side, SILENT=silent, SARA=sara, H2=h2

  ;;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'guess = uves_headid(header, side, /SILENT) [v1.1]'
      return, -1
  endif 

  nhead = n_elements(head)

  ;; Sara naming convention
  if keyword_set(SARA) then begin
      ;; Assumes Raw/
      case strmid(sara,4,5) of
          'carc_': guess = 'ARC'
          'cbias': guess = 'ZRO'
          'cflat': guess = 'TFLT'
          else: begin
              if strmid(sara,4,1) NE 'c' then guess = 'OBJ' else guess = 'UNK'
          end
      endcase

  endif else begin
      obj = sxpar(head, 'OBJECT')
      if strpos(obj, 'FLAT') GE 0 then guess = 'TFLT'
      if strpos(obj, 'WAVE') GE 0 then guess = 'ARC'
      if not keyword_set(GUESS) then guess = 'OBJ'
  endelse

  ;; Modify the Header

  ;; Lamps
  mt = where(strmid(head,0,27) EQ 'HIERARCH ESO INS LAMP1 NAME',nmt1)
  if nmt1 NE 0 then head[mt] = 'LAMPNAME='+strmid(head[mt],31)
  mt = where(strmid(head,0,27) EQ 'HIERARCH ESO INS LAMP2 NAME',nmt2)
  if nmt2 NE 0 then head[mt] = 'LAMPNAME='+strmid(head[mt],31)
  mt = where(strmid(head,0,27) EQ 'HIERARCH ESO INS LAMP3 NAME',nmt3)
  if nmt3 NE 0 then head[mt] = 'LAMPNAME='+strmid(head[mt],31)
  mt = where(strmid(head,0,27) EQ 'HIERARCH ESO INS LAMP4 NAME',nmt4)
  if nmt4 NE 0 then head[mt] = 'LAMPNAME='+strmid(head[mt],31)
  mt = where(strmid(head,0,27) EQ 'HIERARCH ESO INS LAMP5 NAME',nmt5)
  if nmt5 NE 0 then head[mt] = 'LAMPNAME='+strmid(head[mt],31)
  if (nmt1 + nmt2 + nmt3 + nmt4 + nmt5) EQ 0 then $
    head[nhead-5] = 'LAMPNAME= ''NONE'''
  
  if side EQ 1 then begin ; BLUE
      ;; SLIT WIDTH 
      mt = where(strmid(head,0,26) EQ 'HIERARCH ESO INS SLIT2 WID',nmt1)
      if nmt1 EQ 1 then head[mt] = 'SLITWID ='+strmid(head[mt],31) $
      else head[nhead-6] = 'SLITWID = 0.'

      ;; Central wavelength
      mt = where(strmid(head,0,27) EQ 'HIERARCH ESO INS GRAT1 WLEN',nmt1)
      head[mt] = 'WLEN1   ='+strmid(head[mt],31)
      
      ;; Binning (Blue)
      mt = where(strmid(head,0,26) EQ 'HIERARCH ESO DET WIN1 BINX',nmt1)
      head[mt] = 'BINX    ='+strmid(head[mt],30)
      mt = where(strmid(head,0,26) EQ 'HIERARCH ESO DET WIN1 BINY',nmt1)
      head[mt] = 'BINY    ='+strmid(head[mt],30)
      
      ;; Pre/Over scan
      mt = where(strmid(head,0,27) EQ 'HIERARCH ESO DET OUT1 PRSCX',nmt1)
      head[mt] = 'PRESCAN ='+strmid(head[mt],31)
      mt = where(strmid(head,0,27) EQ 'HIERARCH ESO DET OUT1 OVSCX',nmt1)
      head[mt] = 'OSCAN   ='+strmid(head[mt],31)

      ;; Gain, RN
      mt = where(strmid(head,0,25) EQ 'HIERARCH ESO DET OUT1 RON',nmt1)
      head[mt] = 'CCDRN   ='+strmid(head[mt],30)
      mt = where(strmid(head,0,26) EQ 'HIERARCH ESO DET OUT1 GAIN',nmt1)
      head[mt] = 'CCDGAIN ='+strmid(head[mt],30)

  endif

  if side EQ 2 then begin ; RED
      ;; SLIT WIDTH 
      mt = where(strmid(head,0,26) EQ 'HIERARCH ESO INS SLIT3 WID',nmt1)
      if nmt1 EQ 1 then head[mt] = 'SLITWID ='+strmid(head[mt],31) $
      else head[nhead-6] = 'SLITWID = 0.'

      ;; Header 2
      if not keyword_set(h2) then h2 = head  ;; Kludge
      
      ;; Binning (Blue)
;      mt = where(strmid(head,0,26) EQ 'HIERARCH ESO DET WIN1 BINX',nmt1)
      head[mt+10] = 'BINX    ='+strtrim(round(sxpar(h2,'cdelt1')),2)
;      mt = where(strmid(head,0,26) EQ 'HIERARCH ESO DET WIN1 BINY',nmt1)
      head[mt+11] = 'BINY    ='+strtrim(round(sxpar(h2,'cdelt2')),2)
      
      ;; Pre/Over scan
      mt = where(strmid(h2,0,27) EQ 'HIERARCH ESO DET OUT1 PRSCX',nmt1)
      head[mt+12] = 'PRESCAN ='+strmid(h2[mt],31)
      mt = where(strmid(h2,0,27) EQ 'HIERARCH ESO DET OUT1 OVSCX',nmt1)
      head[mt+13] = 'OSCAN   ='+strmid(h2[mt],31)

      ;; Gain, RN
      mt = where(strmid(h2,0,25) EQ 'HIERARCH ESO DET OUT1 RON',nmt1)
      head[mt+14] = 'CCDRN   ='+strmid(h2[mt],30)
      mt = where(strmid(h2,0,26) EQ 'HIERARCH ESO DET OUT1 GAIN',nmt1)
      head[mt+15] = 'CCDGAIN ='+strmid(h2[mt],30)
  endif
  
  return, guess

end
    
