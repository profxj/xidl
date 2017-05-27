;+
; NAME:
;   read_lampfile
;
; PURPOSE:
;   Read in lamp list for with arcimage fitting routines
;   
; CALLING SEQUENCE:
;   lamps = read_lampfile(lampfilename)
;
; INPUTS:
;   lampfilename  - filename with full path.  For example, 
;        lampfilename = '~dfink/idlspec2d/etc/lamphgcdne.dat'
;  we should probably check this in to the deep repository under spec2d/etc
;
; REVISION HISTORY:
;   2001-Sep-18  Written by Finkbeiner
;   2002-May-24  read element MD

function read_lampfile, lampfilename

; Read lamp file into a structure
  readcol, lampfilename, lampwave, lampinten, lampquality, element, $
        format='D,F,A,A', /silent
  lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0B, $
                      element:'', quality:''}
  lamps = replicate(lamps, N_elements(lampwave))
  lamps.lambda = lampwave
  lamps.loglam = alog10(lampwave)
  lamps.intensity = lampinten
  lamps.good = strupcase(lampquality) EQ 'GOOD' AND lampinten GT 0
  lamps.element = element
  lamps.quality=strupcase(lampquality)


;hacked version to read simpler table
; readcol, lampfilename, lampwave, lampinten, format='D,F', /sil
;  lamps = {lambda: 0.0d0, loglam: 0.0d0, intensity: 0.0d0, good: 0B}
;  lamps = replicate(lamps, N_elements(lampwave))
;  lamps.lambda = lampwave
;  lamps.loglam = alog10(lampwave)
;  lamps.intensity = lampinten
;  lamps.good = 1B

  return, lamps
end

