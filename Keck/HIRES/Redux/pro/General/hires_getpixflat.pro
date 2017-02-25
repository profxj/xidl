;+ 
; NAME:
; hires_getpixflat   
;     Version 1.1
;
; PURPOSE:
;  Identify the Pixelflat closest in UT in the HIRES_CALIBS directory
;
; CALLING SEQUENCE:
;   
;  rslt = hires_getpixflat(hires, flat, IFLAT=, FIL=)
;
; INPUTS:
;   hires -- HIRES structure
;
; RETURNS:
;
; OUTPUTS:
;   FLAT=  -- Data from pixel flat
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   IFLAT=  -- Inverse variance of pixel flat
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------

function hires_getpixflat, hires, flat, IFLAT=iflat, FIL=gdfil


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = hires_getpixflat(hires, pflat]) [v1.1]'
      return, -1
  endif 

  ;; Check to see if the flag is set
  if strlen(getenv('HIRES_CALIBS')) EQ 0 then begin
      print, 'hires_getpixflat: You need to set HIRES_CALIBS !!'
      return, -1
  endif

  ;; Root
  case hires.chip of
     -1: return, -1  ;; Nothing for the Single chip
      1: cnm = 'B'
      2: cnm = 'G'
      3: cnm = 'R'
  endcase
  
  root = 'hires_pix'+hires.cross+strtrim(hires.colbin,2)+'x'+ $
    strtrim(hires.rowbin,2)+cnm+'_'
  root = strcompress(root, /rem)
  ;; Files
  fil = findfile(getenv('HIRES_CALIBS')+'/FLATS/'+root+'*', count=nfil)
  if nfil EQ 0 then begin
      print, 'hires_getpixflat: No flats matching the setup! '
      stop
      return, -1
  endif

  ;; Files
  if nfil GT 1 then begin
      mjd = lonarr(nfil)
      ;; Parse date
      for qq=0L,nfil-1 do begin
          pos = strpos(fil[qq], '.fits')
          dt = strmid(fil[qq], pos-9,9)
          mjd[qq] = x_setjdate(dt)
      endfor
      mn = min(abs(mjd-hires.date),imn)
      gdfil = fil[imn]
  endif else gdfil = fil[0]

  if arg_present(flat) then flat = xmrdfits(gdfil, /silent)
  if arg_present(iflat) then iflat = xmrdfits(gdfil, 1, /silent)
      
    
  return, 1
end

