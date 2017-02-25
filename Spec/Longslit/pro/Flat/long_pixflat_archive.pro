;+ 
; NAME:
; long_pixflat_archive   
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
;   26-Jan-2010 Written by JXP
;-
;------------------------------------------------------------------------------

pro long_pixflat_archive, strct, flatnm, pixct, ROOT=root


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'long_pixflat_archive, strct, flatnm [v1.1]'
      return
  endif 

  if not keyword_set(ROOT) then root = ''

  ;; Check to see if the flag is set
  if strlen(getenv('LONGSLIT_DIR')) EQ 0 then begin
      print, 'long_pixflat_archive: You need to set LONGSLIT_DIR !!'
      stop
  endif
  pixct = 0L
  flatnm = ''

  case strct.instrument of
      'LRISBLUE': begin
         path = getenv('LONGSLIT_DIR')+'calib/flats/LRIS/'
         ipos = strpos(strct.grating,'/')
         dispers = 'B'+strmid(strct.grating,0,ipos)
         head = xheadfits(ROOT+strct.filename)
         binning = long(strsplit(sxpar(head, 'BINNING'), ',', /extract))
         cbin = strtrim(binning[0],2)+'x'+strtrim(binning[1],2)

         ;; Grab the files
         pixfil = findfile(path+'pixflat_'+dispers+'_'+cbin+'_*', count=nfil)
         case nfil of
             0: return
             1: begin
                 flatnm = pixfil[0]
                 pixct = 1
                 return
             end
             else: begin
                 mjd = lonarr(nfil)
                 ;; Parse date
                 for qq=0L,nfil-1 do begin
                     pos = strpos(pixfil[qq], '.fits')
                     dt = strmid(pixfil[qq], pos-9,9)
                     mjd[qq] = x_setjdate(dt)
                 endfor
                 mjds = x_setjdate(strmid(sxpar(head,'DATE'),0,10))
                 mn = min(abs(mjd-mjds),imn)
                 flatnm = pixfil[imn]
                 pixct = 1
                 return
             end
         endcase
      end
      else: return
  endcase

  return
end

