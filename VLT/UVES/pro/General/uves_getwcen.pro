;+ 
; NAME:
; uves_getwcen
;     Version 1.1
;
; PURPOSE:
;   Pass back the central wavelength
;
; CALLING SEQUENCE:
;   
;  wcen = uves_getfil(fil)
;
; INPUTS:
;   [setup]   -  Setup identifier 
;
; RETURNS:
;  Structure, image, name, etc.
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /name   -- Only return resolved name (string)
;   CHKFIL  -- Value equal to the number of files matching name
;   SZ      -- Image size
;   SUBFIL  -- Image name generally used to parse the root name of the
;              image  (e.g.  'Arcs/arc_mb0539.fits').  Required in
;              many cases.
;   WCEN    -- Specify camera: Blue (1) or Red (2)
;   INDX    -- Image extension in the fits file (generally 0, 1, or 2)
;
; OPTIONAL OUTPUTS:
;   FIL_NM  -- Filename of the file 
;   HEAD    -- Image header
;
; COMMENTS:
;
; EXAMPLES:
;   ordr_str = uves_getfil('ordr_str', 1, WCEN=340)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------

function uves_getwcen, fil, uves, FRAME=frame, SARA=sara, AFIL=afil

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'wcen = uves_getfil(fil) [v1.0]'
      return, -1
  endif 

  if keyword_set(SARA) or keyword_set(AFIL) then begin
      ;; Header
      head = xheadfits(fil)
      
      ;; Parse
      wcen = round(sxpar(head,'WLEN1'))
      if wcen LE 0 then wcen = round(sxpar(head,'WLEN2'))

      if keyword_set(AFIL) then begin
          prs = strsplit(fil, '_',/extract)
          wcen = long(prs[1])
      endif
      
      if arg_present(frame) then begin
          prs = strsplit(fil, '_',/extract)
          nprs = n_elements(prs)
          pos = strpos(prs[nprs-1],'.fits')
          frame = fix(strmid(prs[nprs-1],0,pos))
      endif
  endif else begin
      a = where(strmatch(uves.rootpth+uves.img_root,fil),na)
      if na NE 1 then stop
      wcen = round(uves[a].xdangl)
      frame = uves[a].frame
  endelse

  return, wcen
end

