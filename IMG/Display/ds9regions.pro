;+
; NAME:
;   ds9regions
;
; PURPOSE:
;   Creates a region file readable by ds9, 
;   with region positions specified by RA and 
;   DEC in either decimal degrees, or sexigesimal 
;   notation.  **This program will only be useful if
;   you have an image with a WCS header, which you can
;   display in ds9.**
;
; CALLING SEQUENCE:
;   ds9regions, filename, ra, dec, shape=shape, $
;   boxWidth=boxWidth, boxHeight=boxHeight, boxAngle=boxAngle, $
;   circleRadius=circleRadius, color=color, text=text, $
;   /MOVEABLE, /RESIZEABLE, /SEXIGESIMAL, /B1950, /APPEND
;
; INPUTS:
;   filename   - String specifying the name of the region file to be created.
;   ra, dec    - RA and DEC in decimal degrees (double).
;                Default is J2000 epoch coordinates.
;                To specify more than one region, ra and dec should be arrays.
;                RA and DEC can be in sexigesimal format if the /SEXIGESIMAL flag is set.
;                Then: RR:RR:RR.R -DD:DD:DD.D
;
; OPTIONAL INPUTS:
;   shape        - String (or string array) specifying ds9 region shape.
;                  Currently supports the following region shapes: 
;                  'box', 'circle', 'circle point', 'x point', 'cross point'
;                  If no shape is given, defaults to circle point.
;   boxWidth     - If shape is 'box', then boxWidth specifies the width of the box in arcsec.
;                  Default boxWidth is 1".
;   boxHeight    ---    Specifies height of box in arcsec.  Default boxHeight is 100".
;   boxAngle     ---    Specifies orientation of box in degrees (N to E).  Default orientation is 0 degrees.
;   circleRadius - If shape is 'circle', then circleRadius specifies the radius of circle in arcsec.  Default is 5".
;   color        - String specifying color of each region.  If no color is given, defaults to green.
;                  All ds9 colors are supported: black, white, red, green, blue, cyan, magenta, yellow
;   text         - String of text to be written next to each region.
;
; OPTIONAL KEYWORDS:
;   MOVEABLE     - The regions will all be moveable in ds9 when loaded.  Otherwise, regions cannot be deleted, resized, moved, etc.
;   RESIZEABLE   - The regions will all be resizeable in ds9 when loaded.
;   SEXIGESIMAL  - Allows RA and DEC to be given as RR:RR:RR.R -DD:DD:DD.D strings.
;   B1950        - RA and DEC specified in B1950 coordinates rather than J2000 coords.
;   APPEND       - Instead of writing a new ds9 region file, appends additional
;                  regions to existing file named 'filename.'
;
; OUTPUT:
;   Writes to a new file named 'filename.'
;
; COMMENTS:
;   If more than one position is specified, but optional inputs are of length 1, 
;   then all positions use the same options.  (See Examples)
;
;   Use ds9 to load the region file on top of a fits image like this:
;   ds9 image.fits -regions filename &
;
; EXAMPLES:
;   ds9regions, 'filename', [273.0599, 273.1922], [-3.20993, -3.21107], shape='circle', $
;   circleRadius=3, color=['yellow','red'], text=['text1', 'text2'], /MOVEABLE
;
;   ds9regions, 'filename', [273.0599, 273.1922], [-3.20993, -3.21107], shape='box', $
;   boxHeight=[175,200]
;
;   ds9regions, 'filename', ['12:42:18.02', '12:49:18.399'], ['-2:15:47.11', '-2:19:21.87'], /sexigesimal, $
;   shape=['x point', 'cross point']
;
; BUGS:
;   -- If using B1950 coordinates, make sure not to mix with J2000 coordinates when using /APPEND.
;
; PROCEDURES CALLED:
;   x_radec()
;
; REVISION HISTORY:
;   14-Mar-2007  Written by L. Pollack
;-

PRO ds9regions, filename, ra, dec, shape=shape, boxWidth=boxWidth, boxHeight=boxHeight, boxAngle=boxAngle, circleRadius=circleRadius, color=color, text=text, B1950=b1950, APPEND=append, SEXIGESIMAL=sexigesimal, RESIZEABLE=resizeable, MOVEABLE=moveable

  if (N_params() lt 3) then begin
     print, 'Syntax - ' + $
            "ds9regions, filename, ra, dec, shape=shape, boxWidth=boxWidth, boxHeight=boxHeight, boxAngle=boxAngle, $"
     print, "         circleRadius=circleRadius, color=color, text=text, /B1950, /APPEND, /SEXIGESIMAL, /RESIZEABLE, /MOVEABLE"
     return
  endif
  
  if not(keyword_set(APPEND)) then openw, 1, filename else openw, 1, filename, /append
  
  ;; Default regions are green, and cannot be edited or deleted once loaded into ds9.
  if not(keyword_set(APPEND)) then begin
     printf, 1, '# Region file format: DS9 version 4.0'
     printf, 1, 'global color=green font="helvetica 10 normal" select=1 highlite=1 edit=0 move=0 delete=0 rotate=0 include=1 fixed=0 source'
     if keyword_set(B1950) then printf, 1, 'fk4' else printf, 1, 'fk5'
  endif

  for i=0, n_elements(ra)-1 do begin
     if keyword_set(SEXIGESIMAL) then begin
        ra_sex=ra[i]
        dec_sex=dec[i]
        x_radec, ra_sex, dec_sex, ra1, dec1
        ra1=strcompress(string(ra1), /rem)
        dec1=strcompress(string(dec1), /rem)
     endif else begin
        ra1=strcompress(string(double(ra[i])), /rem)
        dec1=strcompress(string(double(dec[i])), /rem)
     endelse

     if keyword_set(SHAPE) then begin
        shape1 = (n_elements(SHAPE) lt n_elements(ra)) ?  shape[0] : shape[i]
     endif else begin
        shape1 = 'circle point'  ;defaults to circle point region shape
     endelse

     if shape1 eq 'circle' then begin
        if keyword_set(circleRadius) then begin
           circleRadius1 = (n_elements(circleRadius) lt n_elements(ra)) ? circleRadius[0] : circleRadius[i]
           circleRadius1 = strcompress(string(circleRadius1), /rem)
        endif else begin
           circleRadius1 = '5'  ;if circle shape is specified with no size, defaults to 5" circle radius
        endelse
     endif

     if shape1 eq 'box' then begin
        if keyword_set(boxWidth) then begin
           boxWidth1 = (n_elements(boxWidth) lt n_elements(ra)) ? boxWidth[0] : boxWidth[i]
           boxWidth1 = strcompress(string(boxWidth1), /rem)
        endif else begin
           boxWidth1 = '1'  ;if box shape is specified with no width, defaults to 1" width
        endelse
        if keyword_set(boxHeight) then begin
           boxHeight1 = (n_elements(boxHeight) lt n_elements(ra)) ? boxHeight[0] : boxHeight[i]
           boxHeight1 = strcompress(string(boxHeight1), /rem)
        endif else begin
           boxHeight1 = '100'  ;if box shape is specified with no height, defaults to 100" tall
        endelse
        if keyword_set(boxAngle) then begin
           boxAngle1 = (n_elements(boxAngle) lt n_elements(ra)) ? boxAngle[0] : boxAngle[i]
           boxAngle1 = strcompress(string(boxAngle1), /rem)
        endif else begin
           boxAngle1 = '0'  ;if box shape is specified with no orientation, defaults to an angle of 0 degrees.
        endelse
     endif

     if keyword_set(COLOR) then begin
        color1 = (n_elements(COLOR) lt n_elements(ra)) ? color[0] : color[i]
     endif
     if keyword_set(TEXT) then begin
        text1 = (n_elements(TEXT) lt n_elements(ra)) ? text[0] : text[i]
     endif

     general_shape = (strpos(shape1, 'point') ne -1) ? 'point' : shape1
     if general_shape eq 'point' then shape1 = strmid(shape1, 0, strpos(shape1, 'point')-1)

     command=general_shape+'('+ra1+','+dec1
     if general_shape eq 'box' then command=command+','+boxWidth1+'",'+boxHeight1+'",'+boxAngle1
     if general_shape eq 'circle' then command=command+','+circleRadius1+'"'
     command=command+') # '

     if general_shape eq 'point' then command=command+'point='+shape1+' '
     if keyword_set(COLOR) then command=command+'color='+color1+' '
     if keyword_set(TEXT) then command=command+'text={'+text1+'} '
     if keyword_set(RESIZEABLE) then command=command+'edit=1 '
     if keyword_set(MOVEABLE) then command=command+'move=1'

     doit="printf, 1, '"+command+"'"
     void=execute(doit)

endfor
close, 1

END
