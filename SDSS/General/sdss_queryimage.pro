;+
; NAME:
; sdss_queryimage
;
;PURPOSE
;	to query SDSS image
;SYNTAX
;	sdss_queryimage, ra, dec [, scale , xs=xs, ys=ys, /gridoff, /labeloff]
;INPUTS
;	ra: right ascension (J2000 decimal degrees)
;	dec: declination (J2000 decimal degrees)
;	scale: plate scale (arcsec/pixel) [default to .3961]
;KEYWORDS
;	filename: name of filt to be written to [default to temp.jpeg]
;	xs: xsize (in pixels) [default to 512]
;	ys: y size (in pixels) [default to 512]
;	labeloff: if set leaves off label
;	gridoff: if set leaves off Gid
; RETURNS:  jpeg of the image
;
;OPTIONAL OUTPUTS
;	FILENAME= writes a file to the filename given by the keyword
;	
;NOTES
;	by default gives you pretty picture of blue spiral galaxy (default SDSS image)
;	with ra=18.87667, dec=-0.86083
;Written by R. da Silva, UCSC, 3-11-09
;-
;----------------------------------------------------------------------

function sdss_queryimage, ra1, dec1, scale, filename=filename, $
  label=label,  grid=grid, xs=xs, ys=ys, GRAY=gray, INVERT=invert


  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
            'img = sdss_queryimage, RA, DEC, [scale], FILENAME=, ' + $
            '/GRAY, /LABEL, /GRID [v1.1]'
      return, -1
  endif 

  ;--- Setting Defaults
;  if NOT keyword_set(ra1) then ra1=18.87667
;  if NOT keyword_set(dec1) then dec1=-0.86083
  if NOT keyword_set(scale) then scale=.3961
  if NOT keyword_set(xs) then xs=512
  if NOT keyword_set(ys) then ys=512

;---Making the HTML
;name1='http://casjobs.sdss.org/ImgCutoutDR7/'
name1='http://skyservice.pha.jhu.edu/DR10/ImgCutout/'
name='getjpeg.aspx?ra='
name+=strcompress(string(ra1), /remove_all) 	;setting the ra
name+='&dec='
name+=strcompress(string(dec1), /remove_all) 	;setting the declination
name+='&scale='
name+=strcompress(string(scale), /remove_all) 	;setting the scale
name+='&width='
name+=strcompress(string(xs), /remove_all)	;setting the width
name+='&height='
name+=strcompress(string(ys), /remove_all)	;setting the height

;------ Options
options = ''
if keyword_set(grid) then options+='G'
if keyword_set(label) then options+='L'
if keyword_set(invert) then options+='I'
if strlen(strtrim(options,2)) GT 0 then name+='&opt='+options

name+='&query='

;--- Spawning the command 
spawn, '\wget'+ ' "' +name1+name + '"', temp
if NOT keyword_set(filename) then filename='temp.jpeg'
;--- Renaming file
spawn, '\mv "' + name +'" ' +filename

read_jpeg, filename, img, GRAYSCALE=gray

return, img

end
