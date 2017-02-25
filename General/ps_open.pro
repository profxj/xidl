;
;+
; Name:
;
;   PS_Open
;
; Purpose:
;
;   Set the plotting device to be PostScript, 
;   and save the current plotting device (that can thus be restored).
;   Use in conjunction with ps_close.

; Calling sequence:
;
;   PS_Open
;
; Example:
;   
;   IDL> ps_open,file='test.ps',/por
;   IDL> plot,findgen(10)
;   IDL> ps_close,/nop
;
; Keywords:
;
;   PORTRAIT - Normally, PS_Open defines the plotting area as landscape. If
;     this keyword is set, the plotting area will be portrait.
;
;   FILENAME = <string> - Normally, PS_Open uses '/tmp/idl-$USER.ps' (on
;     Sparx) for the PostScript file, but the keyword FILENAME can be used
;     to override that value.  ($USER is retrieved with a getenv('USER')).
;
;   SQUARE - if set, defines the plotting area to be square.  to also get
;     a square data window, use the position keyword when making the plot.
;
;   FONT - set !P.font to value of FONT
;
;   ENCAPSULATED - produces an encapsulated PostScript plot.
;   
;   BPP = <n> - Set the number of bits per pixel; default value is 8,
;     valid values are 8, 4, OR 2. A higher no produces a better
;     resolution but also a bigger PostScript file.
;
;   COLOR - produces a color PostScript plot.
;
;   LEDGER - selects ledger format (17x11; usefull mostly for color printer)
;
;   XSIZE - the x size in inches of the plotting region
;
;   YSIZE - the y size in inches of the plotting region
;
;   MAXS -- Maximize the plot on a standard 8.5 x 11 page (JXP)
;
;   CMYK -- generate colors in CMYK for publication
;
; See also:
;    PS_Close
;
; History:
;   23-Aug-93: added bpp=[2,4,8] option, with 8 bpp as default
;   24-Aug-93: added support for /color, 
;              and reduced the page size for QMS color printer
;   16-Sep-93: added color flag to common block (see also PS_Close)
;    9-Apr-94: added file name to info message
;   11-Apr-94: added LEDGER support (see also PS_Close)
;   28-Apr-94: added SQUARE, and fixed color portrait offset
;    5-May-94: got /SQUARE to actually work
;   25-Jun-94: added xsize and ysize keywords
;    3-Dec-07  added CMYK, KLC
;-
;
PRO ps_open, portrait = orient, filename = fn, font = font,  $
             square=square_plot,  encapsulated = eps, $
             color = cps, xsize = xs, ysize = ys, $
             bpp = nbpp, ledger = ldgr, MAXS=maxs, CMYK=cmyk, silent=silent

  COMMON ps_common, old_dname, old_pfont, file_name, opened, color, ledger
;
  if keyword_set( MAXS ) then begin
      xs= 10.8
      ys= 8.2
  endif
  IF n_elements(opened) EQ 0 THEN opened = 0
;
  IF NOT keyword_set(orient) THEN orient = 0
  IF keyword_set(ldgr) THEN ledger = 1 ELSE ledger = 0
;
  IF ledger EQ 0 THEN $
    sizes = [7.5, 10.] $
  ELSE BEGIN
    sizes = [16., 10.]
    orient = 1 - orient
  ENDELSE
;
  xoffset = 0.5
  yoffset = 0.5
  CASE orient OF
    0: BEGIN
      xsize = sizes(1)
      ysize = sizes(0)
    END
    1: BEGIN
      xsize = sizes(0)
      ysize = sizes(1)
    END
    Else: BEGIN
      print, 'Invalid orientation: ', orient
      return
    END
  ENDCASE
;
  IF opened EQ 1 THEN BEGIN
     if ~keyword_set(silent) then message, /info, $
                                           'WARNING: device already opened to PS, closing it first.' 
     if ~keyword_set(silent) then message, /info, $
                                           '         any plot in progress will be lost.'
  ENDIF ELSE BEGIN
    old_dname = !D.name
    old_pfont = !P.font
  ENDELSE
  set_plot, 'PS', /copy
  device, /close
  opened = 1
  IF NOT keyword_set(fn) THEN fn = '/tmp/idl-'+getenv('USER')+'.ps'
  IF keyword_set(eps) THEN device, /encapsulated, filename = fn, cmyk = cmyk ELSE $
    device, filename = fn, cmyk = cmyk
;
  file_name = fn
;
  color = 0
  IF keyword_set(cps) THEN BEGIN
    device, /color
    IF orient EQ 0 THEN BEGIN
      xsize = xsize - 1.0 
      xoffset = xoffset + 0.8
;      xoffset = xoffset + 3.3
    ENDIF ELSE BEGIN
      ysize = ysize - 1.0
      yoffset = yoffset + 0.8
    ENDELSE
    color = 1
  ENDIF
;
  IF keyword_set(square_plot) THEN BEGIN
    xsize =  min(sizes)
    ysize =  xsize
  ENDIF
;
  IF orient EQ 0 THEN device, /landscape, $
    /inches, xsize = xsize, ysize = ysize, yoffset=9.4 $ ;,  yoffset=0.5 $
  ELSE device, /portrait, $
    /inches, xsize = xsize, ysize = ysize, xoffset = xoffset, yoffset = yoffset
;
  IF keyword_set(xs) and not keyword_set(ys) then ys = ysize
  IF keyword_set(ys) and not keyword_set(xs) then xs = xsize 
  IF keyword_set(xs) and keyword_set(ys) and not keyword_set(enc) then $
    IF orient EQ 0 THEN device, /landscape, /inches, xsize = xs, ysize = ys, $
      yoffset = 5.5 + .5*xs, xoffset = 4.25 - .5*ys $
    ELSE device, /portrait, /inches, xsize = xs, ysize = ys, $
      xoffset = 4.25 - .5*xs, yoffset = 5.5 - .5*ys
;
  IF keyword_set(nbpp) THEN device, bits_per_pixel = nbpp $
  ELSE device, bits_per_pixel = 8
;
  IF keyword_set(font) THEN !p.font = font
;    
  if ~keyword_set(silent) then message, /info, 'output redirect to PostScript file ' + fn
  IF keyword_set(square_plot) THEN begin
     if ~keyword_set(silent) then message, /info, 'using a square area'
  endif
;
  return
END
