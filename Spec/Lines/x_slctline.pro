;+ 
; NAME:
; x_slctline   
;   Version 1.1
;
; PURPOSE:
;    Launches a GUI and allows the user to select a line.  The code
;    then returns the rest wavelength of that line.
;
; CALLING SEQUENCE:
;   line = x_slctline, lines, /ISM, NLIN=, ILIN=, XOFFSET=, YOFFSET=, /GAL,
;                    PDMENU=
;
; INPUTS:
;   lines  - Line list structure
;
; RETURNS:
;   line - Wavelength (double)
;
; OUTPUTS:
;  line = unit number of line
;
; OPTIONAL KEYWORDS:
;   /ISM     - ISM style
;   /GAL     - Use the Galaxy style
;   YOFFSET= - Placement of the GUI [default: 200]
;   XOFFSET= - Placement of the GUI [default: 200]
;   PDMENU=  - String array formatted for cw_pdmenu
;   NLIN=    - Grab only the first NLIN lines (only for non-ISM or
;              GAL)
;
; OPTIONAL OUTPUTS:
;   ILIN= - Index of the selected line in the line list structure
;
; COMMENTS:
;
; EXAMPLES:
;   line = x_slctline( lines )
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;   29-Nov-2001 Significantly revised
;   08-Dec-2001 Allows ism or name+wave
;-
;------------------------------------------------------------------------------

function x_slctline_setpdmenu_ism, lines
;

  nlin = n_elements(lines)
  wavnm = strmid(strtrim(lines.wave, 2), 0, 7)
  trmnms = strtrim(lines.name, 2)

; Find all unique elements first

  imrk = strpos(lines.name, 'I', 1) 
  vmrk = strpos(lines.name, 'V', 1)
  xmrk = strpos(lines.name, 'X', 1)

  imrk[where(imrk LT 0)] = 999L
  vmrk[where(vmrk LT 0)] = 999L
  xmrk[where(xmrk LT 0)] = 999L

  mrk = lonarr(nlin)

  ;; Find the elm
  for i=0L, nlin-1 do mrk[i] = (imrk[i] < vmrk[i]) < xmrk[i]

  nms = strarr(nlin)
  for i=0L,nlin-1 do nms[i] = strtrim(strmid(lines[i].name, 0, mrk[i]),2)
      
  delvarx, imrk, vmrk, mrk

  elem = x_uniqstr(nms, count=nelem, /sort)

; Create pdmenu

  for ww=0,nelem-1 do begin

      ; First string
      if ww NE nelem-1 then tmps = '1\'+elem[ww] else  tmps = '3\'+elem[ww] 
      if ww EQ 0 then pd = [tmps] else pd = [pd, tmps]

      ; Find all ions
      lenelm = strlen(elem[ww])

      allion = where( nms EQ elem[ww] AND $
                      (strmid(trmnms, lenelm, 1) EQ 'I' $
                       OR strmid(trmnms, lenelm, 1) EQ 'V' $
                       OR strmid(trmnms, lenelm, 1) EQ 'X'), $
                      nallion)
      mrk = strpos( trmnms[allion], ' ', 2)

      ; Parse out the Ion
      allion_nm = strarr(nallion)
      for i=0L,nallion-1 do allion_nm[i] = $
        strmid(trmnms[allion[i]], lenelm, mrk[i]-lenelm)

      ; Save the wavelengths

      svwav = wavnm[allion]
      ions  = x_uniqstr(allion_nm, count=nion, /sort)

      ; Loop on the Ions
      for nn=0L,nion-1 do begin
          
          ; Add the ion entry
          if nn EQ nion-1 then pd = [pd, '3\'+ions[nn]] else $
            pd = [pd, '1\'+ions[nn]] 

          ; Find all wavelengths
          allwv = where(allion_nm EQ ions[nn], nwv)

          for qq=0,nwv-1 do begin
              if qq EQ nwv-1 then pd = [pd, '2\'+svwav[allwv[qq]]] $
              else pd = [pd, '0\'+svwav[allwv[qq]]] 
          endfor

      endfor
  endfor

  return, pd

end

;cccccccccccccccccccccccccccccccccccccccccccccccccc

function x_slctline_setpdmenu_gal, lines
;

  nlin = n_elements(lines)
  wavnm = strmid(strtrim(lines.wave, 2), 0, 7)

; Find all unique elements first


  nms = strtrim(lines.elm, 2)
  elem = x_uniqstr(nms, count=nelem, /sort)

; Create pdmenu

  for ww=0,nelem-1 do begin

      ; First string
      if ww NE nelem-1 then tmps = '1\'+elem[ww] else  tmps = '3\'+elem[ww] 
      if ww EQ 0 then pd = [tmps] else pd = [pd, tmps]

      ; Find all ions
      lenelm = strlen(elem[ww])
      allion = where( nms EQ elem[ww],  nallion)
      allion_nm = strtrim(lines[allion].ion, 2) 
      ions  = x_uniqstr(allion_nm, count=nion, /sort)

      ; Save the wavelengths
      svwav = wavnm[allion]

      ; Loop on the Ions
      for nn=0L,nion-1 do begin
          
          ; Add the ion entry
          if nn EQ nion-1 then pd = [pd, '3\'+ions[nn]] else $
            pd = [pd, '1\'+ions[nn]] 

          ; Find all wavelengths
          allwv = where(allion_nm EQ ions[nn], nwv)

          for qq=0,nwv-1 do begin
              if qq EQ nwv-1 then pd = [pd, '2\'+svwav[allwv[qq]]] $
              else pd = [pd, '0\'+svwav[allwv[qq]]] 
          endfor

      endfor
  endfor

  return, pd

end
      
;;;;;;;;;;;;;;;
; STD  [name, wave]

function x_slctline_setpdmenu, lines, NLIN=nlin

;

  if not keyword_set( NLIN ) then nlin = n_elements(lines)

  wavnm = strmid(strtrim(lines[0:nlin-1].wave, 2), 0, 7)
  nms = strtrim(lines[0:nlin-1].name, 2)

; Find all unique elements first
  elem = x_uniqstr(nms, count=nelem, /sort)

; Create pdmenu
  for ww=0,nelem-1 do begin
      ; First string
      if ww NE nelem-1 then tmps = '1\'+elem[ww] else  tmps = '3\'+elem[ww] 
      if ww EQ 0 then pd = [tmps] else pd = [pd, tmps]

      allwav = where( nms EQ elem[ww], nallwav)
      ; Do 50 at a time
      for i=0L,100 do begin
          i1 = i*40
          i2 = ((i+1)*40 - 1) < (nallwav-1)
          if i2 EQ nallwav-1 then pd = [pd, '3\'+wavnm[allwav[i1]]+$
                                        ' : '+wavnm[allwav[i2]]] $
          else pd = [pd, '1\'+wavnm[allwav[i1]]+$
                     ' : '+wavnm[allwav[i2]]]
          for qq=i1,i2 do begin
              if qq EQ i2 then pd = [pd, '2\'+wavnm[allwav[qq]]] $
              else pd = [pd, '0\'+wavnm[allwav[qq]]] 
          endfor
          if i2 EQ nallwav-1 then break
      endfor
  endfor

  ; FREE YOUR MIND
  delvarx, allwav, nms, wavnm
  return, pd

end
      
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_slctline, lines, ISM=ism, NLIN=nlin, ILIN=ilin, $
                     XOFFSET=xoffset, YOFFSET=yoffset, GAL=gal, $
                     PDMENU=pdmenu

;common x_slctline_ans

;
  if  N_params() NE 1  then begin 
    print,'Syntax - ' + $
             'line = x_slctline(lines, /ISM, NLIN=, XOFFSET=, YOFFSET='
    print,'        ) [V1.1]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( XOFFSET ) then xoffset = 0
  if not keyword_set( YOFFSET ) then yoffset = 0
  if not keyword_set( FONT ) then font = '5x7'

;    
;  x_slctline_initcommon

  ;; Multiple screens?
  oInfo = OBJ_NEW('IDLsysMonitorInfo')
  numMons = oinfo->GetNumberOfMonitors()
  names = oinfo->GetMonitorNames()
  rects = oInfo->GetRectangles()
  primaryIndex = oInfo->GetPrimaryMonitorIndex()
  rh_screen = 0
  OBJ_DESTROY, oInfo

;    WIDGET
  rh_Rect = rects[*, rh_screen]
  ;stop
  ;primaryRect = rects[*, primaryIndex]
  widgLoc = rh_Rect[0:1] ;+ rh_Rect[2:3] / 2
  base = WIDGET_BASE( title = 'x_slctline', /column, $
                      display_name=names[rh_screen], $ ;; Aims for right hand screen JXP 22 Aug 2014
                      XOFFSET=widgLoc[0], $
                      YOFFSET=widgLoc[1])

; Create pd menu

  if not keyword_set( PDMENU ) then begin
      flg = 0
      if keyword_set( ISM ) then flg = 1
      if keyword_set( GAL ) then flg = 2
      case flg of 
          1: menu_desc = x_slctline_setpdmenu_ism(lines) 
          2: menu_desc = x_slctline_setpdmenu_gal(lines) 
          else: menu_desc = x_slctline_setpdmenu(lines, NLIN=nlin) 
      endcase
      pdmenu = menu_desc
  endif else menu_desc = pdmenu

  
  ;; PD Lists
  list_id = cw_pdmenu(base, menu_desc, $
                      font=font, $
                      /help, $
                      /return_name, $
                      uvalue = 'LIST')
; Realize
  WIDGET_CONTROL, base, /realize
  
; LOOP
  svwav = 0.d
  repeat begin
      ev = widget_event(base)
      svwav = double(ev.value)
  end until svwav NE 0.d
  
  WIDGET_CONTROL, /destroy, base

; Send to the xmanager
;  xmanager, 'x_slctline', base, NO_BLOCK=noblock

; FIND THE RIGHT WAV
  diff = abs(lines.wave - svwav)
  tmp = min(diff, imn)
  delvarx, diff, tmp

  if arg_present( ILIN ) then ilin = imn

  return, lines[imn].wave
end
