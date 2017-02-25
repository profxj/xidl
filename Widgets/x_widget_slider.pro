;+ 
; NAME:
; xgetx_plt
;    Version 1.1
;
; PURPOSE:
;  Returns the x-value corresponding to the x-pixel on a GUI.
;   Requires pos, xymnx or the appropriate Structure.
;
; CALLING SEQUENCE:
;  xval = xgetx_plt(xtmp, pos, xymnx, size, /STRCT)
;   
; INPUTS:
;  xtmp -- x-pixel value
;  pos  -- Fraction of plot window covered by plot (2 element array)
;  xymnx -- x-y limits of plot window (4 element array: x0,y0,x1,y1)
;
; RETURNS:
;   xval -- 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STRCT -- xtmp contains a structure with the relevant tags
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uniq = x_uniqstr( lbls, count=count)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_widget_slider, strct, indx, INIT=init, REINIT=reinit, DRAG=drag

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x = x_widget_slider, strct, [indx], /INIT, /REINIT, /DRAG [v1.0]'
    return
  endif 

  if not keyword_set(INDX) then indx = 0 else strct.curslide = indx

  ;; INIT
  if keyword_set(INIT) then begin
      ;; Drop list
      slen = strlen(strct.droplist)
      gd = lindgen(strct.nslide) 
      strct.drop_id = widget_combobox(strct.base_id, $
                                      VALUE=strct.droplist[gd], $
                                      UVALUE=strct.uname+'_DROP', $
                                     font=strct.fonts.big_font)
      ;; Slider
      strct.slider_id = cw_fslider(strct.base_id, $
                                   MAX=strct.max[indx], $
                                   MIN=strct.min[indx], $
                                   xsize=strct.xsize[indx], $
                                   /DOUBLE, $
                                   DRAG=drag, $
                                   VALUE=strct.min[indx], $
                                   UVALUE=strct.uname)
      ;; Min/Max
      mnx_base = widget_base(strct.base_id, /row, /align_center)
      strct.min_id = cw_field(mnx_base, title='Min', $
                              value=strct.min[indx], /floating, $
                              /row, xsize=8, /return_events, $
                              uvalue=strct.uname+'_MIN',$
                             font=strct.fonts.small_font)
      strct.max_id = cw_field(mnx_base, title='Max', $
                              value=strct.max[indx], /floating, $
                              /row, xsize=8, /return_events, $
                              uvalue=strct.uname+'_MAX',$
                             font=strct.fonts.small_font)
  endif

  ;; REINIT
  if keyword_set(REINIT) then begin
      if not keyword_set(INDX) then indx = strct.curslide
      ;; Reset the values
      widget_control, strct.slider_id, set_value=[strct.min[indx]>$
                                                      strct.value[indx]<$
                                                      strct.max[indx], $
                                                      strct.min[indx],$
                                                      strct.max[indx]]
      widget_control, strct.min_id, set_value=strct.min[indx]
      widget_control, strct.max_id, set_value=strct.max[indx]
  endif

  return
end
