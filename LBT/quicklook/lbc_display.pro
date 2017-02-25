;;procedure that gets a lbtc mosaic in inputs and display a final fits

;imagename: fits file as it comes from the telescope
;mosaic:    if set to a variable, in output a single array full mosaic
;           (not included fractional rotation)

;written by MF Sept 2009



PRO lbc_display, imagename,BLUE=BLUE,RED=RED,BIASRED=BIASRED $
                 ,BIASBLUE=BIASBLUE


IF NOT KEYWORD_SET(LOWER) THEN LOWER = 2.0
IF NOT KEYWORD_SET(UPPER) THEN UPPER = 8.0

common atv_color, r_vector, g_vector, b_vector
common atv_state, state, tot_state, nimg, flg_img

path='/home/LBTB/'
pid = call_external(!dlm_path+'/libidl.' +idlutils_so_ext() $
                  , 'getpid', /cdecl)
IF KEYWORD_SET(BLUE) THEN BEGIN
    color='BLUE'
    print,'rm -f  '+ path + COLOR + '_PID_*'
    spawn,'rm -f  '+ path + COLOR + '_PID_*'
    spawn,'touch ' + path + 'BLUE_PID_' + strcompress(string(pid),/rem)
    IF KEYWORD_SET(BIASBLUE) THEN BIASFRAME=BIASBLUE $
    ELSE BIASFRAME='/TestRepo/20090920_LBTB/lbcb.20090920.031527.fits'
ENDIF ELSE IF KEYWORD_SET(RED) THEN BEGIN
    color='RED'
    print,'rm -f  '+ path + COLOR + '_PID_*'
    spawn,'rm -f  '+ path + COLOR + '_PID_*'
    spawn,'touch ' + path  +'RED_PID_' + strcompress(string(pid),/rem)
    IF KEYWORD_SET(BIASRED) THEN BIASFRAME=BIASRED $
    ELSE BIASFRAME='/TestRepo/20090920_LBTB/lbcr.20090920.031400.fits'
ENDIF

;size chips 2304*4608
;prescan 0:50, overscan 2099:2304
;data [51:2098,1:4608]

;;check if science frame
science=1
hed0=headfits(imagename,exten=0)
filter=sxpar(hed0,'FILTER')
imgname=sxpar(hed0,'FILENAME')
title=strcompress(imgname,/rem) + '---' + strcompress(filter,/rem)

numext=sxpar(hed0,'NEXTEND')
IF (numext EQ 1) THEN science=0 ELSE science=1

;no science
IF (science EQ 0) THEN BEGIN 
    mosaic=xmrdfits(imagename,1,/fscale,/silent)
    IF KEYWORD_SET(BIASFRAME) THEN $
       bias = xmrdfits(biasframe, 1, /fscale, /silent) $
    ELSE bias = 0.0d
    mosaic=mosaic-bias
    sky, mosaic, mean_sky, mean_sig, /SILENT
ENDIF ELSE BEGIN
    fit=make_array(2048,4608,/double)
    pres=make_array(50,4608,/double)
    overs=make_array(206,4608,/double)
    hea=make_array(1D4,/string)
    imgstruc={header:hea,data:fit,prescan:pres,overscan:overs}  
    imgstruc=replicate(imgstruc,5)
;;iterate over levels
    for i = 0L, 4L do begin
;open file and load in structure
        fits=xmrdfits(imagename,i,hea,/silent,/fscale)
        imgstruc[i].header=hea
        IF (i NE 0) THEN BEGIN
           IF KEYWORD_SET(BIASFRAME) THEN $
              bias = xmrdfits(biasframe, i, /silent, /fscale) $
           ELSE bias = 0.0*fits
            imgstruc[i].data=fits[50:2097,*]-bias[50:2097,*]
            imgstruc[i].prescan=fits[0:49,*]
            imgstruc[i].overscan=fits[2098:2303,*]
        ENDIF
    endfor
;;create fits mosaic (6178*6673)
    mosaic=make_array(6178,6673)
    
;;fill in third chip
    skymode=fltarr(4)
    skysig=fltarr(4)
    FOR chip=1L,4L DO BEGIN
        sky, imgstruc[chip].data, skymode1, skysig1, /SILENT
        skymode[chip-1] = skymode1
        skysig[chip-1]  = skysig1
    ENDFOR
    mean_sky = total(skymode)/4.0D
    mean_sig = total(skysig)/4.0D

    mosaic[0:2047,0:4607]=imgstruc[3].data
;fill in second chip, leaving 18 pixx for chip gap
    mosaic[2065:4112,0:4607]=imgstruc[2].data
;fill in first chip, leaving 18 pixx for chip gap
    mosaic[4130:6177,0:4607]=imgstruc[1].data
;rotate and fill in 4th chip,  leaving 18 pixy for chip gap
    chip4=rotate(imgstruc[4].data,1)
    mosaic[770:5377,4625:6672]=chip4
ENDELSE

;;display image
min = mean_sky - lower*mean_sig
max = mean_sky + upper*mean_sig
atv, mosaic,min=min,max=max,header=hed0

;; Set title
widget_control, state.base_id, tlb_set_title = title
;state.title_extras=title
;atv_settitle

;; Recenter image
x_cen = 3089
y_cen = 3336
state.centerpix = [x_cen, y_cen]
;; Invert colormap
IF state.invert_colormap EQ 0 THEN BEGIN
    state.invert_colormap = abs(state.invert_colormap - 1)
    r_vector = reverse(r_vector)
    g_vector = reverse(g_vector)
    b_vector = reverse(b_vector)
ENDIF
;; Add something here to modify the zoom (Subaru mask code)

;; Toggle brightness and contrast
IF KEYWORD_SET(BRIGHTNESS) THEN state.brightness = brightness $
ELSE state.brightness = 0.30
IF KEYWORD_SET(CONTRAST) THEN state.contrast = contrast $
ELSE state.contrast = 0.15
atv_stretchct, state.brightness, state.contrast


;; Resize ATV window
resize=1
IF KEYWORD_SET(RESIZE) THEN BEGIN
widget_control, state.base_id, tlb_get_size=tmp_event
tmp_event = [900, 1200]
window = (state.base_min_size > tmp_event)
newbase = window - state.base_pad
newxsize = (tmp_event[0] - state.base_pad[0]) > $
  (state.base_min_size[0] - state.base_pad[0]) 
newysize = (tmp_event[1] - state.base_pad[1]) > $
  (state.base_min_size[1] - state.base_pad[1])
widget_control, state.draw_widget_id, $
  scr_xsize = newxsize, scr_ysize = newysize
widget_control, state.colorbar_widget_id, $
  scr_xsize = newxsize, scr_ysize = state.colorbar_height
state.draw_window_size = [newxsize, newysize]
atv_colorbar
widget_control, state.base_id, /clear_events
atv_refresh
ENDIF

nzoom=3L
FOR j = 0L, nzoom-1L DO atv_zoom, 'out'
;;if(ndzoom ge 1)then FOR j = 0L, ndzoom-1L DO atv_zoom, 'out'

atv_displayall
atv_refresh
;;atv, mosaic,min=min,max=max,/block
while 1 do begin & wait, 0.1 & void = widget_event(/NOWAIT) & endwhile


spawn,'kill -9 ' + strcompress(string(pid),/rem)
spawn,'rm -f ' + path +  COLOR + '_PID_*'

    
END
