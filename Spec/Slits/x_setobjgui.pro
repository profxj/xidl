;+
; NAME:
;       x_setobjgui
;   Version 1.1
; 
; PURPOSE: 
;       Interactive display of slit mask and slit positions and object
;       traces to allow verificaiton of the identificaitons and
;       traces.  Largely based on ATV
;
; CATEGORY: 
;       Image display.
;
; CALLING SEQUENCE:
;       x_setobjgui img, slitstr, objstr, [,min = min_value] 
;           [,max=max_value] 
;           [,/linear] [,/log] [,/histeq] [,/block]
;           [,/align] [,/stretch] [,header = header]
;
; REQUIRED INPUTS:
;   img  -- FITS file
;   slitstr -- Slit structure
;   objstr  -- Object sturcture
;
; OPTIONAL INPUTS:
;       array_name: a 2-D data array to display
;          OR
;       fits_file:  a fits file name, enclosed in single quotes
;
; KEYWORDS:
;       min:        minimum data value to be mapped to the color table
;       max:        maximum data value to be mapped to the color table
;       linear:     use linear stretch
;       log:        use log stretch 
;       histeq:     use histogram equalization
;       block:      block IDL command line until ATV terminates
;       align:      align image with previously displayed image
;       stretch:    keep same min and max as previous image
;       header:     FITS image header (string array) for use with data array
;       
; OUTPUTS:
;       None.  
; 
; COMMON BLOCKS:
;       x_setobjgui_state:  contains variables describing the display state
;       x_setobjgui_images: contains the internal copies of the display image
;       x_setobjgui_color:  contains colormap vectors
;       x_setobjgui_pdata:  contains plot and text annotation information
;
; RESTRICTIONS:
;       Requires IDL version 5.1 or greater.
;       Requires Craig Markwardt's cmps_form.pro routine.
;       Requires the GSFC IDL astronomy user's library routines.
;       Some features may not work under all operating systems.
;
; SIDE EFFECTS:
;       Modifies the color table.
;
; EXAMPLE:
;       To start x_setobjgui running, just enter the command 'x_setobjgui' at the idl
;       prompt, either with or without an array name or fits file name 
;       as an input.  Only one x_setobjgui window will be created at a time,
;       so if one already exists and another image is passed to x_setobjgui
;       from the idl command line, the new image will be displayed in 
;       the pre-existing x_setobjgui window.
;
; MODIFICATION HISTORY:
;       Written by Aaron J. Barth, with contributions by 
;       Douglas Finkbeiner, Michael Liu, David Schlegel, and
;       Wesley Colley.  First released 17 December 1998.
;
;       This version is 1.3, last modified 28 November 2000.
;       For the most current version, revision history, instructions,
;       list of known bugs, and further information, go to:
;              http://cfa-www.harvard.edu/~abarth/x_setobjgui/x_setobjgui.html
;       Revised by JXP Aug-2001
;
;-
;----------------------------------------------------------------------
;        x_setobjgui startup, initialization, and shutdown routines
;----------------------------------------------------------------------

pro x_setobjgui_initcommon

; Routine to initialize the x_setobjgui common blocks.  Use common blocks so
; that other IDL programs can access the x_setobjgui internal data easily.

common x_setobjgui_state, state, flg_img, sobj_slit, sobj_obj
common x_setobjgui_color, r_vector, g_vector, b_vector
common x_setobjgui_pdata, nplot, maxplot, plot_ptr
common x_setobjgui_images, $
  main_image, $
  display_image, $
  scaled_image,  $
  wave_image, $
  sig_image, $
  pan_image


state = {                   $
          version: '1.0', $              ; version # of this release
          head_ptr: ptr_new(), $         ; pointer to image header
          astr_ptr: ptr_new(), $         ; pointer to astrometry info structure
          wcstype: 'none', $             ; coord info type (none/angle/lambda)
          equinox: 'J2000', $            ; equinox of coord system
          display_coord_sys: 'RA--', $   ; coord system displayed
          display_equinox: 'J2000', $    ; equinox of displayed coords
          display_base60: 1B, $          ; Display RA,dec in base 60?
          imagename: '', $               ; image file name
          title_extras: '', $            ; extras for image title
          bitdepth: 8, $                 ; 8 or 24 bit color mode?
          screen_ysize: 1000, $          ; vertical size of screen
          base_id: 0L, $                 ; id of top-level base
          base_min_size: [800L, 600L], $ ; min size for top-level base
          draw_base_id: 0L, $            ; id of base holding draw window
          draw_window_id: 0L, $          ; window id of draw window
          draw_widget_id: 0L, $          ; widget id of draw widget
          track_window_id: 0L, $         ; widget id of tracking window
          pan_widget_id: 0L, $           ; widget id of pan window
          pan_window_id: 0L, $           ; window id of pan window
          active_window_id: 0L, $        ; user's active window outside x_setobjgui
          info_base_id: 0L, $            ; id of base holding info bars
          frame_id: 0L, $                ; id of Frame number
          location_bar_id: 0L, $         ; id of (x,y,value) label
          wavesig_bar_id: 0L, $          ; id of (sig, wave) label
          wavesig_flg: 0, $              ; flg for wave, sigma images
          wcs_bar_id: 0L, $              ; id of WCS label widget
          min_text_id: 0L,  $            ; id of min= widget
          max_text_id: 0L, $             ; id of max= widget
          menu_ids: lonarr(35), $        ; list of top menu items
          colorbar_base_id: 0L, $        ; id of colorbar base widget
          colorbar_widget_id: 0L, $      ; widget id of colorbar draw widget
          colorbar_window_id: 0L, $      ; window id of colorbar
          colorbar_height: 6L, $         ; height of colorbar in pixels
          ncolors: 0B, $                 ; image colors (!d.table_size - 9)
          box_color: 2, $                ; color for pan box and zoom x
          brightness: 0.5, $             ; initial brightness setting
          contrast: 0.5, $               ; initial contrast setting
          keyboard_text_id: 0L, $        ; id of keyboard input widget
          image_min: 0.0, $              ; min(main_image)
          image_max: 0.0, $              ; max(main_image)
          min_value: 0.0, $              ; min data value mapped to colors
          max_value: 0.0, $              ; max data value mapped to colors
          flg_flip: 0, $                 ; Flip 0=none, 1=x, 2=y, 3=both
          flipcoord: [0L,0L], $          ; Flipped coordinates 
          draw_window_size: [712L, 712L], $    ; size of main draw window
          track_window_size: 121L, $     ; size of tracking window
          pan_window_size: 121L, $       ; size of pan window
          pan_scale: 0.0, $              ; magnification of pan image
          image_size: [0L,0L], $         ; size of main_image
          invert_colormap: 0L, $         ; 0=normal, 1=inverted
          coord: [0L, 0L],  $            ; cursor position in image coords
          scaling: 0L, $                 ; 0=linear, 1=log, 2=histeq
          offset: [0L, 0L], $            ; offset to viewport coords
          base_pad: [0L, 0L], $          ; padding around draw base
          zoom_level: 0L, $              ; integer zoom level, 0=normal
          zoom_factor: 1.0, $            ; magnification factor = 2^zoom_level
          centerpix: [0L, 0L], $         ; pixel at center of viewport
          cstretch: 0B, $                ; flag = 1 while stretching colors
          pan_offset: [0L, 0L], $        ; image offset in pan window
          frame: 1L, $                   ; put frame around ps output?
          framethick: 6, $               ; thickness of frame
          lineplot_widget_id: 0L, $      ; id of lineplot widget
          lineplot_window_id: 0L, $      ; id of lineplot window
          lineplot_base_id: 0L, $        ; id of lineplot top-level base
          lineplot_size: [600L, 450L], $ ; size of lineplot window
          lineplot_pad: [0L, 0L], $      ; padding around lineplot window
          cursorpos: lonarr(2), $        ; cursor x,y for photometry & stats
          centerpos: fltarr(2), $        ; centered x,y for photometry
          cursorpos_id: 0L, $            ; id of cursorpos widget
          centerpos_id: 0L, $            ; id of centerpos widget
          centerbox_id: 0L, $            ; id of centeringboxsize widget
          radius_id: 0L, $               ; id of radius widget
          innersky_id: 0L, $             ; id of inner sky widget
          outersky_id: 0L, $             ; id of outer sky widget
          magunits: 0, $                 ; 0=counts, 1=magnitudes
          skytype: 0, $                  ; 0=idlphot,1=median,2=no sky subtract
          photzpt: 25.0, $               ; magnitude zeropoint
          skyresult_id: 0L, $            ; id of sky widget
          photresult_id: 0L, $           ; id of photometry result widget
          fwhm_id: 0L, $                 ; id of fwhm widget
          radplot_widget_id: 0L, $       ; id of radial profile widget
          radplot_window_id: 0L, $       ; id of radial profile window
          photzoom_window_id: 0L, $      ; id of photometry zoom window
          photzoom_size: 190L, $         ; size in pixels of photzoom window
          showradplot_id: 0L, $          ; id of button to show/hide radplot
          photwarning_id: 0L, $          ; id of photometry warning widget
          photwarning: ' ', $            ; photometry warning text
          centerboxsize: 5L, $           ; centering box size
          r: 5L, $                       ; aperture photometry radius
          innersky: 10L, $               ; inner sky radius
          outersky: 20L, $               ; outer sky radius
          headinfo_base_id: 0L, $        ; headinfo base widget id
          stats_base_id: 0L, $           ; base widget for image stats
          statboxsize: 11L, $            ; box size for computing statistics
          statbox_id: 0L, $              ; widget id for stat box size 
          stat_npix_id: 0L, $            ; widget id for # pixels in stats box
          statxcenter_id: 0L, $          ; widget id for stat box x center
          statycenter_id: 0L, $          ; widget id for stat box y center
          statbox_min_id: 0L, $          ; widget id for stat min box
          statbox_max_id: 0L, $          ; widget id for stat max box
          statbox_mean_id: 0L, $         ; widget id for stat mean box
          statbox_median_id: 0L, $       ; widget id for stat median box
          statbox_stdev_id: 0L, $        ; widget id for stat stdev box
          statzoom_size: 300, $          ; size of statzoom window
          statzoom_widget_id: 0L, $      ; widget id for stat zoom window
          statzoom_window_id: 0L, $      ; window id for stat zoom window
          showstatzoom_id: 0L, $         ; widget id for show/hide button
          pan_pixmap: 0L, $              ; window id of pan pixmap
          default_autoscale: 1, $        ; autoscale images by default?
          current_dir: '', $             ; current readfits directory
          graphicsdevice: '', $          ; screen device
          nobj: 0, $                     ;;;;; NEW STUFF ;;;;;;;;
          obj: lonarr(500), $                     
          ydelt: 0., $                     
          sz: lonarr(2), $                     
          newrefresh: 0, $               ; refresh since last blink?
          blinks: 0B $                   ; remembers which images are blinked
        }


nplot = 0
maxplot = 5000
plot_ptr = ptrarr(maxplot+1)  ; The 0th element isn't used.


end

;---------------------------------------------------------------------

pro x_setobjgui_startup

; This routine initializes the x_setobjgui internal variables, and creates and
; realizes the window widgets.  It is only called by the x_setobjgui main
; program level, when there is no previously existing x_setobjgui window.

common x_setobjgui_state
common x_setobjgui_color

; Read in a color table to initialize !d.table_size
; As a bare minimum, we need the 8 basic colors used by ATV_ICOLOR(),
; plus 2 more for a color map.

loadct, 0, /silent
if (!d.table_size LT 12) then begin
    message, 'Too few colors available for color table'
    x_setobjgui_shutdown
endif

; Initialize the common blocks
x_setobjgui_initcommon

state.ncolors = !d.table_size - 9
if (!d.n_colors LE 256) then begin
    state.bitdepth = 24
;    state.bitdepth = 8
endif else begin
    state.bitdepth = 24
    device, decomposed=0
endelse

state.graphicsdevice = !d.name

state.screen_ysize = (get_screen_size())[1]

; Get the current window id
x_setobjgui_getwindow


; Define the widgets.  For the widgets that need to be modified later
; on, save their widget ids in state variables

base = widget_base(title = 'x_setobjgui', $
                   /column, /base_align_right, $
                   app_mbar = top_menu, $
                   uvalue = 'x_setobjgui_base', $
                   /tlb_size_events)
state.base_id = base

tmp_struct = {cw_pdmenu_s, flags:0, name:''}

top_menu_desc = [ $
                  {cw_pdmenu_s, 1, 'File'}, $ ; file menu
                  {cw_pdmenu_s, 0, 'ReadFits'}, $
                  {cw_pdmenu_s, 0, 'WritePS'},  $
                  {cw_pdmenu_s, 0, 'WriteTiff'}, $
                  {cw_pdmenu_s, 2, 'Quit'}, $
                  {cw_pdmenu_s, 1, 'ColorMap'}, $ ; color menu
                  {cw_pdmenu_s, 0, 'Grayscale'}, $
                  {cw_pdmenu_s, 0, 'Blue-White'}, $
                  {cw_pdmenu_s, 0, 'Red-Orange'}, $
                  {cw_pdmenu_s, 0, 'Green-White'}, $
                  {cw_pdmenu_s, 0, 'Rainbow'}, $
                  {cw_pdmenu_s, 0, 'BGRY'}, $
                  {cw_pdmenu_s, 0, 'Stern Special'}, $
                  {cw_pdmenu_s, 2, 'ATV Special'}, $
                  {cw_pdmenu_s, 1, 'Scaling'}, $ ; scaling menu
                  {cw_pdmenu_s, 0, 'Linear'}, $
                  {cw_pdmenu_s, 0, 'Log'}, $
                  {cw_pdmenu_s, 2, 'HistEq'}, $
                  {cw_pdmenu_s, 1, 'Labels'}, $ ; labels menu
                  {cw_pdmenu_s, 0, 'TextLabel'}, $
                  {cw_pdmenu_s, 0, 'Contour'}, $
                  {cw_pdmenu_s, 0, 'Compass'}, $
                  {cw_pdmenu_s, 0, 'ScaleBar'}, $
                  {cw_pdmenu_s, 0, 'EraseLast'}, $
                  {cw_pdmenu_s, 2, 'EraseAll'}, $
                  {cw_pdmenu_s, 1, 'Flip'}, $
                  {cw_pdmenu_s, 0, 'Flip-0'}, $
                  {cw_pdmenu_s, 0, 'Flip-x'}, $
                  {cw_pdmenu_s, 0, 'Flip-y'}, $
                  {cw_pdmenu_s, 2, 'Flip-xy'}, $
                  {cw_pdmenu_s, 1, 'ImageInfo'}, $    ;info menu
                  {cw_pdmenu_s, 0, 'Photometry'}, $
                  {cw_pdmenu_s, 0, 'Statistics'}, $
                  {cw_pdmenu_s, 0, 'ImageHeader'}, $
                  {cw_pdmenu_s, 0, '--------------'}, $
                  {cw_pdmenu_s, 0, 'RA,dec (J2000)'}, $
                  {cw_pdmenu_s, 0, 'RA,dec (B1950)'}, $
                  {cw_pdmenu_s, 0, '--------------'}, $
                  {cw_pdmenu_s, 0, 'RA,dec (J2000) deg'}, $
                  {cw_pdmenu_s, 0, 'Galactic'}, $
                  {cw_pdmenu_s, 0, 'Ecliptic (J2000)'}, $
                  {cw_pdmenu_s, 2, 'Native'}, $
                  {cw_pdmenu_s, 1, 'Help'}, $         ; help menu
                  {cw_pdmenu_s, 2, 'ATV Help'} $
                ]

top_menu = cw_pdmenu(top_menu, top_menu_desc, $
                     ids = state.menu_ids, $
                     /mbar, $
                     /help, $
                     /return_name, $
                     uvalue = 'top_menu')

track_base =    widget_base(base, /row)
state.info_base_id = widget_base(track_base, /column, /base_align_right)
buttonbar_base = widget_base(base, column=2, /base_align_center)

state.draw_base_id = widget_base(base, $
                                 /column, /base_align_left, $
                                 uvalue = 'draw_base', $
                                 frame = 2, /tracking_events)

state.colorbar_base_id = widget_base(base, $
                                     uvalue = 'cqolorbar_base', $
                                     /column, /base_align_left, $
                                     frame = 2)

minmax_fram_base = widget_base(state.info_base_id, /row, /base_align_right)


;state.frame_id = cw_bgroup(minmax_fram_base, ['0','1','2','3'], $
;                           column=4, uvalue='frame_bg', /exclusive,$
;                          /no_release, /frame, LABEL_TOP='Frame', $
;                          set_value=0) 

min_base = widget_base(minmax_fram_base, /column, /base_align_right)

state.min_text_id = cw_field(min_base, $
                             uvalue = 'min_text', $
                             /floating,  $
                             title = 'Min=', $
                             value = state.min_value,  $
                             /return_events, $
                             xsize = 12)

state.max_text_id = cw_field(min_base, $
                             uvalue = 'max_text', $
                             /floating,  $
                             title = 'Max=', $
                             value = state.max_value, $
                             /return_events, $
                             xsize = 12)

tmp_string = string(1000, 1000, 1.0d-10, $
                    format = '("(",i4,",",i4,") ",g14.7)' )

state.location_bar_id = widget_label (state.info_base_id, $
                                      value = tmp_string,  $
                                      uvalue = 'location_bar',  frame = 1)

tmp_string = string(1.0e-10, 1.0d-10, $
                    format = '(g12.5,3x,g14.7)' )
state.wavesig_bar_id = widget_label (state.info_base_id, $
                                      value = tmp_string,  $
                                      uvalue = 'wavesig_bar',  frame = 1)

tmp_string = string(12, 12, 12.001, -60, 60, 60.01, ' J2000', $
        format = '(i2,":",i2,":",f6.3,"  ",i3,":",i2,":",f5.2," ",a6)' )
    
;state.wcs_bar_id = widget_label (state.info_base_id, $
;                                 value = tmp_string,  $
;                                 uvalue = 'wcs_bar',  frame = 1)

state.pan_widget_id = widget_draw(track_base, $
                                  xsize = state.pan_window_size, $
                                  ysize = state.pan_window_size, $
                                  frame = 2, uvalue = 'pan_window', $
                                  /button_events, /motion_events)

track_window = widget_draw(track_base, $
                           xsize=state.track_window_size, $
                           ysize=state.track_window_size, $
                           frame=2, uvalue='track_window')

;modebase = widget_base(buttonbar_base, /row, /base_align_center)
;modelist = ['Color', 'Zoom', 'Blink', 'ImExam']
;mode_droplist_id = widget_droplist(modebase, $
;                                   frame = 1, $
;                                   title = 'MouseMode:', $
;                                   uvalue = 'mode', $
;                                   value = modelist)

button_base = widget_base(buttonbar_base, row=2, /base_align_right)

invert_button = widget_button(button_base, $
                              value = 'Invert', $
                              uvalue = 'invert')

restretch_button = widget_button(button_base, $
                             value = 'Restretch', $
                             uvalue = 'restretch_button')

autoscale_button = widget_button(button_base, $
                                 uvalue = 'autoscale_button', $
                                 value = 'AutoScale')

fullrange_button = widget_button(button_base, $
                                 uvalue = 'full_range', $
                                 value = 'FullRange')

state.keyboard_text_id = widget_text(button_base, $
                                     /all_events, $
                                     scr_xsize = 1, $
                                     scr_ysize = 1, $
                                     units = 0, $
                                     uvalue = 'keyboard_text', $
                                     value = '')

zoomin_button = widget_button(button_base, $
                              value = 'ZoomIn', $
                              uvalue = 'zoom_in')

zoomout_button = widget_button(button_base, $
                               value = 'ZoomOut', $
                               uvalue = 'zoom_out')

zoomone_button = widget_button(button_base, $
                               value = 'Zoom1', $
                               uvalue = 'zoom_one')

center_button = widget_button(button_base, $
                              value = 'Center', $
                              uvalue = 'center')

done_button = widget_button(button_base, $
                            value = 'Done', $
                            uvalue = 'done')

; Set widget y size for small screens
state.draw_window_size[1] = state.draw_window_size[1] < $
  (state.screen_ysize - 300)

state.draw_widget_id = widget_draw(state.draw_base_id, $
                                   uvalue = 'draw_window', $
                                   /motion_events,  /button_events, $
                                   scr_xsize = state.draw_window_size[0], $
                                   scr_ysize = state.draw_window_size[1]) 

state.colorbar_widget_id = widget_draw(state.colorbar_base_id, $
                                       uvalue = 'colorbar', $
                                       scr_xsize = state.draw_window_size[0], $
                                       scr_ysize = state.colorbar_height)

; Create the widgets on screen

widget_control, base, /realize
widget_control, state.pan_widget_id, draw_motion_events = 0

; get the window ids for the draw widgets

widget_control, track_window, get_value = tmp_value
state.track_window_id = tmp_value
widget_control, state.draw_widget_id, get_value = tmp_value
state.draw_window_id = tmp_value
widget_control, state.pan_widget_id, get_value = tmp_value
state.pan_window_id = tmp_value
widget_control, state.colorbar_widget_id, get_value = tmp_value
state.colorbar_window_id = tmp_value

; set the event handlers

widget_control, top_menu, event_pro = 'x_setobjgui_topmenu_event'
widget_control, state.draw_widget_id, event_pro = 'x_setobjgui_draw_color_event'
widget_control, state.draw_base_id, event_pro = 'x_setobjgui_draw_base_event'
widget_control, state.keyboard_text_id, event_pro = 'x_setobjgui_keyboard_event'
widget_control, state.pan_widget_id, event_pro = 'x_setobjgui_pan_event'

; Find window padding sizes needed for resizing routines.
; Add extra padding for menu bar, since this isn't included in 
; the geometry returned by widget_info.
; Also add extra padding for margin (frame) in draw base.

basegeom = widget_info(state.base_id, /geometry)
drawbasegeom = widget_info(state.draw_base_id, /geometry)

state.base_pad[0] = basegeom.xsize - drawbasegeom.xsize $
  + (2 * basegeom.margin)
state.base_pad[1] = basegeom.ysize - drawbasegeom.ysize + 30 $
  + (2 * basegeom.margin)

state.base_min_size = [state.base_pad[0] + state.base_min_size[0], $
                       state.base_pad[1] + 100]

; Initialize the vectors that hold the current color table.
; See the routine x_setobjgui_stretchct to see why we do it this way.

r_vector = bytarr(state.ncolors)
g_vector = bytarr(state.ncolors)
b_vector = bytarr(state.ncolors)

x_setobjgui_getct, 0
state.invert_colormap = 0

; Create a pixmap window to hold the pan image
window, /free, xsize=state.pan_window_size, ysize=state.pan_window_size, $
  /pixmap
state.pan_pixmap = !d.window
x_setobjgui_resetwindow

x_setobjgui_colorbar

end

;--------------------------------------------------------------------

pro x_setobjgui_colorbar

; Routine to tv the colorbar at the bottom of the x_setobjgui window

common x_setobjgui_state

x_setobjgui_setwindow, state.colorbar_window_id

xsize = (widget_info(state.colorbar_widget_id, /geometry)).xsize

b = congrid( findgen(state.ncolors), xsize) + 8
c = replicate(1, state.colorbar_height)
a = b # c

tv, a

x_setobjgui_resetwindow

end

;--------------------------------------------------------------------

pro x_setobjgui_shutdown, windowid

; routine to kill the x_setobjgui window(s) and clear variables to conserve
; memory when quitting x_setobjgui.  The windowid parameter is used when
; x_setobjgui_shutdown is called automatically by the xmanager, if x_setobjgui is
; killed by the window manager.

common x_setobjgui_images
common x_setobjgui_state
common x_setobjgui_color
common x_setobjgui_pdata

; Kill top-level base if it still exists
if (xregistered ('x_setobjgui')) then widget_control, state.base_id, /destroy

; Destroy all pointers to plots and their heap variables
;if (nplot GT 0) then begin
;    x_setobjguierase, /norefresh
;endif

if (size(state, /tname) EQ 'STRUCT') then begin
;    if (!d.name EQ state.graphicsdevice) then wdelete, state.pan_pixmap
    if (ptr_valid(state.head_ptr)) then ptr_free, state.head_ptr
    if (ptr_valid(state.astr_ptr)) then ptr_free, state.astr_ptr
endif

delvarx, plot_ptr
delvarx, main_image
delvarx, display_image
delvarx, wave_image
delvarx, sig_image
delvarx, scaled_image
delvarx, pan_image
delvarx, r_vector
delvarx, g_vector
delvarx, b_vector
delvarx, state

return    
end

;--------------------------------------------------------------------
;                  main x_setobjgui event loops
;--------------------------------------------------------------------

pro x_setobjgui_topmenu_event, event

; Event handler for top menu

common x_setobjgui_state
common x_setobjgui_images

widget_control, event.id, get_uvalue = event_name

if (!d.name NE state.graphicsdevice and event_name NE 'Quit') then return
if (state.bitdepth EQ 24) then true = 1 else true = 0

; Need to get active window here in case mouse goes to menu from top
; of x_setobjgui window without entering the main base
x_setobjgui_getwindow

case event_name of
    
; File menu options:
    'ReadFits': begin
        x_setobjgui_readfits, newimage=newimage
        if (newimage EQ 1) then begin
            x_setobjgui_getstats
            x_setobjgui_settitle
            state.zoom_level =  0
            state.zoom_factor = 1.0
            if (state.default_autoscale EQ 1) then x_setobjgui_autoscale
            x_setobjgui_set_minmax
            x_setobjgui_displayall
        endif
    end
    'WritePS' : x_setobjgui_writeps
    'WriteTiff': x_setobjgui_writetiff
    'Quit':     x_setobjgui_shutdown
; ColorMap menu options:            
    'Grayscale': x_setobjgui_getct, 0
    'Blue-White': x_setobjgui_getct, 1
    'Red-Orange': x_setobjgui_getct, 3
    'BGRY': x_setobjgui_getct, 4
    'Rainbow': x_setobjgui_getct, 13
    'Stern Special': x_setobjgui_getct, 15
    'Green-White': x_setobjgui_getct, 8
    'ATV Special': x_setobjgui_makect, event_name
; Scaling options:
    'Linear': begin
        state.scaling = 0
        x_setobjgui_displayall
    end
    'Log': begin
        state.scaling = 1
        x_setobjgui_displayall
    end

    'HistEq': begin
        state.scaling = 2
        x_setobjgui_displayall
    end

; Label options:
    'TextLabel': x_setobjgui_textlabel
    'Contour': x_setobjgui_oplotcontour
    'Compass': x_setobjgui_setcompass
    'ScaleBar': x_setobjgui_setscalebar
    'EraseLast': x_setobjguierase, 1
    'EraseAll': x_setobjguierase
; Blink options:
    'SetBlink1': begin   
        x_setobjgui_setwindow, state.draw_window_id
        blink_image1 = tvrd(true = true) 
    end
    'SetBlink2': begin   
        x_setobjgui_setwindow, state.draw_window_id
        blink_image2 = tvrd(true = true)
    end
    'SetBlink3': begin   
        x_setobjgui_setwindow, state.draw_window_id
        blink_image3 = tvrd(true = true)
    end

; Flip options:
    'Flip-0': begin
        state.flg_flip = 0
        x_setobjgui_scaleimage
        x_setobjgui_makepan
        x_setobjgui_refresh
    end
    'Flip-x': begin
        state.flg_flip = 1
        x_setobjgui_scaleimage
        x_setobjgui_makepan
        x_setobjgui_refresh
    end
    'Flip-y': begin
        state.flg_flip = 2
        x_setobjgui_scaleimage
        x_setobjgui_makepan
        x_setobjgui_refresh
    end
    'Flip-xy': begin
        state.flg_flip = 3
        x_setobjgui_scaleimage
        x_setobjgui_makepan
        x_setobjgui_refresh
    end
; Info options:
    'Photometry': x_setobjgui_apphot
    'ImageHeader': x_setobjgui_headinfo
    'Statistics': x_setobjgui_showstats

; Help options:            
    'ATV Help': x_setobjgui_help
    
    else: print, 'Unknown event in file menu!'
endcase

; Need to test whether x_setobjgui is still alive, since the quit option
; might have been selected.        
if (xregistered('x_setobjgui', /noshow)) then x_setobjgui_resetwindow

end

;--------------------------------------------------------------------

pro x_setobjgui_draw_color_event, event

; Event handler for color mode

common x_setobjgui_state
common x_setobjgui_images

if (!d.name NE state.graphicsdevice) then return

case event.type of
    0: begin           ; button press
        if (event.press EQ 1) then begin
            state.cstretch = 1
            x_setobjgui_stretchct, event.x, event.y, /getmouse
            x_setobjgui_colorbar
        endif else begin
            x_setobjgui_zoom, 'none', /recenter
        endelse
    end
    1: begin
        state.cstretch = 0  ; button release
        if (state.bitdepth EQ 24) then x_setobjgui_refresh
        x_setobjgui_draw_motion_event, event
    end
    2: begin                ; motion event
        if (state.cstretch EQ 1) then begin
            x_setobjgui_stretchct, event.x, event.y, /getmouse 
            if (state.bitdepth EQ 24) then x_setobjgui_refresh, /fast
        endif else begin 
            x_setobjgui_draw_motion_event, event
        endelse
    end 
endcase

widget_control, state.keyboard_text_id, /sensitive, /input_focus

end

;--------------------------------------------------------------------

pro x_setobjgui_draw_zoom_event, event

; Event handler for zoom mode

common x_setobjgui_state
 
if (!d.name NE state.graphicsdevice) then return

if (event.type EQ 0) then begin 
    case event.press of
        1: x_setobjgui_zoom, 'in', /recenter
        2: x_setobjgui_zoom, 'none', /recenter
        4: x_setobjgui_zoom, 'out', /recenter
    endcase
endif

if (event.type EQ 2) then x_setobjgui_draw_motion_event, event

widget_control, state.keyboard_text_id, /sensitive, /input_focus

end

;---------------------------------------------------------------------

pro x_setobjgui_draw_blink_event, event

; Event handler for blink mode

common x_setobjgui_state
common x_setobjgui_images

if (!d.name NE state.graphicsdevice) then return
if (state.bitdepth EQ 24) then true = 1 else true = 0

case event.type of
    0: begin                    ; button press
        x_setobjgui_setwindow, state.draw_window_id
                                ; define the unblink image if needed
        if ((state.newrefresh EQ 1) AND (state.blinks EQ 0)) then begin
            unblink_image = tvrd(true = true)
            state.newrefresh = 0
        endif
        
        case event.press of
            1: if n_elements(blink_image1) GT 1 then $
              tv, blink_image1, true = true
            2: if n_elements(blink_image2) GT 1 then $
              tv, blink_image2, true = true
            4: if n_elements(blink_image3) GT 1 then $
              tv, blink_image3, true = true  
            else: event.press = 0 ; in case of errors
        endcase
        state.blinks = (state.blinks + event.press) < 7
    end
    
    1: begin                    ; button release
        if (n_elements(unblink_image) EQ 0) then return ; just in case
        x_setobjgui_setwindow, state.draw_window_id
        state.blinks = (state.blinks - event.release) > 0
        case state.blinks of
            0: tv, unblink_image, true = true
            1: if n_elements(blink_image1) GT 1 then $
              tv, blink_image1, true = true else $
              tv, unblink_image, true = true
            2: if n_elements(blink_image2) GT 1 then $
              tv, blink_image2, true = true else $
              tv, unblink_image, true = true
            3: if n_elements(blink_image1) GT 1 then begin
                tv, blink_image1, true = true
            endif else if n_elements(blink_image2) GT 1 then begin
                tv, blink_image2, true = true
            endif else begin
                tv, unblink_image, true = true
            endelse
            4: if n_elements(blink_image3) GT 1 then $
              tv, blink_image3, true = true $
            else tv, unblink_image, true = true
            5: if n_elements(blink_image1) GT 1 then begin
                tv, blink_image1, true = true 
            endif else if n_elements(blink_image3) GT 1 then begin
                tv, blink_image3, true = true
            endif else begin
                tv, unblink_image, true = true
            endelse 
            6: if n_elements(blink_image2) GT 1 then begin
                tv, blink_image2, true = true
            endif else if n_elements(blink_image4) GT 1 then begin
                tv, blink_image4, true = true
            endif else begin
                tv, unblink_image, true = true
            endelse
            else: begin         ; check for errors
                state.blinks = 0
                tv, unblink_image, true = true
            end
        endcase
    end
    2: x_setobjgui_draw_motion_event, event ; motion event
endcase

widget_control, state.keyboard_text_id, /sensitive, /input_focus
x_setobjgui_resetwindow

end

;--------------------------------------------------------------------

pro x_setobjgui_draw_motion_event, event

; Event handler for motion events in draw window

common x_setobjgui_state

if (!d.name NE state.graphicsdevice) then return

tmp_event = [event.x, event.y]            
state.coord = $
  round( (0.5 >  ((tmp_event / state.zoom_factor) + state.offset) $
          < (state.image_size - 0.5) ) - 0.5)

state.flipcoord = state.coord

; Allow for Flipping (JXP)
case state.flg_flip of
    0 : 
    1 : state.flipcoord[0] = state.image_size[0]-1-state.coord[0] 
    2 : state.flipcoord[1] = state.image_size[1]-1-state.coord[1] 
    3 : begin
        state.flipcoord[0] = state.image_size[0]-1-state.coord[0] 
        state.flipcoord[1] = state.image_size[1]-1-state.coord[1] 
    end
endcase

x_setobjgui_gettrack

end

;--------------------------------------------------------------------

pro x_setobjgui_draw_base_event, event

; event handler for exit events of main draw base.  There's no need to
; define enter events, since as soon as the pointer enters the draw
; window the motion event will make the text widget sensitive again.
; Enter/exit events are often generated incorrectly, anyway.

common x_setobjgui_state

if (event.enter EQ 0) then begin
    widget_control, state.keyboard_text_id, sensitive = 0
endif

end

;----------------------------------------------------------------------

pro x_setobjgui_keyboard_event, event

; Event procedure for keyboard input when the cursor is in the 
; main draw window.

common x_setobjgui_state

eventchar = string(event.ch)

if (!d.name NE state.graphicsdevice and eventchar NE 'q') then return

case eventchar of
    '1': x_setobjgui_move_cursor, eventchar
    '2': x_setobjgui_move_cursor, eventchar
    '3': x_setobjgui_move_cursor, eventchar
    '4': x_setobjgui_move_cursor, eventchar
    '6': x_setobjgui_move_cursor, eventchar
    '7': x_setobjgui_move_cursor, eventchar
    '8': x_setobjgui_move_cursor, eventchar
    '9': x_setobjgui_move_cursor, eventchar
    'r': x_setobjgui_rowplot
    'c': x_setobjgui_colplot
;    's': x_setobjgui_surfplot
;    't': x_setobjgui_contourplot
;    'p': x_setobjgui_apphot
    'i': x_setobjgui_showstats
    'q': x_setobjgui_shutdown
;  Big pan
    '{': x_setobjgui_bigpan, 'up'
    '}': x_setobjgui_bigpan, 'down'
    '[': x_setobjgui_bigpan, 'left'
    ']': x_setobjgui_bigpan, 'right'
;  Modified 6-23-01 JXP
    'z': x_setobjgui_zoom, 'in', /recenter
    'Z': x_setobjgui_zoom, 'out', /recenter
;  Object Traces
    'm': x_setobjgui_moveobj  ; Move existing
    'd': x_setobjgui_delobj  ; Delete existing
    'n': x_setobjgui_addobj  ; Create serendip
;  SLITS
    't': x_setobjgui_sngshift
    'b': x_setobjgui_sngshift, /bottom
;  Query slit
    '?': x_setobjgui_query
    else:  ;any other key press does nothing
endcase

if (xregistered('x_setobjgui', /noshow)) then $
  widget_control, state.keyboard_text_id, /clear_events

end

;--------------------------------------------------------------------

pro x_setobjgui_pan_event, event

; event procedure for moving the box around in the pan window

common x_setobjgui_state

if (!d.name NE state.graphicsdevice) then return

case event.type of
    0: begin                     ; button press
        widget_control, state.pan_widget_id, draw_motion_events = 1
        x_setobjgui_pantrack, event
    end
    1: begin                     ; button release
        widget_control, state.pan_widget_id, draw_motion_events = 0
        widget_control, state.pan_widget_id, /clear_events
        x_setobjgui_pantrack, event
        x_setobjgui_refresh
    end
    2: begin
        x_setobjgui_pantrack, event     ; motion event
        widget_control, state.pan_widget_id, /clear_events
    end
    else:
endcase

end

;--------------------------------------------------------------------

pro x_setobjgui_event, event

; Main event loop for x_setobjgui top-level base, and for all the buttons.

common x_setobjgui_state
common x_setobjgui_images
common x_setobjgui_color

widget_control, event.id, get_uvalue = uvalue

if (!d.name NE state.graphicsdevice and uvalue NE 'done') then return

; Get currently active window
x_setobjgui_getwindow

case uvalue of

    'x_setobjgui_base': begin     
        c = where(tag_names(event) EQ 'ENTER', count)
        if (count EQ 0) then begin       ; resize event
            x_setobjgui_resize
            x_setobjgui_refresh
        endif
    end

    'mode': case event.index of
        0: widget_control, state.draw_widget_id, $
          event_pro = 'x_setobjgui_draw_color_event'
        1: widget_control, state.draw_widget_id, $
          event_pro = 'x_setobjgui_draw_zoom_event'
        2: widget_control, state.draw_widget_id, $
          event_pro = 'x_setobjgui_draw_blink_event'
        3: widget_control, state.draw_widget_id, $
          event_pro = 'x_setobjgui_draw_phot_event'
        else: print, 'Unknown mouse mode!'
    endcase

    'invert': begin                  ; invert the color table
        state.invert_colormap = abs(state.invert_colormap - 1)

        r_vector = reverse(r_vector)
        g_vector = reverse(g_vector)
        b_vector = reverse(b_vector)

        x_setobjgui_stretchct, state.brightness, state.contrast
        if (state.bitdepth EQ 24) then x_setobjgui_refresh
    end
    
    'restretch_button': x_setobjgui_restretch

    'min_text': begin     ; text entry in 'min = ' box
        x_setobjgui_get_minmax, uvalue, event.value
        x_setobjgui_displayall
    end

    'max_text': begin     ; text entry in 'max = ' box
        x_setobjgui_get_minmax, uvalue, event.value
        x_setobjgui_displayall
    end

    'autoscale_button': begin   ; autoscale the image
        x_setobjgui_autoscale
        x_setobjgui_displayall
    end

    'full_range': begin    ; display the full intensity range
        state.min_value = state.image_min
        state.max_value = state.image_max
        if state.min_value GE state.max_value then begin
            state.min_value = state.max_value - 1
            state.max_value = state.max_value + 1
        endif
        x_setobjgui_set_minmax
        x_setobjgui_displayall
    end
    
    'zoom_in':  x_setobjgui_zoom, 'in'         ; zoom buttons
    'zoom_out': x_setobjgui_zoom, 'out'
    'zoom_one': x_setobjgui_zoom, 'one'

    'center': begin   ; center image and preserve current zoom level
        state.centerpix = round(state.image_size / 2.)
        x_setobjgui_refresh
    end

;   Frame Number

    'frame_bg' : x_setobjgui_setframe, event.value

    'done':  x_setobjgui_shutdown

    else:  print, 'No match for uvalue....'  ; bad news if this happens

endcase
end

;----------------------------------------------------------------------

pro x_setobjgui_message, msg_txt, msgtype=msgtype, window=window

; Routine to display an error or warning message.  Message can be
; displayed either to the IDL command line or to a popup window,
; depending on whether /window is set.
; msgtype must be 'warning', 'error', or 'information'.

common x_setobjgui_state

if (n_elements(window) EQ 0) then window = 0

if (window EQ 1) then begin  ; print message to popup window
    case msgtype of
        'warning': t = dialog_message(msg_txt, dialog_parent = state.base_id)
        'error': t = $
          dialog_message(msg_txt,/error,dialog_parent=state.base_id)
        'information': t = $
          dialog_message(msg_txt,/information,dialog_parent=state.base_id)
        else: 
    endcase
endif else begin           ;  print message to IDL console
    message = strcompress(strupcase(msgtype) + ': ' + msg_txt)
    print, message
endelse

end

;-----------------------------------------------------------------------
;      main x_setobjgui routines for scaling, displaying, cursor tracking...
;-----------------------------------------------------------------------

pro x_setobjgui_displayall

; Call the routines to scale the image, make the pan image, and
; re-display everything.  Use this if the scaling changes (log/
; linear/ histeq), or if min or max are changed, or if a new image is
; passed to x_setobjgui.  If the display image has just been moved around or
; zoomed without a change in scaling, then just call x_setobjgui_refresh
; rather than this routine.

x_setobjgui_scaleimage
x_setobjgui_makepan
x_setobjgui_refresh

end

;---------------------------------------------------------------------

pro x_setobjgui_refresh, fast = fast

; Make the display image from the scaled_image, and redisplay the pan
; image and tracking image. 
; The /fast option skips the steps where the display_image is
; recalculated from the main_image.  The /fast option is used in 24
; bit color mode, when the color map has been stretched but everything
; else stays the same.

common x_setobjgui_state
common x_setobjgui_images

x_setobjgui_getwindow
if (not(keyword_set(fast))) then begin
    x_setobjgui_getoffset
    x_setobjgui_getdisplay
    x_setobjgui_displaymain
    x_setobjgui_plotall
    x_setobjgui_plotobj
    x_setobjgui_plotslit
endif else begin
    x_setobjgui_displaymain
endelse

; redisplay the pan image and plot the boundary box
x_setobjgui_setwindow, state.pan_pixmap
erase
tv, pan_image, state.pan_offset[0], state.pan_offset[1]
x_setobjgui_resetwindow

x_setobjgui_setwindow, state.pan_window_id
if (not(keyword_set(fast))) then erase
tv, pan_image, state.pan_offset[0], state.pan_offset[1]
x_setobjgui_resetwindow
x_setobjgui_drawbox, /norefresh

if (state.bitdepth EQ 24) then x_setobjgui_colorbar

; redisplay the tracking image
if (not(keyword_set(fast))) then x_setobjgui_gettrack

x_setobjgui_resetwindow

state.newrefresh = 1


end

;--------------------------------------------------------------------

pro x_setobjgui_getdisplay

; make the display image from the scaled image by applying the zoom
; factor and matching to the size of the draw window, and display the
; image.

common x_setobjgui_state
common x_setobjgui_images

widget_control, /hourglass   

display_image = bytarr(state.draw_window_size[0], state.draw_window_size[1])

view_min = round(state.centerpix - $
                  (0.5 * state.draw_window_size / state.zoom_factor))
view_max = round(view_min + state.draw_window_size / state.zoom_factor)

view_min = (0 > view_min < (state.image_size - 1)) 
view_max = (0 > view_max < (state.image_size - 1)) 

newsize = round( (view_max - view_min + 1) * state.zoom_factor) > 1
startpos = abs( round(state.offset * state.zoom_factor) < 0)

tmp_image = congrid(scaled_image[view_min[0]:view_max[0], $
                                 view_min[1]:view_max[1]], $
                    newsize[0], newsize[1])

xmax = newsize[0] < (state.draw_window_size[0] - startpos[0])
ymax = newsize[1] < (state.draw_window_size[1] - startpos[1])

; Deal with flipping

display_image[startpos[0], startpos[1]] = tmp_image[0:xmax-1, 0:ymax-1]
delvarx, tmp_image

end

;-----------------------------------------------------------------------

pro x_setobjgui_displaymain

; Display the main image and overplots

common x_setobjgui_state
common x_setobjgui_images

x_setobjgui_setwindow, state.draw_window_id
tv, display_image
x_setobjgui_resetwindow

end

;--------------------------------------------------------------------

pro x_setobjgui_getoffset
common x_setobjgui_state

; Routine to calculate the display offset for the current value of
; state.centerpix, which is the central pixel in the display window.

state.offset = $
  round( state.centerpix - $
         (0.5 * state.draw_window_size / state.zoom_factor) )

end

;----------------------------------------------------------------------


pro x_setobjgui_makepan

; Make the 'pan' image that shows a miniature version of the full image.

common x_setobjgui_state
common x_setobjgui_images

sizeratio = state.image_size[1] / state.image_size[0]

if (sizeratio GE 1) then begin
    state.pan_scale = float(state.pan_window_size) / float(state.image_size[1])
endif else begin
    state.pan_scale = float(state.pan_window_size) / float(state.image_size[0])
endelse

tmp_image = $
  scaled_image[0:state.image_size[0]-1, 0:state.image_size[1]-1]

pan_image = $
  congrid(tmp_image, round(state.pan_scale * state.image_size[0])>1, $
          round(state.pan_scale * state.image_size[1])>1 )

state.pan_offset[0] = round((state.pan_window_size - (size(pan_image))[1]) / 2)
state.pan_offset[1] = round((state.pan_window_size - (size(pan_image))[2]) / 2)

end

;----------------------------------------------------------------------


pro x_setobjgui_move_cursor, direction

; Use keypad arrow keys to step cursor one pixel at a time.
; Get the new track image, and update the cursor position.

common x_setobjgui_state

i = 1L

case direction of
    '2': state.coord[1] = max([state.coord[1] - i, 0])
    '4': state.coord[0] = max([state.coord[0] - i, 0])
    '8': state.coord[1] = min([state.coord[1] + i, state.image_size[1] - i])
    '6': state.coord[0] = min([state.coord[0] + i, state.image_size[0] - i])
    '7': begin
        state.coord[1] = min([state.coord[1] + i, state.image_size[1] - i])
        state.coord[0] = max([state.coord[0] - i, 0])
    end
    '9': begin
        state.coord[1] = min([state.coord[1] + i, state.image_size[1] - i])
        state.coord[0] = min([state.coord[0] + i, state.image_size[0] - i])
    end
    '3': begin
        state.coord[1] = max([state.coord[1] - i, 0])
        state.coord[0] = min([state.coord[0] + i, state.image_size[0] - i])
    end
    '1': begin
        state.coord[1] = max([state.coord[1] - i, 0])
        state.coord[0] = max([state.coord[0] - i, 0])
    end

endcase

newpos = (state.coord - state.offset + 0.5) * state.zoom_factor

x_setobjgui_setwindow,  state.draw_window_id
tvcrs, newpos[0], newpos[1], /device
x_setobjgui_resetwindow

x_setobjgui_gettrack

; Prevent the cursor move from causing a mouse event in the draw window
widget_control, state.draw_widget_id, /clear_events

x_setobjgui_resetwindow

end

;----------------------------------------------------------------------

pro x_setobjgui_set_minmax

; Updates the min and max text boxes with new values.

common x_setobjgui_state

widget_control, state.min_text_id, set_value = string(state.min_value)
widget_control, state.max_text_id, set_value = string(state.max_value)

end

;----------------------------------------------------------------------

pro x_setobjgui_get_minmax, uvalue, newvalue

; Change the min and max state variables when user inputs new numbers
; in the text boxes. 

common x_setobjgui_state

case uvalue of
    
    'min_text': begin
        if (newvalue LT state.max_value) then begin
            state.min_value = newvalue
        endif
    end

    'max_text': begin
        if (newvalue GT state.min_value) then begin
            state.max_value = newvalue
        endif
    end
        
endcase

x_setobjgui_set_minmax

end

;--------------------------------------------------------------------

pro x_setobjgui_zoom, zchange, recenter = recenter
common x_setobjgui_state

; Routine to do zoom in/out and recentering of image.  The /recenter
; option sets the new display center to the current cursor position.

case zchange of
    'in':    state.zoom_level = (state.zoom_level + 1) < 12
    'out':   state.zoom_level = (state.zoom_level - 1) > (-12) 
    'one':   state.zoom_level =  0
    'none':  ; no change to zoom level: recenter on current mouse position
    else:  print,  'problem in x_setobjgui_zoom!'
endcase

; JXP
if state.zoom_level LE 0 then state.zoom_factor = 1./(1.+abs(state.zoom_level)*0.5)
if state.zoom_level GT 0 then state.zoom_factor = 2.^state.zoom_level

if keyword_set( RECENTER ) then begin
    state.centerpix = state.coord
    x_setobjgui_getoffset
endif

x_setobjgui_refresh

if keyword_set( RECENTER ) then begin
    newpos = (state.coord - state.offset + 0.5) * state.zoom_factor
    x_setobjgui_setwindow,  state.draw_window_id
    tvcrs, newpos[0], newpos[1], /device 
    x_setobjgui_resetwindow
    x_setobjgui_gettrack
endif

x_setobjgui_resetwindow

end

; BIG PAN
;--------------------------------------------------------------------

pro x_setobjgui_bigpan, direct
common x_setobjgui_state

; Routine to do zoom in/out and recentering of image.  The /recenter
; option sets the new display center to the current cursor position.

  case direct of
      'down': offset = [0.,-0.5*state.draw_window_size[1]/state.zoom_factor]
      'up': offset = [0.,0.5*state.draw_window_size[1]/state.zoom_factor]
      'right': offset = [0.5*state.draw_window_size[0]/state.zoom_factor,0.]
      'left': offset = [-0.5*state.draw_window_size[0]/state.zoom_factor,0.]
      else:  print,  'problem in x_setobjgui_bigpan!'
  endcase

  state.centerpix = state.centerpix+offset
  state.coord = state.centerpix
  x_setobjgui_getoffset
  
  x_setobjgui_refresh
  
  return

end

;-----------------------------------------------------------------------

pro x_setobjgui_autoscale

; Routine to auto-scale the image.  

common x_setobjgui_state 
common x_setobjgui_images

widget_control, /hourglass

if (n_elements(main_image) LT 5.e5) then begin
    med = median(main_image)
    sig = stddev(main_image)
endif else begin   ; resample big images before taking median, to save memory
    boxsize = 10
    rx = state.image_size[0] mod boxsize
    ry = state.image_size[1] mod boxsize
    nx = state.image_size[0] - rx
    ny = state.image_size[1] - ry
    tmp_img = rebin(main_image[0: nx-1, 0: ny-1], $
                    nx/boxsize, ny/boxsize, /sample)
    med = median(tmp_img)
    sig = stddev(temporary(tmp_img))
endelse

state.max_value = (med + (2 * sig)) < state.image_max
state.min_value = (med - (2 * sig))  > state.image_min

if (finite(state.min_value) EQ 0) then state.min_value = state.image_min
if (finite(state.max_value) EQ 0) then state.max_value = state.image_max

if (state.min_value GE state.max_value) then begin
    state.min_value = state.min_value - 1
    state.max_value = state.max_value + 1
endif

x_setobjgui_set_minmax

end  

;--------------------------------------------------------------------

pro x_setobjgui_restretch

; Routine to restretch the min and max to preserve the display
; visually but use the full color map linearly.  Written by DF, and
; tweaked and debugged by AJB.  It doesn't always work exactly the way
; you expect (especially in log-scaling mode), but mostly it works fine.

common x_setobjgui_state

sx = state.brightness
sy = state.contrast

if state.scaling EQ 2 then return ; do nothing for hist-eq mode

IF state.scaling EQ 0 THEN BEGIN 
    sfac = (state.max_value-state.min_value)
    state.max_value = sfac*(sx+sy)+state.min_value
    state.min_value = sfac*(sx-sy)+state.min_value
ENDIF 

IF state.scaling EQ 1 THEN BEGIN

    offset = state.min_value - $
      (state.max_value - state.min_value) * 0.01

    sfac = alog10((state.max_value - offset) / (state.min_value - offset))
    state.max_value = 10.^(sfac*(sx+sy)+alog10(state.min_value - offset)) $
      + offset
    state.min_value = 10.^(sfac*(sx-sy)+alog10(state.min_value - offset)) $
      + offset
    
ENDIF 

; do this differently for 8 or 24 bit color, to prevent flashing
if (state.bitdepth EQ 8) then begin
    x_setobjgui_set_minmax
    x_setobjgui_displayall
    state.brightness = 0.5      ; reset these
    state.contrast = 0.5
    x_setobjgui_stretchct, state.brightness, state.contrast
endif else begin
    state.brightness = 0.5      ; reset these
    state.contrast = 0.5
    x_setobjgui_stretchct, state.brightness, state.contrast
    x_setobjgui_set_minmax
    x_setobjgui_displayall
endelse

end

;----------------------------------------------------------------------

function x_setobjgui_wavestring

; function to return string with wavelength info for spectral images

common x_setobjgui_state

cd = (*state.astr_ptr).cd[0,0]
crpix = (*state.astr_ptr).crpix[0]
crval = (*state.astr_ptr).crval[0]

cunit = sxpar(*state.head_ptr, 'cunit1')
cunit = strcompress(string(cunit), /remove_all)
if (cunit NE '0') then begin
    cunit = strcompress(strupcase(strmid(cunit,0,1)) + strmid(cunit,1), $
                        /remove_all)
endif else begin
    cunit = ''
endelse

shifta = float(sxpar(*state.head_ptr, 'SHIFTA1'))

wavelength = crval + ((state.coord[0] - crpix) * cd) + (shifta * cd)
wstring = string(wavelength, format='(F8.2)')

wavestring = strcompress('Wavelength:  ' + wstring + ' ' + cunit)

return, wavestring

end

;--------------------------------------------------------------------

pro x_setobjgui_gettrack

; Create the image to display in the track window that tracks
; cursor movements.  Also update the coordinate display and the
; (x,y) and pixel value.

common x_setobjgui_state
common x_setobjgui_images

; Get x and y for center of track window

zcenter = (0 > state.coord < state.image_size)

track = bytarr(11,11)
boxsize=5
xmin = 0 > (zcenter[0] - boxsize)
xmax = (zcenter[0] + boxsize) < (state.image_size[0] - 1) 
ymin = 0 > (zcenter[1] - boxsize) 
ymax = (zcenter[1] + boxsize) < (state.image_size[1] - 1)

startx = abs( (zcenter[0] - boxsize) < 0 )
starty = abs( (zcenter[1] - boxsize) < 0 ) 

track[startx,starty] = scaled_image[xmin:xmax,ymin:ymax]
track_image = rebin(track, $
                    state.track_window_size, state.track_window_size, $
                    /sample)

x_setobjgui_setwindow, state.track_window_id
tv, track_image

; Overplot an X on the central pixel in the track window, to show the
; current mouse position

; Changed central x to be green always
plots, [0.46, 0.54], [0.46, 0.54], /normal, color = state.box_color, psym=0
plots, [0.46, 0.54], [0.54, 0.46], /normal, color = state.box_color, psym=0

; update location bar with x, y, and pixel value

loc_string = $
  string(state.flipcoord[0], $
         state.flipcoord[1], $
         main_image[state.flipcoord[0], $
                    state.flipcoord[1]], $
         format = '("(",i4,",",i4,") ",g14.7)') 
widget_control, state.location_bar_id, set_value = loc_string

; update wavelength and sigma

if state.wavesig_flg NE 0 then begin
    wavsig_string = string(0.d, wave_image[state.flipcoord[0], $
                                           state.flipcoord[1]], $
                           format = '(g12.5,3x,g14.7)' )
    widget_control, state.wavesig_bar_id, set_value = wavsig_string
endif

; Update coordinate display

if (state.wcstype EQ 'angle') then begin
    xy2ad, state.coord[0], state.coord[1], *(state.astr_ptr), lon, lat

    wcsstring = x_setobjgui_wcsstring(lon, lat, (*state.astr_ptr).ctype,  $
                              state.equinox, state.display_coord_sys, $
                              state.display_equinox, state.display_base60)

;    widget_control, state.wcs_bar_id, set_value = wcsstring

endif    

if (state.wcstype EQ 'lambda') then begin
    wavestring = x_setobjgui_wavestring()
;    widget_control, state.wcs_bar_id, set_value = wavestring
endif

x_setobjgui_resetwindow

end

;----------------------------------------------------------------------

pro x_setobjgui_drawbox, norefresh=norefresh

; routine to draw the box on the pan window, given the current center
; of the display image.

common x_setobjgui_state
common x_setobjgui_images

x_setobjgui_setwindow, state.pan_window_id

view_min = round(state.centerpix - $
        (0.5 * state.draw_window_size / state.zoom_factor)) 
view_max = round(view_min + state.draw_window_size / state.zoom_factor) - 1

; Create the vectors which contain the box coordinates

box_x = float((([view_min[0], $
                 view_max[0], $
                 view_max[0], $
                 view_min[0], $
                 view_min[0]]) * state.pan_scale) + state.pan_offset[0]) 

box_y = float((([view_min[1], $
                 view_min[1], $
                 view_max[1], $
                 view_max[1], $
                 view_min[1]]) * state.pan_scale) + state.pan_offset[1]) 

; Redraw the pan image and overplot the box
if (not(keyword_set(norefresh))) then $
    device, copy=[0,0,state.pan_window_size, state.pan_window_size, 0, 0, $
                  state.pan_pixmap]

plots, box_x, box_y, /device, color = state.box_color, psym=0

x_setobjgui_resetwindow

end

;----------------------------------------------------------------------

pro x_setobjgui_pantrack, event

; routine to track the view box in the pan window during cursor motion

common x_setobjgui_state

; get the new box coords and draw the new box

tmp_event = [event.x, event.y] 

newpos = state.pan_offset > tmp_event < $
  (state.pan_offset + (state.image_size * state.pan_scale))

state.centerpix = round( (newpos - state.pan_offset ) / state.pan_scale)

x_setobjgui_drawbox
x_setobjgui_getoffset

end

;----------------------------------------------------------------------

pro x_setobjgui_resize

; Routine to resize the draw window when a top-level resize event
; occurs.

common x_setobjgui_state

widget_control, state.base_id, tlb_get_size=tmp_event

drawpad = (widget_info(state.draw_base_id,/geometry)).xsize - $
  state.draw_window_size[0]

window = (state.base_min_size > tmp_event)

newbase = window - state.base_pad

; Note: don't know why the (-4) is needed below, but it seems
; to work... without it the base becomes the wrong size when the
; colorbar draw widget is resized.

widget_control, state.colorbar_base_id, $
  xsize = newbase[0]-4, ysize = state.colorbar_height + 6

widget_control, state.draw_base_id, $
  xsize = newbase[0], ysize = newbase[1]

newxsize = (widget_info(state.draw_base_id,/geometry)).xsize - drawpad
newysize = (widget_info(state.draw_base_id,/geometry)).ysize - drawpad

widget_control, state.draw_widget_id, $
  xsize = newxsize, ysize = newysize
widget_control, state.colorbar_widget_id, $
  xsize = newxsize, ysize = state.colorbar_height

state.draw_window_size = [newxsize, newysize]

x_setobjgui_colorbar

end

;----------------------------------------------------------------------

pro x_setobjgui_scaleimage

; Create a byte-scaled copy of the image, scaled according to
; the state.scaling parameter.  Add a padding of 5 pixels around the
; image boundary, so that the tracking window can always remain
; centered on an image pixel even if that pixel is at the edge of the
; image.    

common x_setobjgui_state
common x_setobjgui_images

; Since this can take some time for a big image, set the cursor 
; to an hourglass until control returns to the event loop.

widget_control, /hourglass

delvarx, scaled_image 

; JXP
; Allow for flipping

if state.flg_flip EQ 0 then begin

case state.scaling of
    0: scaled_image = $                 ; linear stretch
      bytscl(main_image, $
             /nan, $
             min=state.min_value, $
             max=state.max_value, $
             top = state.ncolors - 1) + 8
    
    1: begin                            ; log stretch
        offset = state.min_value - $
          (state.max_value - state.min_value) * 0.01

        scaled_image = $        
          bytscl( alog10(main_image - offset), $
                  min=alog10(state.min_value - offset), /nan, $
                  max=alog10(state.max_value - offset),  $
                  top=state.ncolors - 1) + 8   
    end
    

    2: scaled_image = $                 ; histogram equalization
      bytscl(hist_equal(main_image, $
                        minv = state.min_value, $    
                        maxv = state.max_value), $
             /nan, top = state.ncolors - 1) + 8
    
endcase

endif else begin  ; FLIP!

case state.flg_flip of 
    1: tmp_image = rotate(main_image, 5)  ; Flip x
    2: tmp_image = rotate(main_image, 7)  ; Flip y
    3: tmp_image = rotate(main_image, 2)  ; Flip x,y
    else: 
endcase

case state.scaling of
    0: scaled_image = $                 ; linear stretch
      bytscl(tmp_image, $
             /nan, $
             min=state.min_value, $
             max=state.max_value, $
             top = state.ncolors - 1) + 8
    
    1: begin                            ; log stretch
        offset = state.min_value - $
          (state.max_value - state.min_value) * 0.01

        scaled_image = $        
          bytscl( alog10(tmp_image - offset), $
                  min=alog10(state.min_value - offset), /nan, $
                  max=alog10(state.max_value - offset),  $
                  top=state.ncolors - 1) + 8   
    end
    

    2: scaled_image = $                 ; histogram equalization
      bytscl(hist_equal(tmp_image, $
                        minv = state.min_value, $    
                        maxv = state.max_value), $
             /nan, top = state.ncolors - 1) + 8
    
endcase
delvarx, tmp_image
endelse


end

;----------------------------------------------------------------------

pro x_setobjgui_getstats, align=align

; Get basic image stats: min and max, and size.
; set slign keyword to preserve alignment of previous image

common x_setobjgui_state
common x_setobjgui_images

; this routine operates on main_image, which is in the
; x_setobjgui_images common block

widget_control, /hourglass

state.image_size = [ (size(main_image))[1], (size(main_image))[2] ]

state.image_min = min(main_image, max=maxx, /nan)
state.image_max = maxx

if (state.min_value GE state.max_value) then begin
    state.min_value = state.min_value - 1
    state.max_value = state.max_value + 1
endif

; zero the current display position on the center of the image,
; unless user selected /align keyword

state.coord = round(state.image_size / 2.)
IF NOT keyword_set(align) THEN state.centerpix = round(state.image_size / 2.)
x_setobjgui_getoffset

; Clear all plot annotations
x_setobjguierase, /norefresh  

end

;-------------------------------------------------------------------

pro x_setobjgui_setwindow, windowid

; replacement for wset.  Reads the current active window first.
; This should be used when the currently active window is an external
; (i.e. non-x_setobjgui) idl window.  Use x_setobjgui_setwindow to set the window to
; one of the x_setobjgui window, then display something to that window, then
; use x_setobjgui_resetwindow to set the current window back to the currently
; active external window.  Make sure that device is not set to
; postscript, because if it is we can't display anything.

common x_setobjgui_state

if (!d.name NE 'PS') then begin
    state.active_window_id = !d.window
    wset, windowid
endif

end

;---------------------------------------------------------------------

pro x_setobjgui_resetwindow

; reset to current active window

common x_setobjgui_state

; The empty command used below is put there to make sure that all
; graphics to the previous x_setobjgui window actually get displayed to screen
; before we wset to a different window.  Without it, some line
; graphics would not actually appear on screen.

if (!d.name NE 'PS') then begin
    empty
    wset, state.active_window_id
endif

end

;------------------------------------------------------------------

pro x_setobjgui_getwindow

; get currently active window id

common x_setobjgui_state

if (!d.name NE 'PS') then begin
    state.active_window_id = !d.window
endif
end

;--------------------------------------------------------------------
;    Fits file reading routines
;--------------------------------------------------------------------

pro x_setobjgui_readfits, fitsfilename=fitsfilename, newimage=newimage

; Read in a new image when user goes to the File->ReadFits menu.
; Do a reasonable amount of error-checking first, to prevent unwanted
; crashes. 

common x_setobjgui_state
common x_setobjgui_images

newimage = 0
cancelled = 0
if (n_elements(fitsfilename) EQ 0) then window = 1 else window = 0

; If fitsfilename hasn't been passed to this routine, get filename
; from dialog_pickfile.
if (n_elements(fitsfilename) EQ 0) then begin
    fitsfile = $
      dialog_pickfile(filter = '*.fits', $
                      group = state.base_id, $
                      /must_exist, $
                      /read, $
                      path = state.current_dir, $
                      get_path = tmp_dir, $
                      title = 'Select Fits Image')        
    if (tmp_dir NE '') then state.current_dir = tmp_dir
    if (fitsfile EQ '') then return ; 'cancel' button returns empty string
endif else begin
    fitsfile = fitsfilename
endelse

; Get fits header so we know what kind of image this is.
head = headfits(fitsfile)

; Check validity of fits file header 
if (n_elements(strcompress(head, /remove_all)) LT 2) then begin
    x_setobjgui_message, 'File does not appear to be a valid fits image!', $
      window = window, msgtype = 'error'
    return
endif
if (!ERR EQ -1) then begin
    x_setobjgui_message, $
      'Selected file does not appear to be a valid FITS image!', $
      msgtype = 'error', window = window
    return
endif

; Two system variable definitions are needed in order to run fits_info
defsysv,'!TEXTOUT',1
defsysv,'!TEXTUNIT',0

; Find out if this is a fits extension file, and how many extensions
; JXP kludge to deal with error in fits_info
numext = 0
;fits_info, fitsfile, n_ext = numext, /silent
instrume = strcompress(string(sxpar(head, 'INSTRUME')), /remove_all)
origin = strcompress(sxpar(head, 'ORIGIN'), /remove_all)
naxis = sxpar(head, 'NAXIS')

; Make sure it's not a 1-d spectrum
if (numext EQ 0 AND naxis LT 2) then begin
    x_setobjgui_message, 'Selected file is not a 2-d FITS image!', $
      window = window, msgtype = 'error'
    return
endif

state.title_extras = ''

; Now call the subroutine that knows how to read in this particular
; data format:

if ((numext GT 0) AND (instrume NE 'WFPC2')) then begin
    x_setobjgui_fitsext_read, fitsfile, numext, head, cancelled
endif else if ((instrume EQ 'WFPC2') AND (naxis EQ 3)) then begin
    x_setobjgui_wfpc2_read, fitsfile, head, cancelled
endif else if ((naxis EQ 3) AND (origin EQ '2MASS')) then begin
    x_setobjgui_2mass_read, fitsfile, head, cancelled
endif else begin
    x_setobjgui_plainfits_read, fitsfile, head, cancelled
endelse

if (cancelled EQ 1) then return

; Make sure it's a 2-d image
if ( (size(main_image))[0] NE 2 ) then begin
    x_setobjgui_message, 'Selected file is not a 2-D fits image!', $
      msgtype = 'error', window = window
    main_image = fltarr(512, 512)
    newimage = 1
    return
endif

widget_control, /hourglass

state.imagename = fitsfile
x_setobjgui_setheader, head
newimage = 1

end

;----------------------------------------------------------
;  Subroutines for reading specific data formats
;---------------------------------------------------------------

pro x_setobjgui_fitsext_read, fitsfile, numext, head, cancelled

; Fits reader for fits extension files

common x_setobjgui_state
common x_setobjgui_images

numlist = ''
for i = 1, numext do begin
    numlist = strcompress(numlist + string(i) + '|', /remove_all)
endfor

numlist = strmid(numlist, 0, strlen(numlist)-1)

droptext = strcompress('0, droplist, ' + numlist + $
                       ', label_left=Select Extension:, set_value=0')

formdesc = ['0, button, Read Primary Image, quit', $
            '0, label, OR:', $
            droptext, $
            '0, button, Read Fits Extension, quit', $
            '0, button, Cancel, quit']

textform = cw_form(formdesc, /column, $
                   title = 'Fits Extension Selector')

if (textform.tag4 EQ 1) then begin  ; cancelled 
    cancelled = 1
    return                         
endif

if (textform.tag3 EQ 1) then begin   ;extension selected
    extension = long(textform.tag2) + 1
endif else begin
    extension = 0               ; primary image selected
endelse

; Make sure it's not a fits table: this would make mrdfits crash
head = headfits(fitsfile, exten=extension)
xten = strcompress(sxpar(head, 'XTENSION'), /remove_all)
if (xten EQ 'BINTABLE') then begin
    x_setobjgui_message, 'File appears to be a FITS table, not an image.', $
      msgtype='error', /window
    cancelled = 1
    return
endif

if (extension GE 1) then begin
    state.title_extras = strcompress('Extension ' + string(extension))
endif else begin
    state.title_extras = 'Primary Image'
endelse

; Read in the image
delvarx, main_image
main_image = mrdfits(fitsfile, extension, head, /silent, /fscale) 

end

;----------------------------------------------------------------

pro x_setobjgui_plainfits_read, fitsfile, head, cancelled

common x_setobjgui_state
common x_setobjgui_images

; Fits reader for plain fits files, no extensions.

delvarx, main_image
main_image = mrdfits(fitsfile, 0, head, /silent, /fscale) 

end

;----------------------------------------------------------------------

pro x_setobjgui_writeps

; Writes an encapsulated postscript file of the current display.
; Calls cmps_form to get postscript file parameters.

common x_setobjgui_state
common x_setobjgui_images
common x_setobjgui_color

widget_control, /hourglass

view_min = round(state.centerpix - $
                  (0.5 * state.draw_window_size / state.zoom_factor))
view_max = round(view_min + state.draw_window_size / state.zoom_factor)

xsize = (state.draw_window_size[0] / state.zoom_factor) > $
  (view_max[0] - view_min[0] + 1)
ysize = (state.draw_window_size[1] / state.zoom_factor) > $
  (view_max[1] - view_min[1] + 1)


aspect = float(ysize) / float(xsize)
fname = strcompress(state.current_dir + 'x_setobjgui.ps', /remove_all)

tvlct, rr, gg, bb, 8, /get

forminfo = cmps_form(cancel = canceled, create = create, $
                     aspect = aspect, parent = state.base_id, $
                     /preserve_aspect, $
                     xsize = 6.0, ysize = 6.0 * aspect, $
                     /color, /encapsulated, $
                     /nocommon, papersize='Letter', $
                     bits_per_pixel=8, $
                     filename = fname, $
                     button_names = ['Create PS File'])

if (canceled) then return
if (forminfo.filename EQ '') then return

tvlct, rr, gg, bb, 8

tmp_result = findfile(forminfo.filename, count = nfiles)

result = ''
if (nfiles GT 0) then begin
    mesg = strarr(2)
    mesg[0] = 'Overwrite existing file:'
    tmp_string = strmid(forminfo.filename, rstrpos(forminfo.filename, '/') + 1)
    mesg[1] = strcompress(tmp_string + '?', /remove_all)
    result =  dialog_message(mesg, $
                             /default_no, $
                             dialog_parent = state.base_id, $
                             /question)                 
endif

if (strupcase(result) EQ 'NO') then return
    
widget_control, /hourglass

screen_device = !d.name

; In 8-bit mode, the screen color table will have fewer than 256
; colors.  Stretch out the existing color table to 256 colors for the
; postscript plot.

set_plot, 'ps'

device, _extra = forminfo

tvlct, rr, gg, bb, 8, /get

rn = congrid(rr, 248)
gn = congrid(gg, 248)
bn = congrid(bb, 248)

tvlct, temporary(rn), temporary(gn), temporary(bn), 8

; Make a full-resolution version of the display image, accounting for
; scalable pixels in the postscript output

newdisplay = bytarr(xsize, ysize)

startpos = abs(round(state.offset) < 0)

view_min = (0 > view_min < (state.image_size - 1)) 
view_max = (0 > view_max < (state.image_size - 1)) 

dimage = bytscl(scaled_image[view_min[0]:view_max[0], $
                                 view_min[1]:view_max[1]], $
                    top = 247, min=8, max=(!d.table_size-1)) + 8


newdisplay[startpos[0], startpos[1]] = temporary(dimage)

; if there's blank space around the image border, keep it black
tv, newdisplay
x_setobjgui_plotall

if (state.frame EQ 1) then begin    ; put frame around image
    plot, [0], [0], /nodata, position=[0,0,1,1], $
      xrange=[0,1], yrange=[0,1], xstyle=5, ystyle=5, /noerase
    boxx = [0,0,1,1,0,0]
    boxy = [0,1,1,0,0,1]
    oplot, boxx, boxy, color=0, thick=state.framethick
endif

tvlct, temporary(rr), temporary(gg), temporary(bb), 8


device, /close
set_plot, screen_device


end

;----------------------------------------------------------------------
;       routines for defining the color maps
;----------------------------------------------------------------------

pro x_setobjgui_stretchct, brightness, contrast,  getmouse = getmouse

; routine to change color stretch for given values of 
; brightness and contrast.
; Complete rewrite 2000-Sep-21 - Doug Finkbeiner
; This routine is now shorter and easier to understand.  

common x_setobjgui_state
common x_setobjgui_color

; if GETMOUSE then assume mouse positoin passed; otherwise ignore
; inputs

if (keyword_set(getmouse)) then begin 
   state.brightness = brightness/float(state.draw_window_size[0])
   state.contrast = contrast/float(state.draw_window_size[1])
endif

x = state.brightness*(state.ncolors-1)
y = state.contrast*(state.ncolors-1) > 2   ; Minor change by AJB 
high = x+y & low = x-y
diff = (high-low) > 1

slope = float(state.ncolors-1)/diff ;Scale to range of 0 : nc-1
intercept = -slope*low
p = long(findgen(state.ncolors)*slope+intercept) ;subscripts to select
tvlct, r_vector[p], g_vector[p], b_vector[p], 8

end

;------------------------------------------------------------------

pro x_setobjgui_initcolors

; Load a simple color table with the basic 8 colors in the lowest 
; 8 entries of the color table.  Also set top color to white.

common x_setobjgui_state

rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
tvlct, 255*rtiny, 255*gtiny, 255*btiny

tvlct, [255],[255],[255], !d.table_size-1

end

;--------------------------------------------------------------------

pro x_setobjgui_getct, tablenum

; Read in a pre-defined color table, and invert if necessary.

common x_setobjgui_color
common x_setobjgui_state
common x_setobjgui_images


loadct, tablenum, /silent,  bottom=8
tvlct, r, g, b, 8, /get

x_setobjgui_initcolors

r = r[0:state.ncolors-2]
g = g[0:state.ncolors-2]
b = b[0:state.ncolors-2]

if (state.invert_colormap EQ 1) then begin
r = reverse(r)
g = reverse(g)
b = reverse(b)
endif

r_vector = r
g_vector = g
b_vector = b

x_setobjgui_stretchct, state.brightness, state.contrast
if (state.bitdepth EQ 24 AND (n_elements(pan_image) GT 10) ) then $
  x_setobjgui_refresh

end

;--------------------------------------------------------------------


function x_setobjgui_polycolor, p

; Routine to return an vector of length !d.table_size-8,
; defined by a 5th order polynomial.   Called by x_setobjgui_makect
; to define new color tables in terms of polynomial coefficients.

common x_setobjgui_state

x = findgen(256)

y = p[0] + x * p[1] + x^2 * p[2] + x^3 * p[3] + x^4 * p[4] + x^5 * p[5]

w = where(y GT 255, nw)
if (nw GT 0) then y(w) = 255

w =  where(y LT 0, nw)
if (nw GT 0) then y(w) = 0

z = congrid(y, state.ncolors)

return, z
end

;----------------------------------------------------------------------

pro x_setobjgui_makect, tablename

; Define new color tables here.  Invert if necessary.

common x_setobjgui_state
common x_setobjgui_color

case tablename of
    'ATV Special': begin
        r = x_setobjgui_polycolor([39.4609, $
                           -5.19434, $
                           0.128174, $
                           -0.000857115, $
                           2.23517e-06, $
                           -1.87902e-09])
        
        g = x_setobjgui_polycolor([-15.3496, $
                           1.76843, $
                           -0.0418186, $
                           0.000308216, $
                           -6.07106e-07, $
                           0.0000])
        
        b = x_setobjgui_polycolor([0.000, $ 
                           12.2449, $
                           -0.202679, $
                           0.00108027, $
                           -2.47709e-06, $
                           2.66846e-09])

   end

; add more color table definitions here as needed...
    else: return

endcase

if (state.invert_colormap EQ 1) then begin
r = reverse(r)
g = reverse(g)
b = reverse(b)
endif

r_vector = temporary(r)
g_vector = temporary(g)
b_vector = temporary(b)

x_setobjgui_stretchct, state.brightness, state.contrast
if (state.bitdepth EQ 24) then x_setobjgui_refresh

end

;----------------------------------------------------------------------

function x_setobjgui_icolor, color

; Routine to reserve the bottom 8 colors of the color table
; for plot overlays and line plots.

if (n_elements(color) EQ 0) then return, 1

ncolor = N_elements(color)

; If COLOR is a string or array of strings, then convert color names
; to integer values
if (size(color,/tname) EQ 'STRING') then begin ; Test if COLOR is a string
    
; Detemine the default color for the current device
    if (!d.name EQ 'X') then defcolor = 7 $ ; white for X-windows
    else defcolor = 0           ; black otherwise
    
    icolor = 0 * (color EQ 'black') $
      + 1 * (color EQ 'red') $
      + 2 * (color EQ 'green') $
      + 3 * (color EQ 'blue') $
      + 4 * (color EQ 'cyan') $
      + 5 * (color EQ 'magenta') $
      + 6 * (color EQ 'yellow') $
      + 7 * (color EQ 'white') $
      + defcolor * (color EQ 'default')
    
endif else begin
    icolor = long(color)
endelse

return, icolor
end 
 
;---------------------------------------------------------------------
;    routines dealing with image header, title,  and related info
;--------------------------------------------------------------------

pro x_setobjgui_settitle

; Update title bar with the image file name

common x_setobjgui_state

if (state.imagename EQ '') then begin
    widget_control, state.base_id, tlb_set_title = 'x_setobjgui'
endif else begin
    slash = rstrpos(state.imagename, '/')
    ; inserted untested code for MacOS and Windows delimiters
    if (slash EQ -1) then slash = rstrpos(state.imagename, '\')
    if (slash EQ -1) then slash = rstrpos(state.imagename, ':')

    if (slash NE -1) then name = strmid(state.imagename, slash+1) $
      else name = state.imagename
    title = strcompress('x_setobjgui:  '+ name + '  ' + state.title_extras)
    widget_control, state.base_id, tlb_set_title = title
endelse

end

;----------------------------------------------------------------------

pro x_setobjgui_setheader, head

; Routine to keep the image header using a pointer to a 
; heap variable.  If there is no header (i.e. if x_setobjgui has just been
; passed a data array rather than a filename), then make the
; header pointer a null pointer.  Get astrometry info from the 
; header if available.  If there's no astrometry information, set 
; state.astr_ptr to be a null pointer.

common x_setobjgui_state

; Kill the header info window when a new image is read in

if (xregistered('x_setobjgui_headinfo')) then begin
    widget_control, state.headinfo_base_id, /destroy
endif

if (xregistered('x_setobjgui_stats')) then begin
    widget_control, state.stats_base_id, /destroy
endif

if (n_elements(head) LE 1) then begin
; If there's no image header...
    state.wcstype = 'none'
    ptr_free, state.head_ptr
    state.head_ptr = ptr_new()
    ptr_free, state.astr_ptr
    state.astr_ptr = ptr_new()
;    widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    return
endif

ptr_free, state.head_ptr
state.head_ptr = ptr_new(head)

; Get astrometry information from header, if it exists
ptr_free, state.astr_ptr        ; kill previous astrometry info
state.astr_ptr = ptr_new()
extast, head, astr, noparams

; No valid astrometry in header
if (noparams EQ -1) then begin 
;    widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    state.wcstype = 'none'
    return
endif

; coordinate types that we can't use:
if ( (strcompress(string(astr.ctype[0]), /remove_all) EQ 'PIXEL') $
     or (strcompress(string(astr.ctype[0]), /remove_all) EQ '') ) then begin
;    widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    state.wcstype = 'none'
    return
endif

; Image is a 2-d calibrated spectrum (probably from stis):
if (astr.ctype[0] EQ 'LAMBDA') then begin
    state.wcstype = 'lambda'
    state.astr_ptr = ptr_new(astr)
;    widget_control, state.wcs_bar_id, set_value = '                 '
    return
endif

; Good astrometry info in header:
state.wcstype = 'angle'
;widget_control, state.wcs_bar_id, set_value = '                 '

; Check for GSS type header  
if strmid( astr.ctype[0], 5, 3) EQ 'GSS' then begin
    hdr1 = head
    gsss_STDAST, hdr1
    extast, hdr1, astr, noparams
endif

; Create a pointer to the header info
state.astr_ptr = ptr_new(astr)

; Get the equinox of the coordinate system
equ = get_equinox(head, code)
if (code NE -1) then begin
    if (equ EQ 2000.0) then state.equinox = 'J2000'
    if (equ EQ 1950.0) then state.equinox = 'B1950'
    if (equ NE 2000.0 and equ NE 1950.0) then $
      state.equinox = string(equ, format = '(f6.1)')
endif else begin
    IF (strmid(astr.ctype[0], 0, 4) EQ 'GLON') THEN BEGIN 
        state.equinox = 'J2000' ; (just so it is set)
    ENDIF ELSE BEGIN                          
        ptr_free, state.astr_ptr    ; clear pointer
        state.astr_ptr = ptr_new()
        state.equinox = 'J2000'
        state.wcstype = 'none'
;        widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    ENDELSE 
endelse

; Set default display to native system in header
state.display_equinox = state.equinox
state.display_coord_sys = strmid(astr.ctype[0], 0, 4)

end

;---------------------------------------------------------------------


pro x_setobjgui_headinfo

common x_setobjgui_state

; If there's no header, kill the headinfo window and exit this
; routine.
if (not(ptr_valid(state.head_ptr))) then begin
    if (xregistered('x_setobjgui_headinfo')) then begin
        widget_control, state.headinfo_base_id, /destroy
    endif

    x_setobjgui_message, 'No header information available for this image!', $
      msgtype = 'error', /window
    return
endif


; If there is header information but not headinfo window,
; create the headinfo window.
if (not(xregistered('x_setobjgui_headinfo', /noshow))) then begin

    headinfo_base = $
      widget_base(/floating, $
                  /base_align_right, $
                  group_leader = state.base_id, $
                  /column, $
                  title = 'x_setobjgui image header information', $
                  uvalue = 'headinfo_base')
    state.headinfo_base_id = headinfo_base

    h = *(state.head_ptr)

    headinfo_text = widget_text(headinfo_base, $
                            /scroll, $
                            value = h, $
                            xsize = 85, $
                            ysize = 24)
    
    headinfo_done = widget_button(headinfo_base, $
                              value = 'Done', $
                              uvalue = 'headinfo_done')

    widget_control, headinfo_base, /realize
    xmanager, 'x_setobjgui_headinfo', headinfo_base, /no_block

endif


end

;---------------------------------------------------------------------

pro x_setobjgui_headinfo_event, event

common x_setobjgui_state

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'headinfo_done': widget_control, event.top, /destroy
    else:
endcase

end

;----------------------------------------------------------------------
;             routines to do plot overlays
;----------------------------------------------------------------------

pro x_setobjgui_plot1plot, iplot
common x_setobjgui_pdata
common x_setobjgui_state

; Plot a point or line overplot on the image

x_setobjgui_setwindow, state.draw_window_id

widget_control, /hourglass

oplot, [(*(plot_ptr[iplot])).x], [(*(plot_ptr[iplot])).y], $
  _extra = (*(plot_ptr[iplot])).options

x_setobjgui_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

pro x_setobjgui_plot1text, iplot
common x_setobjgui_pdata
common x_setobjgui_state

; Plot a text overlay on the image
x_setobjgui_setwindow, state.draw_window_id

widget_control, /hourglass

xyouts, (*(plot_ptr[iplot])).x, (*(plot_ptr[iplot])).y, $
  (*(plot_ptr[iplot])).text, _extra = (*(plot_ptr[iplot])).options

x_setobjgui_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

pro x_setobjgui_plot1contour, iplot
common x_setobjgui_pdata
common x_setobjgui_state

; Overplot contours on the image

x_setobjgui_setwindow, state.draw_window_id
widget_control, /hourglass

xrange = !x.crange
yrange = !y.crange

; The following allows for 2 conditions, depending upon whether X and Y
; are set

dims = size( (*(plot_ptr[iplot])).z,/dim )

if (size( (*(plot_ptr[iplot])).x,/N_elements ) EQ dims[0] $
    AND size( (*(plot_ptr[iplot])).y,/N_elements) EQ dims[1] ) then begin
    
   contour, (*(plot_ptr[iplot])).z, (*(plot_ptr[iplot])).x, $
      (*(plot_ptr[iplot])).y, $
      position=[0,0,1,1], xrange=xrange, yrange=yrange, $
      xstyle=5, ystyle=5, /noerase, $
      _extra = (*(plot_ptr[iplot])).options
    
endif else begin
    
    contour, (*(plot_ptr[iplot])).z, $
      position=[0,0,1,1], xrange=xrange, yrange=yrange, $
      xstyle=5, ystyle=5, /noerase, $
      _extra = (*(plot_ptr[iplot])).options
          
endelse

x_setobjgui_resetwindow
state.newrefresh=1
end

;---------------------------------------------------------------------

pro x_setobjgui_plot1compass, iplot

; Uses idlastro routine arrows to plot compass arrows.

common x_setobjgui_pdata
common x_setobjgui_state

x_setobjgui_setwindow, state.draw_window_id

widget_control, /hourglass

arrows, *(state.head_ptr), $
  (*(plot_ptr[iplot])).x, $
  (*(plot_ptr[iplot])).y, $
  thick = (*(plot_ptr[iplot])).thick, $
  charsize = (*(plot_ptr[iplot])).charsize, $
  arrowlen = (*(plot_ptr[iplot])).arrowlen, $
  color = (*(plot_ptr[iplot])).color, $
  notvertex = (*(plot_ptr[iplot])).notvertex, $
  /data

x_setobjgui_resetwindow
state.newrefresh=1
end

;---------------------------------------------------------------------

pro x_setobjgui_plot1scalebar, iplot

; uses modified version of idlastro routine arcbar to plot a scalebar

common x_setobjgui_pdata
common x_setobjgui_state

x_setobjgui_setwindow, state.draw_window_id
widget_control, /hourglass

; routine arcbar doesn't recognize color=0, because it uses 
; keyword_set to check the color.  So we need to set !p.color = 0
; to get black if the user wants color=0

!p.color = 0

x_setobjgui_arcbar, *(state.head_ptr), $
  (*(plot_ptr[iplot])).arclen, $
  position = (*(plot_ptr[iplot])).position, $
  thick = (*(plot_ptr[iplot])).thick, $
  size = (*(plot_ptr[iplot])).size, $
  color = (*(plot_ptr[iplot])).color, $
  seconds = (*(plot_ptr[iplot])).seconds, $
  /data

x_setobjgui_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

pro x_setobjgui_plotwindow
common x_setobjgui_state

x_setobjgui_setwindow, state.draw_window_id

; Set plot window
xrange=[state.offset[0], $
 state.offset[0] + state.draw_window_size[0] / state.zoom_factor] - 0.5
yrange=[state.offset[1], $
 state.offset[1] + state.draw_window_size[1] / state.zoom_factor] - 0.5

plot, [0], [0], /nodata, position=[0,0,1,1], $
 xrange=xrange, yrange=yrange, xstyle=5, ystyle=5, /noerase

x_setobjgui_resetwindow
end

;----------------------------------------------------------------------

pro x_setobjgui_plotall
common x_setobjgui_state
common x_setobjgui_pdata

; Routine to overplot all line, text, and contour plots

if (nplot EQ 0) then return

x_setobjgui_plotwindow

for iplot = 1, nplot do begin
    case (*(plot_ptr[iplot])).type of
        'points'  : x_setobjgui_plot1plot, iplot
        'text'    : x_setobjgui_plot1text, iplot
        'contour' : x_setobjgui_plot1contour, iplot
        'compass' : x_setobjgui_plot1compass, iplot
        'scalebar': x_setobjgui_plot1scalebar, iplot
        else      : print, 'Problem in x_setobjgui_plotall!'   
    endcase
endfor

end

;----------------------------------------------------------------------

pro x_setobjguiplot, x, y, _extra = options
common x_setobjgui_pdata
common x_setobjgui_state

; Routine to read in line plot data and options, store in a heap
; variable structure, and plot the line plot

if (not(xregistered('x_setobjgui', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 1) then begin
   print, 'Too few parameters for ATVPLOT.'
   return
endif

if (n_elements(options) EQ 0) then options = {color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'color', 'red')
   options.color = x_setobjgui_icolor(options.color)

   pstruct = {type: 'points',   $     ; points
              x: x,             $     ; x coordinate
              y: y,             $     ; y coordinate
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)

   x_setobjgui_plotwindow
   x_setobjgui_plot1plot, nplot

endif else begin
   print, 'Too many calls to ATVPLOT.'
endelse

end

;----------------------------------------------------------------------

pro x_setobjguixyouts, x, y, text, _extra = options
common x_setobjgui_pdata
common x_setobjgui_state

; Routine to read in text overplot string and options, store in a heap
; variable structure, and overplot the text

if (not(xregistered('x_setobjgui', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 3) then begin
   print, 'Too few parameters for ATVXYOUTS'
   return
endif

if (n_elements(options) EQ 0) then options = {color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'color', 'red')
   options.color = x_setobjgui_icolor(options.color)

;  set default font to 1
   c = where(tag_names(options) EQ 'FONT', count)
   if (count EQ 0) then options = create_struct(options, 'font', 1)

   pstruct = {type: 'text',   $       ; type of plot 
              x: x,             $     ; x coordinate
              y: y,             $     ; y coordinate
              text: text,       $     ; text to plot
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)

   x_setobjgui_plotwindow
   x_setobjgui_plot1text, nplot

endif else begin
   print, 'Too many calls to ATVPLOT.'
endelse

end

;----------------------------------------------------------------------

pro x_setobjguicontour, z, x, y, _extra = options
common x_setobjgui_pdata
common x_setobjgui_state

; Routine to read in contour plot data and options, store in a heap
; variable structure, and overplot the contours.  Data to be contoured
; need not be the same dataset displayed in the x_setobjgui window, but it
; should have the same x and y dimensions in order to align the
; overplot correctly.

if (not(xregistered('x_setobjgui', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 1) then begin
   print, 'Too few parameters for ATVCONTOUR.'
   return
endif

if (n_params() EQ 1 OR n_params() EQ 2) then begin
    x = 0
    y = 0
endif

if (n_elements(options) EQ 0) then options = {c_color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'C_COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'c_color', 'red')
   options.c_color = x_setobjgui_icolor(options.c_color)

   pstruct = {type: 'contour',  $     ; type of plot
              z: z,             $     ; z values
              x: x,             $     ; x coordinate
              y: y,             $     ; y coordinate
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)

   x_setobjgui_plotwindow
   x_setobjgui_plot1contour, nplot

endif else begin
   print, 'Too many calls to ATVCONTOUR.'
endelse

end

;----------------------------------------------------------------------

pro x_setobjguierase, nerase, norefresh = norefresh
common x_setobjgui_pdata

; Routine to erase line plots from ATVPLOT, text from ATVXYOUTS, and
; contours from ATVCONTOUR.

if (n_params() LT 1) then begin
    nerase = nplot
endif else begin
    if (nerase GT nplot) then nerase = nplot
endelse

for iplot = nplot - nerase + 1, nplot do begin
    ptr_free, plot_ptr[iplot]
    plot_ptr[iplot] = ptr_new()
endfor

nplot = nplot - nerase

if (NOT keyword_set(norefresh)) then x_setobjgui_refresh

end

;----------------------------------------------------------------------

pro x_setobjgui_textlabel

; widget front end for x_setobjguixyouts

formdesc = ['0, text, , label_left=Text: , width=15', $
            '0, integer, 0, label_left=x: ', $
            '0, integer, 0, label_left=y: ', $
            '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
            '0, float, 2.0, label_left=Charsize: ', $
            '0, integer, 1, label_left=Charthick: ', $
            '0, integer, 0, label_left=Orientation: ', $
            '1, base, , row', $
            '0, button, Cancel, quit', $
            '0, button, DrawText, quit']
            
textform = cw_form(formdesc, /column, $
                   title = 'x_setobjgui text label')

if (textform.tag9 EQ 1) then begin
; switch red and black indices
    case textform.tag3 of
        0: labelcolor = 1
        1: labelcolor = 0
        else: labelcolor = textform.tag3
    endcase

    x_setobjguixyouts, textform.tag1, textform.tag2, textform.tag0, $
      color = labelcolor, charsize = textform.tag4, $
      charthick = textform.tag5, orientation = textform.tag6
endif

end

;---------------------------------------------------------------------

pro x_setobjgui_oplotcontour

; widget front end for x_setobjguicontour

common x_setobjgui_state
common x_setobjgui_images

minvalstring = strcompress('0, float, ' + string(state.min_value) + $
                           ', label_left=MinValue: , width=15 ')
maxvalstring = strcompress('0, float, ' + string(state.max_value) + $
                           ', label_left=MaxValue: , width=15')

formdesc = ['0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
;            '0, float, 1.0, label_left=Charsize: ', $
;            '0, integer, 1, label_left=Charthick: ', $
            '0, droplist, solid|dotted|dashed|dashdot|dashdotdotdot|longdash, label_left=Linestyle: , set_value=0', $
            '0, integer, 1, label_left=LineThickness: ', $
            minvalstring, $
            maxvalstring, $
            '0, integer, 6, label_left=NLevels: ', $
            '1, base, , row,', $
            '0, button, Cancel, quit', $
            '0, button, DrawContour, quit']
            
cform = cw_form(formdesc, /column, $
                   title = 'x_setobjgui text label')


if (cform.tag8 EQ 1) then begin
; switch red and black indices
    case cform.tag0 of
        0: labelcolor = 1
        1: labelcolor = 0
        else: labelcolor = cform.tag0
    endcase

    x_setobjguicontour, main_image, c_color = labelcolor, $
;      c_charsize = cform.tag1, c_charthick = cform.tag2, $
      c_linestyle = cform.tag1, $
      c_thick = cform.tag2, $
      min_value = cform.tag3, max_value = cform.tag4, $, 
      nlevels = cform.tag5
endif

end

;---------------------------------------------------------------------

pro x_setobjgui_setscalebar

; Routine to prompt user for scalebar parameters

common x_setobjgui_state
common x_setobjgui_images
common x_setobjgui_pdata

if (nplot GE maxplot) then begin
    x_setobjgui_message, 'Total allowed number of overplots exceeded.', $
      mesgtype = 'error', /window
    return
endif
    

if (state.wcstype NE 'angle') then begin 
    x_setobjgui_message, 'Cannot get coordinate info for this image!', $
      mesgtype = 'error', /window
    return
endif

view_min = round(state.centerpix - $
        (0.5 * state.draw_window_size / state.zoom_factor)) 
view_max = round(view_min + state.draw_window_size / state.zoom_factor) - 1

xpos = string(round(view_min[0] + 0.75 * (view_max[0] - view_min[0])))
ypos = string(round(view_min[1] + 0.15 * (view_max[1] - view_min[1])))

xposstring = strcompress('0,integer,'+xpos+',label_left=X (left end of bar): ')
yposstring = strcompress('0,integer,'+ypos+',label_left=Y (center of bar): ')

formdesc = [ $
             xposstring, $
             yposstring, $
             '0, float, 5.0, label_left=BarLength: ', $
             '0, droplist, arcsec|arcmin, label_left=Units:,set_value=0', $
             '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
             '0, integer, 1, label_left=LineThickness: ', $
             '0, float, 1, label_left=Charsize: ', $
             '1, base, , row,', $
             '0, button, Cancel, quit', $
             '0, button, DrawScalebar, quit']
            
cform = cw_form(formdesc, /column, $
                   title = 'x_setobjgui scalebar properties')

if (cform.tag8 EQ 1) then return

; switch red and black indices
case cform.tag4 of
    0: labelcolor = 1
    1: labelcolor = 0
    else: labelcolor = cform.tag4
endcase


cform.tag0 = 0 > cform.tag0 < (state.image_size[0] - 1)
cform.tag1 = 0 > cform.tag1 < (state.image_size[1] - 1)
cform.tag3 = abs(cform.tag3 - 1)  ; set default to be arcseconds

arclen = cform.tag2
if (float(round(arclen)) EQ arclen) then arclen = round(arclen)

pstruct = {type: 'scalebar',  $  ; type of plot
           arclen: arclen, $
           seconds: cform.tag3, $
           position: [cform.tag0,cform.tag1], $ 
           color: labelcolor, $
           thick: cform.tag5, $
           size: cform.tag6 $
          }

nplot = nplot + 1
plot_ptr[nplot] = ptr_new(pstruct)

x_setobjgui_plotwindow
x_setobjgui_plot1scalebar, nplot

end

;---------------------------------------------------------------------
;          routines for drawing in the lineplot window
;---------------------------------------------------------------------

pro x_setobjgui_lineplot_init

; This routine creates the window for line plots

common x_setobjgui_state

state.lineplot_base_id = $
  widget_base(/floating, $
              group_leader = state.base_id, $
              /column, $
              /base_align_right, $
              title = 'x_setobjgui plot', $
              /tlb_size_events, $
              uvalue = 'lineplot_base')

state.lineplot_widget_id = $
  widget_draw(state.lineplot_base_id, $
              frame = 0, $
              scr_xsize = state.lineplot_size[0], $
              scr_ysize = state.lineplot_size[1], $
              uvalue = 'lineplot_window')

lbutton_base = $
  widget_base(state.lineplot_base_id, $
              /base_align_bottom, $
              /row)

lineplot_done = $
  widget_button(lbutton_base, $
                value = 'Done', $
                uvalue = 'lineplot_done')

widget_control, state.lineplot_base_id, /realize
widget_control, state.lineplot_widget_id, get_value = tmp_value
state.lineplot_window_id = tmp_value

basegeom = widget_info(state.lineplot_base_id, /geometry)
drawgeom = widget_info(state.lineplot_widget_id, /geometry)

state.lineplot_pad[0] = basegeom.xsize - drawgeom.xsize
state.lineplot_pad[1] = basegeom.ysize - drawgeom.ysize
    
xmanager, 'x_setobjgui_lineplot', state.lineplot_base_id, /no_block

x_setobjgui_resetwindow
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; JXP

pro x_setobjgui_specplot

common x_setobjgui_state
common x_setobjgui_images


; JXP 17-Nov-01
  if state.wavesig_flg MOD 2 NE 1 then return

  state.active_window_id = !d.window

  title = strcompress('Plot of Spectrum at Row ' + $
                    string(state.flipcoord[1]))
  x_splot, wave_image[*,state.flipcoord[1]], $
    main_image[*, state.flipcoord[1]], TITLE=title
  
;  widget_control, state.lineplot_base_id, /clear_events
  
;  x_setobjgui_resetwindow
return
end

;--------------------------------------------------------------------

pro x_setobjgui_rowplot

common x_setobjgui_state
common x_setobjgui_images

;if (not (xregistered('x_setobjgui_lineplot', /noshow))) then begin
;    x_setobjgui_lineplot_init
;endif

;x_setobjgui_setwindow, state.lineplot_window_id
;erase

;plot, main_image[*, state.flipcoord[1]], $
;  xst = 3, yst = 3, psym = 10, $
;  title = strcompress('Plot of row ' + $
;                      string(state.flipcoord[1])), $
;  xtitle = 'Column', $
;  ytitle = 'Pixel Value', $
;  color = 7 

; JXP 17-Nov-01
  state.active_window_id = !d.window

  title = strcompress('Plot of row ' + $
                    string(state.flipcoord[1]))
  x_splot, main_image[*, state.flipcoord[1]], TITLE=title
  
;  widget_control, state.lineplot_base_id, /clear_events
  
;  x_setobjgui_resetwindow
return
end

;--------------------------------------------------------------------

pro x_setobjgui_colplot

common x_setobjgui_state
common x_setobjgui_images

;if (not (xregistered('x_setobjgui_lineplot', /noshow))) then begin
;    x_setobjgui_lineplot_init
;endif

;x_setobjgui_setwindow, state.lineplot_window_id

;erase

;plot, main_image[state.flipcoord[0], *], $
;  xst = 3, yst = 3, psym = 10, $
;  title = strcompress('Plot of column ' + $
;                      string(state.flipcoord[0])), $
;  xtitle = 'Row', $
;  ytitle = 'Pixel Value', $
;  color = 7
;widget_control, state.lineplot_base_id, /clear_events
        
; JXP 17-Nov-01
  state.active_window_id = !d.window

  title = strcompress('Plot of column ' + $
                      string(state.flipcoord[0]))
  x_splot, main_image[state.flipcoord[0], *], TITLE=title
  
end

;--------------------------------------------------------------------

pro x_setobjgui_surfplot

common x_setobjgui_state
common x_setobjgui_images

if (not (xregistered('x_setobjgui_lineplot', /noshow))) then begin
    x_setobjgui_lineplot_init
endif

x_setobjgui_setwindow, state.lineplot_window_id
erase

plotsize = $
  fix(min([50, state.image_size[0]/2., state.image_size[1]/2.]))
center = plotsize > state.coord < (state.image_size - plotsize) 

tmp_string = $
  strcompress('Surface plot of ' + $
              strcompress('['+string(center[0]-plotsize)+ $
                          ':'+string(center[0]+plotsize-1)+ $
                          ','+string(center[1]-plotsize)+ $
                          ':'+string(center[1]+plotsize-1)+ $
                          ']', /remove_all))

surface, $
  main_image[center[0]-plotsize:center[0]+plotsize-1, $
             center[1]-plotsize:center[1]+plotsize-1], $
  title = temporary(tmp_string), $
  xtitle = 'X', ytitle = 'Y', ztitle = 'Pixel Value', $
  color = 7

widget_control, state.lineplot_base_id, /clear_events

x_setobjgui_resetwindow
end

;--------------------------------------------------------------------

pro x_setobjgui_contourplot

common x_setobjgui_state
common x_setobjgui_images

if (not (xregistered('x_setobjgui_lineplot', /noshow))) then begin
    x_setobjgui_lineplot_init
endif

x_setobjgui_setwindow, state.lineplot_window_id
erase

plotsize = $
  fix(min([10, state.image_size[0]/2., state.image_size[1]/2.]))
center = plotsize > state.coord < (state.image_size - plotsize) 

contour_image =  main_image[center[0]-plotsize:center[0]+plotsize-1, $
                            center[1]-plotsize:center[1]+plotsize-1]
if (state.scaling EQ 1) then begin
    contour_image = alog10(contour_image)
    logflag = 'Log'
endif else begin
    logflag = ''
endelse

; JXP -- Centroid!

cntrd, contour_image, plotsize/2., plotsize/2., xcen, ycen, plotsize/4., /silent

tmp_string =  $
  strcompress(logflag + $
              ' Contour plot of ' + $
              strcompress('['+string(round(center[0]-plotsize))+ $
                          ':'+string(round(center[0]+plotsize-1))+ $
                          ','+string(round(center[1]-plotsize))+ $
                          ':'+string(round(center[1]+plotsize-1))+ $
                          ']', /remove_all) +$
              '  Center = ['+$
              strcompress(string(xcen+center[0]-plotsize)+','+$
                          string(ycen+center[1]-plotsize)+']', /remove_all))

contour, temporary(contour_image), $
  nlevels = 10, $
  /follow, $
  title = temporary(tmp_string), $
  xtitle = 'X', ytitle = 'Y', color = 7

widget_control, state.lineplot_base_id, /clear_events
        
x_setobjgui_resetwindow
end

;----------------------------------------------------------------------

pro x_setobjgui_lineplot_event, event

common x_setobjgui_state

widget_control, event.id, get_uvalue = uvalue


case uvalue of
    'lineplot_done': widget_control, event.top, /destroy
    'lineplot_base': begin                       ; Resize event
        x_setobjgui_setwindow, state.lineplot_window_id
        state.lineplot_size = [event.x, event.y]- state.lineplot_pad
        widget_control, state.lineplot_widget_id, $
          xsize = (state.lineplot_size[0] > 100), $
          ysize = (state.lineplot_size[1] > 100)
        x_setobjgui_resetwindow
    end    
else:
endcase

end

;----------------------------------------------------------------------
;                         help window
;---------------------------------------------------------------------

pro x_setobjgui_help
common x_setobjgui_state

h = strarr(110)
i = 0
h[i] =  'x_setobjgui HELP'
i = i + 1
h[i] =  ''
i = i + 1
h[i] =  'MENU BAR:'
i = i + 1
h[i] =  'File->ReadFits:         Read in a new fits image from disk'
i = i + 1
h[i] =  'File->WritePS:          Write a PostScript file of the current display'
i = i + 1
h[i] =  'File->WriteTiff:        Write a tiff image of the current display'
i = i + 1
h[i] =  'File->Quit:             Quits x_setobjgui'
i = i + 1
h[i] =  'ColorMap Menu:          Selects color table'
i = i + 1
h[i] =  'Scaling Menu:           Selects linear, log, or histogram-equalized scaling'
i = i + 1
h[i] =  'Labels->TextLabel:      Brings up a dialog box for text input'
i = i + 1
h[i] =  'Labels->Contour:        Brings up a dialog box for overplotting contours'
i = i + 1
h[i] =  'Labels->Compass:        Draws a compass (requires WCS info in header)'
i = i + 1
h[i] =  'Labels->Scalebar:       Draws a scale bar (requires WCS info in header)'
i = i + 1
h[i] =  'Labels->EraseLast:      Erases the most recent plot label'
i = i + 1
h[i] =  'Labels->EraseAll:       Erases all plot labels'
i = i + 1
h[i] =  'Blink->SetBlink:        Sets the current display to be the blink image'
i = i + 1
h[i] =  '                             for mouse button 1, 2, or 3'
i = i + 1
h[i] =  'ImageInfo->ImageHeader: Display the FITS header, if there is one.'
i = i + 1
h[i] =  'ImageInfo menu also gives a choice of coordinate systems, '
i = i + 1
h[i] =  '    or of native image coordinates (default), for images with a WCS.'
i = i + 1
h[i] =  ''
i = i + 1
h[i] =  'CONTROL PANEL ITEMS:'
i = i + 1
h[i] = 'Min:             shows minimum data value displayed; enter new min value here'
i = i + 1
h[i] = 'Max:             shows maximum data value displayed; enter new max value here'
i = i + 1
h[i] = 'Pan Window:      use mouse to drag the image-view box around'
i = i + 1
h[i] = ''
i = i + 1
h[i] = 'MOUSE MODE SELECTOR:'
i = i + 1
h[i] =  'Color:          sets color-stretch mode:'
i = i + 1
h[i] = '                    With mouse button 1 down, drag mouse to change the color stretch.  '
i = i + 1
h[i] = '                    Move vertically to change contrast, and'
i = i + 1
h[i] = '                         horizontally to change brightness.'
i = i + 1 
h[i] = '                    button 2 or 3: center on current position'
i = i + 2
h[i] = 'BUTTONS:'
i = i + 1
h[i] = 'Invert:          inverts the current color table'
i = i + 1
h[i] = 'Restretch:       sets min and max to preserve display colors while linearizing the color table'
i = i + 1
h[i] = 'AutoScale:       sets min and max to show data values around image median'
i = i + 1
h[i] = 'FullRange:       sets min and max to show the full data range of the image'
i = i + 1
h[i] = 'ZoomIn:          zooms in by x2'
i = i + 1
h[i] = 'ZoomOut:         zooms out by x2'
i = i + 1
h[i] = 'Zoom1:           sets zoom level to original scale'
i = i + 1
h[i] = 'Center:          centers image on display window'
i = i + 1
h[i] = 'Done:            quits x_setobjgui'
i = i + 1
h[i] = ''
i = i + 1
h[i] = 'Keyboard commands in display window:'
i = i + 1
h[i] = '    Numeric keypad (with NUM LOCK on) moves cursor'
i = i + 1
h[i] = '    {: Pan Up'
i = i + 1
h[i] = '    }: Pan down'
i = i + 1
h[i] = '    [: Pan left'
i = i + 1
h[i] = '    ]: Pan right'
i = i + 1
h[i] = '    m: Move closest object trace to cursor location'
i = i + 1
h[i] = '    d: Delete closest non-science object'
i = i + 1
h[i] = '    r: row plot'
i = i + 1
h[i] = '    c: column plot'
i = i + 1
h[i] = '    s: surface plot'
i = i + 1
h[i] = '    t: move top slit'
i = i + 1
h[i] = '    b: move bottom slit'
i = i + 1
h[i] = '    q: quits x_setobjgui'
i = i + 2
h[i] = 'IDL COMMAND LINE HELP:'
i = i + 1
h[i] =  'To pass an array to x_setobjgui:'
i = i + 1
h[i] =  '   x_setobjgui, array_name [, options]'
i = i + 1
h[i] = 'To pass a fits filename to x_setobjgui:'
i = i + 1
h[i] = '    x_setobjgui, fitsfile_name [, options] (enclose filename in single quotes) '
i = i + 1
h[i] = 'Command-line options are: '
i = i + 1
h[i]  = '   [,min = min_value] [,max=max_value] [,/linear] [,/log] [,/histeq]'
i = i + 1
h[i]  = '   [,/block] [,/align] [,/stretch] [,header=header]'
i = i + 2
h[i] = 'To overplot a contour plot on the draw window:'
i = i + 1
h[i] = '    x_setobjguicontour, array_name [, options...]'
i = i + 1
h[i] = 'To overplot text on the draw window: '
i = i + 1
h[i] = '    x_setobjguixyouts, x, y, text_string [, options]  (enclose string in single quotes)'
i = i + 1
h[i] = 'To overplot points or lines on the current plot:'
i = i + 1
h[i] = '    x_setobjguiplot, xvector, yvector [, options]'
i = i + 2
h[i] = 'The options for x_setobjguicontour, x_setobjguixyouts, and x_setobjguiplot are essentially'
i = i + 1
h[i] =  'the same as those for the idl contour, xyouts, and plot commands,'
i = i + 1
h[i] = 'except that data coordinates are always used.' 
i = i + 1
h[i] = 'The default color for overplots is red.'
i = i + 2
h[i] = 'The lowest 8 entries in the color table are:'
i = i + 1
h[i] = '    0 = black'
i = i + 1
h[i] = '    1 = red'
i = i + 1
h[i] = '    2 = green'
i = i + 1
h[i] = '    3 = blue'
i = i + 1
h[i] = '    4 = cyan'
i = i + 1
h[i] = '    5 = magenta'
i = i + 1
h[i] = '    6 = yellow'
i = i + 1
h[i] = '    7 = white'
i = i + 1
h[i] = '    The top entry in the color table is also reserved for white. '
i = i + 2
h[i] = 'Other commands:'
i = i + 1
h[i] = 'x_setobjguierase [, N]:       erases all (or last N) plots and text'
i = i + 1
h[i] = 'x_setobjgui_shutdown:   quits x_setobjgui'
i = i + 1
h[i] = 'NOTE: If x_setobjgui should crash, type x_setobjgui_shutdown at the idl prompt.'
i = i + 5
h[i] = strcompress('ATV.PRO version '+state.version)
i = i + 1
h[i] = 'For full instructions, or to download the most recent version, go to:'
i = i + 1
h[i] = 'http://cfa-www.harvard.edu/~abarth/x_setobjgui/x_setobjgui.html'


if (not (xregistered('x_setobjgui_help', /noshow))) then begin

helptitle = strcompress('x_setobjgui v' + state.version + ' help')

    help_base =  widget_base(/floating, $
                             group_leader = state.base_id, $
                             /column, $
                             /base_align_right, $
                             title = helptitle, $
                             uvalue = 'help_base')

    help_text = widget_text(help_base, $
                            /scroll, $
                            value = h, $
                            xsize = 85, $
                            ysize = 24)
    
    help_done = widget_button(help_base, $
                              value = 'Done', $
                              uvalue = 'help_done')

    widget_control, help_base, /realize
    xmanager, 'x_setobjgui_help', help_base, /no_block
    
endif

end

;----------------------------------------------------------------------

pro x_setobjgui_help_event, event

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'help_done': widget_control, event.top, /destroy
    else:
endcase

end

;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;   Stuff specific to this program

; PLOT
pro x_setobjgui_plotobj

common x_setobjgui_state

  ; Set Window

  x_setobjgui_plotwindow
  x_setobjgui_setwindow, state.draw_window_id
	
; SCIENCE OBJECTS
  sci = where(sobj_obj.obj_id EQ 'a' AND sobj_obj.flg_anly NE 0,nsci)
  ; Points
  oplot, [sobj_obj[sci].xcen], [sobj_obj[sci].ycen], $
    psym=2, color=1
  ; Trace
  for i=0L,nsci-1 do begin
      oplot, findgen(state.sz[0]), sobj_obj[sci[i]].trace[0:state.sz[0]-1],$
        psym=-3, color=2
  endfor

; EXTRA OBJ
  ext = where(sobj_obj.obj_id NE 'a' AND sobj_obj.flg_anly NE 0,next)
  ; Trace
  for i=0L,next-1 do begin
      oplot, findgen(state.sz[0]), sobj_obj[ext[i]].trace[0:state.sz[0]-1],$
        psym=-3, color=3
  endfor

  x_setobjgui_resetwindow
  state.newrefresh=1

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; MOVE THE OBJ
pro x_setobjgui_moveobj

common x_setobjgui_state

  ; Find nearest obj
  gdobj = where(sobj_obj.flg_anly NE 0)
  mn = min(abs(state.coord[1]-sobj_obj[gdobj].trace[state.coord[0]]),imn)

  ; Move to cursor position
  ydelt = state.coord[1]-sobj_obj[gdobj[imn]].trace[state.coord[0]]
  if ydelt GT 50 then begin
      print, 'x_setobjgui_moveobj: Too big a move!'
      return
  endif
  sobj_obj[gdobj[imn]].trace[0:state.sz[0]] =  $
    sobj_obj[gdobj[imn]].trace[0:state.sz[0]] + ydelt
  ; Reset ycen too
  sobj_obj[gdobj[imn]].ycen = sobj_obj[gdobj[imn]].ycen + ydelt 
  
  ; Refresh
  x_setobjgui_refresh
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DELETE THE OBJ
pro x_setobjgui_delobj

common x_setobjgui_state

  ; Find nearest obj
  gdobj = where(sobj_obj.flg_anly NE 0, ngd)
  if ngd EQ 0 then return
  mn = min(abs(state.coord[1]-sobj_obj[gdobj].trace[state.coord[0]]),imn)

  ; Check if Science
  if sobj_obj[gdobj[imn]].obj_id EQ 'a' then begin
      ans = x_guilist(['yes','no'], 'Delete Science Obj?')
      if strtrim(ans,2) EQ 'no' then return
  endif

  ; Move to cursor position
  ydelt = state.coord[1]-sobj_obj[gdobj[imn]].trace[state.coord[0]]
  if ydelt GT 50 then begin
      print, 'x_setobjgui_moveobj: Too far from Obj!'
      return
  endif
  
  ; Delete!
  sobj_obj[gdobj[imn]].flg_anly = 0

  ; Update Obj structure
  a = where(sobj_obj[gdobj[imn]].slit_id EQ sobj_slit.id, na)
  sobj_slit[a[0]].nobj = sobj_slit[a[0]].nobj - 1

  ; Update Obj Name
  if sobj_slit[a[0]].nobj GT 1 then begin
      a = where(sobj_obj.slit_id EQ sobj_obj[gdobj[imn]].slit_id, na)
      for i=1,na-1 do sobj_obj[a[i]].obj_id = x_objnumid(i)
  endif

  ; Refresh
  x_setobjgui_refresh
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ADD AN OBJ
pro x_setobjgui_addobj

common x_setobjgui_state
common x_setobjgui_images

  ; Find Slit
  mnbot = min(abs(state.coord[1]-sobj_slit.yedg_sky[state.coord[0],0]),ibot)
  mntop = min(abs(state.coord[1]-sobj_slit.yedg_sky[state.coord[0],1]),itop)
  if mnbot LT mntop then gd = ibot else gd = itop

  ;; Find science obj
  sci = where(sobj_obj.obj_id EQ 'a' AND sobj_obj.slit_id EQ sobj_slit[gd].id)
  clm = state.coord[0]

  ;; yedg
  yedg = round(sobj_slit[gd].yedg_orig[clm,*])

  ;; Add serendip
  sobj_slit[gd].ypos[sobj_slit[gd].nobj] = state.coord[1]-yedg[0]

  tmp = { specobjstrct }
  tmp.slit_fil = ' '
  tmp.img_fil = ' '
  tmp.spec2d_fil = ' '
  tmp.UT = ' '
  tmp.instr_strct = ' '
  tmp.field = sobj_slit[gd].field
  tmp.flg_anly = 1
  sobj_obj = [sobj_obj,tmp]
  nobj = n_elements(sobj_obj)-1

  ;; Update Obj
  x_slittoobj, sobj_obj, nobj, sobj_slit[gd], clm
  sobj_obj[nobj].ycen = state.coord[1]
  sobj_obj[nobj].aper = [-2.,2.] ; Arbitrary

  ;; ID
  sobj_obj[nobj].obj_id = x_objnumid(sobj_slit[gd].nobj)

  ;; Trace
  sz_img = size(main_image, /dimensions)
  off = sobj_obj[nobj].ycen - sobj_slit[gd].yedg_orig[sobj_obj[nobj].xcen,0]
  sobj_obj[nobj].trace[0:sz_img[0]-1] = sobj_slit[gd].yedg_orig[0:sz_img[0]-1,0]$
        + off

  ;; Increment  
  sobj_slit[gd].nobj = sobj_slit[gd].nobj + 1

  ;; Refresh
  x_setobjgui_refresh
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SLITS

; MOVE
pro x_setobjgui_sngshift, bottom=bottom
common x_setobjgui_state

  ; 
  if keyword_set( BOTTOM ) then begin ; BOTTOM
      mn = min(abs(state.coord[1]-sobj_slit.yedg_sky[state.coord[0],0]),imn)
      ; Shift
      shft = sobj_slit[imn].yedg_sky[state.coord[0],0]-state.coord[1]
      sobj_slit[imn].yedg_sky[*,0] = sobj_slit[imn].yedg_sky[*,0] - shft
  endif else begin ; TOP
      mn = min(abs(state.coord[1]-sobj_slit.yedg_sky[state.coord[0],1]),imn)
      ; Shift
      shft = sobj_slit[imn].yedg_sky[state.coord[0],1]-state.coord[1]
      sobj_slit[imn].yedg_sky[*,1] = sobj_slit[imn].yedg_sky[*,1] - shft
  endelse

  ; Refresh
  x_setobjgui_refresh

  return
end

;;;;;;;;
; PLOT
pro x_setobjgui_plotslit

common x_setobjgui_state

  ; Set Window

  x_setobjgui_plotwindow
  x_setobjgui_setwindow, state.draw_window_id

; Slits
  gdslit = where(sobj_slit.flg_anly NE 0, nslit)
;  nslit = n_elements(sobj_slit)
  ; Points
  for j=0L,nslit-1 do begin
      i = gdslit[j]
      oplot, findgen(state.sz[0]), sobj_slit[i].yedg_sky[0:state.sz[0]-1,1],$
        color=1
      oplot, findgen(state.sz[0]), sobj_slit[i].yedg_sky[0:state.sz[0]-1,0],$
        color=5
  endfor

  x_setobjgui_resetwindow
  state.newrefresh=1

  return
end

;;;;;;;;;
; QUERY SLIT

pro x_setobjgui_query

common x_setobjgui_state

  ;; Find the slit
  mnbot = min(abs(state.coord[1]-sobj_slit.yedg_sky[state.coord[0],0]),ibot)
  mntop = min(abs(state.coord[1]-sobj_slit.yedg_sky[state.coord[0],1]),itop)
  if mnbot LT mntop then gd = ibot else gd = itop

  ;; Update
  tmp_string = string(sobj_slit[gd].id, 0.0,$
                      format = '(i5,3x,g14.7)' )
  widget_control, state.wavesig_bar_id, set_value = tmp_string

  print, 'Slit ID = ', strtrim(sobj_slit[gd].id,2)

  return
end
  


;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;--------------------------------------------------------------------
;    x_setobjgui main program.  needs to be last in order to compile.
;---------------------------------------------------------------------

pro x_setobjgui, image, objstr, slitstr, $
                 min = minimum, $
                 max = maximum, $
                 autoscale = autoscale,  $
                 linear = linear, $
                 log = log, $
                 histeq = histeq, $
                 noblock = noblock, $
                 align = align, $
                 stretch = stretch, $
                 header = header, $
                 wvimg=wvimg

common x_setobjgui_state
common x_setobjgui_images

if (not(keyword_set(noblock))) then block = 1 else block = 0

newimage = 0

;if ( (n_params() LT 3) AND (xregistered('x_setobjgui', /noshow))) then begin
if (n_params() LT 3) then begin
    print, 'USAGE: x_setobjgui, array_name OR fitsfile, objstr, slitstr '
    print, '         [,min = min_value] [,max=max_value] '
    print, '         [,/linear] [,/log] [,/histeq] [,/block]'
    print, '         [,/align] [,/stretch] [,header=header]'
    return
endif

; Set structures
  sobj_obj = objstr
  sobj_slit = slitstr

  x_setobjgui_shutdown, windowid

if (!d.name NE 'X' AND !d.name NE 'WIN' AND !d.name NE 'MAC') then begin
    print, 'x_setobjgui: Graphics device must be set to X, WIN, or MAC for this GUI to work.'
    retall
endif

; Before starting up x_setobjgui, get the user's external window id.  We can't
; use the x_setobjgui_getwindow routine yet because we haven't run x_setobjgui
; startup.  A subtle issue: x_setobjgui_resetwindow won't work the first time
; through because xmanager doesn't get called until the end of this
; routine.  So we have to deal with the external window explicitly in
; this routine.
if (not (xregistered('x_setobjgui', /noshow))) then begin
   userwindow = !d.window
   x_setobjgui_startup
   align = 0B     ; align, stretch keywords make no sense if we are
   stretch = 0B   ; just starting up. 
endif


; If image is a filename, read in the file
if ( (n_params() NE 0) AND (size(image, /tname) EQ 'STRING')) then begin
    ifexists = findfile(image, count=count)
    if (count EQ 0) then begin
        print, 'ERROR: File not found!'
    endif else begin
        x_setobjgui_readfits, fitsfilename=image, newimage=newimage
    endelse
endif



; Check for existence of array
if ( (n_params() NE 0) AND (size(image, /tname) NE 'STRING') AND $
   (size(image, /tname) EQ 'UNDEFINED')) then begin
    print, 'ERROR: Data array does not exist!'
endif

; If user has passed x_setobjgui a data array, read it into main_image.
if ( (n_params() NE 0) AND (size(image, /tname) NE 'STRING') AND $
   (size(image, /tname) NE 'UNDEFINED')) then begin
; Make sure it's a 2-d array
    if ( (size(image))[0] NE 2) then begin
        print, 'ERROR: Input data must be a 2-d array!'    
    endif else begin
        main_image = image
        newimage = 1
        state.imagename = ''
        state.title_extras = ''
        x_setobjgui_setheader, header
    endelse
endif


; WAVE IMAGE
if keyword_set( WVIMG ) then begin
    ; Update flag
    if state.wavesig_flg MOD 2 EQ 0 then state.wavesig_flg = state.wavesig_flg + 1
    wave_image = x_readimg(wvimg, /dscale)
endif

;   Define default startup image 
if (n_elements(main_image) LE 1) then begin
    main_image = cos(((findgen(500)- 250)*2) # ((findgen(500)-250)))
    imagename = ''
    newimage = 1
    x_setobjgui_setheader, ''
endif

  state.sz = size(main_image, /dimensions)

if (newimage EQ 1) then begin  
; skip this part if new image is invalid or if user selected 'cancel'
; in dialog box
    x_setobjgui_settitle
    x_setobjgui_getstats, align=align
    
    delvarx, display_image

    if n_elements(minimum) GT 0 then begin
        state.min_value = minimum
    endif
    
    if n_elements(maximum) GT 0 then begin 
        state.max_value = maximum
    endif
    
    if state.min_value GE state.max_value then begin
        state.min_value = state.max_value - 1.
    endif
    
    if (keyword_set(linear)) then state.scaling = 0
    if (keyword_set(log))    then state.scaling = 1
    if (keyword_set(histeq)) then state.scaling = 2
    
; Only perform autoscale if current stretch invalid or stretch keyword
; not set
    IF (state.min_value EQ state.max_value) OR $
      (keyword_set(stretch) EQ 0) THEN BEGIN 

       if (keyword_set(autoscale) OR $
           ((state.default_autoscale EQ 1) AND (n_elements(minimum) EQ 0) $
            AND (n_elements(maximum) EQ 0)) ) $
         then x_setobjgui_autoscale
    ENDIF 
    x_setobjgui_set_minmax
    
    IF NOT keyword_set(align) THEN BEGIN 
       state.zoom_level = 0
       state.zoom_factor = 1.0
    ENDIF 

    x_setobjgui_displayall
    
    x_setobjgui_resetwindow
endif



; Register the widget with xmanager if it's not already registered
if (not(xregistered('x_setobjgui', /noshow))) then begin
    nb = abs(block - 1)
    xmanager, 'x_setobjgui', state.base_id, no_block = nb, cleanup = 'x_setobjgui_shutdown'
    wset, userwindow
    ; if blocking mode is set, then when the procedure reaches this
    ; line x_setobjgui has already been terminated.  If non-blocking, then
    ; the procedure continues below.  If blocking, then the state
    ; structure doesn't exist any more so don't set active window.
    if (block EQ 0) then state.active_window_id = userwindow
endif

; RESET the OBJ Structure
  gdobj = where(sobj_obj.flg_anly NE 0)
  objstr = temporary(sobj_obj[gdobj])
  slitstr = temporary(sobj_slit)
  if keyword_set(STATE) then delvarx, state

end





