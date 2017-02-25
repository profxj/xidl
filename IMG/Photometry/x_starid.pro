;+ 
; NAME:
; x_starid   
;    Version 1.1
;
; PURPOSE:
;    Allows the user to interactively identify standard stars in a
;    field using lists and the like.  This is a rather quick way of
;    doing the calibrations.
;
; CALLING SEQUENCE:
;   x_starid, img, date, ccd, tel, LST_PATH=, OBJ=, XSIZE=, YSIZE=
;
; INPUTS:
;   img        - Image(s) for Masking
;   date       - EPOCH of the obs (decimal years)
;   ccd        - Name of the CCD  (Tek5, SITe1, LRISR, SITe3)
;   tel        - Name of telescope
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   OBJ       - Object names
;   XSIZE     - Size of gui in screen x-pixels (default = 800)
;   YSIZE     - Size of gui in screen y-pixels (default = 800)
;   OUT_PATH  - Output name for id files (default: 'Photo/')
;   LST_PATH  - Path name for Lists [default: $XIDL_DIR/IMG/Photometry/Lists/]
;   INORIENT  - Number in range [-4,4] to change orientation of 
;                  standards plotted on image.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_starid, img, 2004.2, 'SITe1', 'LCO-40'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Aug-2001 Written by JXP
;   04-Jan-2002 Updated and revised by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Common
;;;;

pro x_starid_initcommon

;

  ; COLOR
  common xcommon_color, r_vector, g_vector, b_vector
  xicmm_colors

  ; IMAGES
  common x_starid_images, $
    main_image, $
    img_size, $ 
    display_image, $
    tv_image, $
    header

end

;;;;
; Events
;;;;

pro x_starid_event, ev

common xcommon_color
common x_starid_images

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
;      'ERRORB' : widget_control, state.error_msg_id, set_value=''
      'CONTRAST': begin
          state.contrast = ev.value
          ximgd_stretchct, state
          x_starid_PlotImg, state
      end
      'BRIGHT': begin
          state.brightness = ev.value
          ximgd_stretchct, state
          x_starid_PlotImg, state
      end
      'FLIP': x_starid_Flip, state
      'SVLIST': x_starid_SvList, state
      'REDO': begin
          print, 'Be sure to choose a new list!'
          x_starid_Reset, state
          x_starid_PlotImg, state
      end
      'LISTS': begin
          state.strmode = 1
          state.nlist = ev.index
          x_starid_PrsList, state
          if state.flg_strs NE 0 then begin
              state.flg_strs = 0
              x_starid_PlotImg, state
          endif
      end
      'DRAW' : begin
          case ev.type of
              0 : begin ; Button press
                  state.press = ev.press
                  case ev.press of
                      1 : begin     ; Establish a star
                          x_starid_SetStr, state
                          if state.strmode NE 0 then x_starid_PlotImg, state
                      end 
                      2 : begin    ; Delete star
                      end
                      4 : x_starid_SetZoom, state, 1  ; Set Zoom region
                  endcase
              end
              1 : begin ; Button release
                  if( state.press EQ 4) then begin
                      x_starid_SetZoom, state, 2 
                      x_starid_UpdDisplay, state
                      x_starid_PlotImg, state
                  endif
              end
              2 : begin ; Motion event
                  state.tv.xcurs = ev.x
                  state.tv.ycurs = ev.y
              end
          endcase
      end
      'CENTER' : begin
          x_starid_CenterAll, state
          x_starid_PlotImg, state
      end
      'FULLIMG' : begin
          state.zoomreg = [0,0,img_size[0]-1,img_size[1]-1]
          x_setgridxy, state, state.zoomreg
          x_starid_UpdDisplay, state
          x_starid_PlotImg, state
      end
      'NEXT' : begin
          if state.strmode NE 0 then x_starid_Output, state
          if(state.curimg EQ state.nimg-1) then begin 
              widget_control, ev.top, /destroy 
              return
          endif else begin      ; Next image
              state.curimg = state.curimg + 1
              if state.flg_flip EQ 1 then begin
                  state.orient = -1*state.orient
                  state.flg_flip = 0
              endif
              widget_control, state.list_id, set_list_select=-1
              x_starid_ReadImg, state
              x_starid_PlotImg, state
          endelse
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_starid_PlotImg, state
  
common x_starid_images
common xcommon_color

; Plot Data

  widget_control, state.draw_id, get_value=wind
  wset, wind

  tv_image[0:state.tv.gridsize[0]-1, 0:state.tv.gridsize[1]-1] = $
    congrid(display_image, state.tv.gridsize[0], state.tv.gridsize[1])

; Edges
  
  if state.tv.gridsize[0] LT state.tv.winsize[0] then $
    tv_image[state.tv.gridsize[0]:state.tv.winsize[0]-1,*] = 0
  if state.tv.gridsize[1] LT state.tv.winsize[1] then $
    tv_image[*,state.tv.gridsize[1]:state.tv.winsize[1]-1,*] = 0

  tv, tv_image

; Points

  if state.flg_strs EQ 1 then x_starid_plotpts, state

end

;----------------------------------------------------------------------

pro x_starid_plotpts, state

  xrange=[state.tv.xymnx[0], state.tv.xymnx[2]]
  yrange=[state.tv.xymnx[1], state.tv.xymnx[3]]
  
  plot, extrac(state.x_str, 0, state.nstrs), $
    extrac(state.y_str, 0, state.nstrs), $
    position=[0.,0.,1., 1.], $
    xrange=xrange, yrange=yrange, xstyle=5, ystyle=5, /noerase, $
    psym=1, color=1

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_starid_Reset, state

; Plotting

  state.brightness = 0.5
  widget_control, state.drag_contr_id, set_value=state.brightness
  state.contrast = 0.5
  widget_control, state.drag_brght_id, set_value=state.contrast
  ximgd_stretchct, state
  state.flg_strs = 0	
  state.strmode = 0	
  widget_control, state.list_id, set_list_select=-1

end

;;;;;;;;;;;;;;;;;;;;
;  ReadImage
;;;;;;;;;;;;;;;;;;;;

pro x_starid_ReadImg, state

common x_starid_images

  widget_control, /hourglass

; Read Fits
  delvarx, main_image
  main_image = xmrdfits( state.img[state.curimg], 0, header, /fscale, /silent)
  img_size = size(main_image, /dimensions)


; Set tv size

  state.zoomreg = [0, 0, img_size[0]-1, img_size[1]-1]
  x_setgridxy, state, state.zoomreg

; Set min max

  mx = max(main_image, min=mn)
  state.imgmin = min(main_image, max=fmax)
  state.imgmax = fmax
  med = median(main_image, /even)
  state.pltmax = (med + 200) < mx
  state.pltmin = (med - 100) > mn

  state.tv.svxymnx = state.tv.xymnx

; Display Image

  delvarx, display_image
  display_image = bytscl(main_image, min=state.pltmin, max=state.pltmax, /nan, $
                        top=state.ncolors-1) + 8B

; Name
  widget_control, state.name_id, set_value=strmid(state.img[state.curimg],0,20)
  widget_control, state.obj_id, set_value=strmid(state.obj[state.curimg],0,20)

; Reset stuff

  state.strmode = 0
  state.flg_strs = 0
  state.nstrs = 0

end

;;;;;;;;;
;  Flip
;;;;;;;;;

pro x_starid_Flip, state

common x_starid_images

;  Flips the display image and clears out any pre-existing stars  

  widget_control, /hourglass

  main_image = rotate(temporary(main_image), 5) ; Flip x

  state.zoomreg = [0,0,img_size[0]-1,img_size[1]-1]
  x_setgridxy, state, state.zoomreg
  if state.flg_strs NE 0 then begin
      state.flg_strs = 0
      x_starid_PlotImg, state
  endif
  x_starid_UpdDisplay, state
  x_starid_PlotImg, state

; Reset orientation

  state.orient = -1*state.orient
  state.flg_flip = abs(state.flg_flip - 1)

return
end

;;;;;;;;;
;  Update Display
;;;;;;;;;

pro x_starid_UpdDisplay, state

common x_starid_images

  widget_control, /hourglass

; Reset Display

  delvarx, display_image
  display_image = bytscl(main_image[state.zoomreg[0]:state.zoomreg[2], $
                                    state.zoomreg[1]:state.zoomreg[3]], $
                         min=state.pltmin, max=state.pltmax, $
                         /nan, top=state.ncolors-1) + 8B
end

;;;;;;;;;
;  Set Zoom 
;;;;;;;;;

pro x_starid_SetZoom, state, flg

common x_starid_images

  if flg EQ 1 then begin
      state.tmpreg[0] = x_tvx(state.tv, /intg)
      state.tmpreg[1] = x_tvy(state.tv, /intg)
  endif else begin
      state.tmpreg[2] = x_tvx(state.tv, /intg)
      state.tmpreg[3] = x_tvy(state.tv, /intg)

      if (state.tmpreg[0] NE state.tmpreg[2] AND $
          state.tmpreg[1] NE state.tmpreg[3]) then begin
          state.zoomreg[0] = state.tmpreg[0] < state.tmpreg[2]
          state.zoomreg[2] = state.tmpreg[0] > state.tmpreg[2]
          state.zoomreg[1] = state.tmpreg[1] < state.tmpreg[3]
          state.zoomreg[3] = state.tmpreg[1] > state.tmpreg[3]

          state.zoomreg[0] = state.zoomreg[0] > 0
          state.zoomreg[2] = state.zoomreg[2] > 0
          state.zoomreg[1] = state.zoomreg[1] < img_size[0]-1
          state.zoomreg[3] = state.zoomreg[3] < img_size[1]-1

          x_setgridxy, state, state.zoomreg
      endif
  endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;
;  Set List
;;;;;;;;;
   
pro x_starid_PrsList, state

; Read the List 
  strfil = strjoin([state.lpath, state.lists[state.nlist]])
  readcol, strfil, FORMAT='A,A,A', nam, ra, dec
  state.nstrs = n_elements(nam)
  state.nam_str[0:state.nstrs-1] = nam
  
; Convert RA, DEC

  for i=0L,state.nstrs-1 do begin
      x_radec, ra[i], dec[i], rad, decd
      ; Precess
      precess, rad, decd, 2000.0, state.date
      state.ra_str[i] = rad
      state.dec_str[i] = decd
      state.cen_str[i] = 0
  endfor

end

;;;;;;;;;
;  Set Stars
;;;;;;;;;

pro x_starid_SetStr, state

common x_starid_images

  widget_control, /hourglass

  case state.strmode of 
      0 : begin
          print, 'You need to chose a list first!'
      end
      1 : begin    ; Find first star!
          x = x_tvx(state.tv)
          y = x_tvy(state.tv)

          cntrd, main_image, x, y, xcen, ycen, 8.0
          if xcen EQ -1 then begin
              print, 'You missed the star, maybe zoom in'
              state.nstrs = state.nstrs - 1
              return
          endif 
          state.flg_strs = 1
          ; First star position
          state.x_str[0] = xcen
          state.y_str[0] = ycen
          for i=1,state.nstrs-1 do begin
              ; Calculate delta arcsec
              dxa = (state.ra_str[i]-state.ra_str[0])*3600.* $
                cos(state.dec_str[i]*!pi/180.)
              dya = (state.dec_str[i]-state.dec_str[0])*3600.
              ; Pixels
              dxp = dxa/state.arcpix
              dyp = dya/state.arcpix
              ; Orient


              case state.orient of 
                  1 : begin  ; standard orientation (N up, E left)
                      state.x_str[i] = state.x_str[0] - dxp
                      state.y_str[i] = state.y_str[0] + dyp
                  end
                  2 : begin  
                      state.x_str[i] = state.x_str[0] + dyp
                      state.y_str[i] = state.y_str[0] + dxp
                  end
                  -2 : begin  ; e.g. SITe1 w/o flip
                      state.x_str[i] = state.x_str[0] - dyp
                      state.y_str[i] = state.y_str[0] + dxp
                  end
                  -3 : begin  ; e.g. Tek5 w/o flip
                      state.x_str[i] = state.x_str[0] - dxp
                      state.y_str[i] = state.y_str[0] - dyp
                  end
                  3 : begin  ; e.g. Tek5 w/ flip, LRISR
                      state.x_str[i] = state.x_str[0] + dxp
                      state.y_str[i] = state.y_str[0] - dyp
                   end
                  4 : begin  
                      state.x_str[i] = state.x_str[0] - dyp
                      state.y_str[i] = state.y_str[0] - dxp
                  end
                  -4 : begin  ; e.g. SITe1 w/o flip
                      state.x_str[i] = state.x_str[0] + dyp
                      state.y_str[i] = state.y_str[0] - dxp
                  end


                  else : print, 'Oops!'
              endcase
          endfor
          state.strmode = 2
      end
      2 : begin    ; Delete nearest and accept value
          xpx = x_tvx(state.tv)
          ypx = x_tvy(state.tv)

          dist = sqrt( (xpx-state.x_str)^2 + (ypx-state.y_str)^2 )
          mdist = min(dist, jmin)

          state.x_str[jmin] = xpx
          state.y_str[jmin] = ypx
      end

      else :
  endcase
end
  
;;;;;;;;;
;  Centroid all the stars on the image
;;;;;;;;;

pro x_starid_CenterAll, state

common x_starid_images


  for i=0,state.nstrs-1 do begin
      ; Check that it is on the image!
      if (state.x_str[i] GT 10. AND $
          state.x_str[i] LT float(img_size[0])-10. AND $
          state.y_str[i] GT 10. AND $
          state.y_str[i] LT float(img_size[1])-10.) then begin

          cntrd, main_image, state.x_str[i], state.y_str[i], xcen, ycen, $
            7.0/state.arcpix
          
          if(xcen NE -1) then begin
              state.x_str[i] = xcen
              state.y_str[i] = ycen
              state.cen_str[i] = 1
          endif else print, state.nam_str[i], ' did not centroid!'
      endif
  endfor
return
end
        


;;;;
; Output 
;;;;

pro x_starid_Output, state

common x_starid_images

  widget_control, /hourglass

; Name
  lastslsh = strpos(state.img[state.curimg],'/',/reverse_search)
  if lastslsh NE -1 then nopth = strmid(state.img[state.curimg], lastslsh+1) $
  else nopth = state.img[state.curimg]

; Take off the xx_
;  undpos = strpos(nopth,'_') 
;  if undpos NE -1 then nopth = strmid(nopth,undpos+1)

; Take off the fits
  lenno = strlen(nopth) 
  nopth = strmid(nopth,0,lenno-5)

; Add the Directory and header info
  outfil = strjoin( [state.opath, 'std_', nopth, '.dat'] ) 


; Write
  close, /all
  openw, 1, outfil

  printf, 1, state.img[state.curimg]
  printf, 1, strjoin([state.lpath, state.lists[state.nlist]])

  for i=0,state.nstrs-1 do begin
      ; Check that it is on the image!
      if (state.x_str[i] GT 10. AND $
          state.x_str[i] LT float(img_size[0])-10. AND $
          state.y_str[i] GT 10. AND $
          state.y_str[i] LT float(img_size[1])-10.) then begin
          if state.flg_flip EQ 0 then $
            printf, 1, FORMAT='(a15,f,f)',state.nam_str[i], $
            state.x_str[i],state.y_str[i] $
          else $
            printf, 1, FORMAT='(a15,f,f)', state.nam_str[i], $
            float(img_size[0]-1)-state.x_str[i], state.y_str[i]
      endif

  endfor

  close, 1

 
return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_starid_SvList, state 

;;
; Saves the current positions of stars into a new list
;  It takes the RA and DEC from the first star and calculates new
;  ones from the current positions for the rest.
;;

  ; Grab the good stars
  if state.cen_str[0] EQ 0 then begin
	print, 'x_starid_SvList: First star needs to be centered!'
	return
  endif	
  gdstrs = where(state.cen_str NE 0, ngd)
  if ngd EQ 0 then begin
	print, 'x_starid_SvList: No stars worth keeping!'
	return
  endif	
  ; Grab the new filename 
  filnm = x_guistring('Enter the list name (full path)')

  ; Convert positions of stars to RA,DEC
  close, 1
  openw, 1, filnm
  ; First one first
  x_radec, ra, dec, state.ra_str[0], state.dec_str[0], /flip
  printf, 1, FORMAT='(A15,1x,A10,1x,A10)', state.nam_str[0], ra, dec
  for j=1L,ngd-1 do begin
     ; Deal with orientation!
     i = gdstrs[j]
     case state.orient of 
      1 : begin  ; standard orientation (N up, E left)
       dxp = state.x_str[0] - state.x_str[i]
       dyp = state.y_str[i] - state.y_str[0]
      end
      2 : begin  ; e.g. SITe1 w/o flip
       dyp = state.x_str[i] - state.x_str[0]
       dxp = state.y_str[i] - state.y_str[0]
      end
      -2 : begin  ; e.g. SITe1 w/o flip
       dyp = state.x_str[0] - state.x_str[i]
       dxp = state.y_str[i] - state.y_str[0]
      end
      -3 : begin  ; e.g. Tek5 w/o flip
       dxp = state.x_str[0] - state.x_str[i]
       dyp = state.y_str[0] - state.y_str[i]
      end
      3 : begin  ; e.g. Tek5 w/ flip, LRIS
       dxp = state.x_str[i] - state.x_str[0]
       dyp = state.y_str[0] - state.y_str[i]
      end
       else : print, 'Oops!'
     endcase
    ; Convert to decimal ra, dec
    dxa = dxp*state.arcpix
    dya = dyp*state.arcpix
    decd_ra  = dxa/cos(state.dec_str[0]*!pi/180.)/3600. + state.ra_str[0]
    decd_dec = dya/3600. + state.dec_str[0]
    x_radec, ra, dec, decd_ra, decd_dec, /flip
    printf, 1, FORMAT='(A15,2x,A10,2x,A10)', state.nam_str[i], ra, dec
  endfor
  close, 1

  ; Reload the lists
  a = findfile(strjoin([state.lpath,'*.lst']), count=nlist)
  state.nlist = nlist
  if state.nlist NE 0 then $
    for i = 0,state.nlist-1 do state.lists[i] = fileandpath(a[i])
  widget_control, state.list_id, set_value=state.lists[0:state.nlist-1]
  widget_control, state.list_id, set_list_select=-1

return

end

; 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_starid, img, date, ccd, tel, LST_PATH=lst_path, XSIZE=xsize, YSIZE=ysize,$
              OUT_PATH=out_path, OBJ=obj, INORIENT=inorient

common xcommon_color
common x_starid_images

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'x_starid, img, date, ccd, tel, LST_PATH=, XSIZE=, YSIZE=, OBJ=,'
    print, '      OUT_PATH= (v1.1)'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then    xsize = 1000
  if not keyword_set( YSIZE ) then    ysize = 1000
  if not keyword_set( OBJ ) then  begin
	  obj = strarr(n_elements(img))
	  obj = '?'
  endif
  if not keyword_set( LST_PATH ) then $
    lst_path = getenv('XIDL_DIR')+'/IMG/Photometry/Lists/'
  if not keyword_set( OUT_PATH ) then begin
      out_path = 'Photo/'
      a = findfile(out_path+'..', count=count)
      if count EQ 0 then file_mkdir, out_path
  endif
  

;    STATE
  ; tv structure
  tmp = {tvstruct}
  tmp.pos = [0., 0., 1., 1.]

  state = {             $
            img: img, $                  ; Important stuff
            obj: obj, $                  ; Important stuff
            nimg: 0, $
            curimg: 0, $
            orient: 0, $ 
            flg_flip: 0, $
            arcpix: 0.d, $
            date: date, $
            nlist: 0, $                  ; Lists
            lists: strarr(500), $
            lpath: lst_path, $
            opath: out_path, $
            flg_strs: 0, $               ; Stars
            strmode: 0, $
            nstrs: 0, $
            nam_str: strarr(1000), $
            ra_str: dblarr(1000), $
            dec_str: dblarr(1000), $
            x_str: fltarr(1000), $
            y_str: fltarr(1000), $
            cen_str: bytarr(1000), $       ; Centroid flag
            imgmin: 0., $                ; Image stuff
            imgmax: 0., $
            zoomreg: lonarr(4), $         
            tmpreg: lonarr(4), $         
            tv: tmp, $              ; TV structure
            ncolors: 0L, $
            brightness: 0.0, $
            contrast: 0.0, $
            pltmax: 0.0, $
            pltmin: 0.0, $
            press: 0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            drag_contr_id: 0L, $
            drag_brght_id: 0L, $
            draw_id: 0L, $
            dx_text_id: 0L, $
            dy_text_id: 0L, $
            list_id: 0L, $
            list_butt_id: 0L, $
            name_id: 0L, $
            obj_id: 0L, $
            dithsv_id: 0L, $
            nodith_id: 0L, $
            error_msg_id: 0L, $
            help_text_id: 0L $
          }

;  Images

  state.nimg = n_elements(img)

;    WIDGET
  base = WIDGET_BASE( title = 'x_starid: ID Stars', /column, $
                      xoffset=100,yoffset=300)
  state.base_id = base
  
;      Toolbar
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
;        Version + Name
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='x_starid', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.1)', /align_center)
  state.name_id = WIDGET_LABEL(labelbase, value=strmid(img[0],0,20), $
                               /align_center)
  state.obj_id = WIDGET_LABEL(labelbase, value=strtrim(strmid(obj[0],0,20),2),$
                              /align_center)

;;;;;;;;;;;;;;
;        List Buttons

;  lbstr = ['New List', 'Load List', 'Save List']
;  state.list_butt_id = cw_bgroup(toolbar, lbstr, uvalue='LSBUTT', $
;                                /COLUMN, /EXCLUSIVE, /NO_RELEASE)

;;;;;;;;;;;;;;
;        Lists
  a = findfile(strjoin([lst_path,'*.lst']), count=nlist)
  state.nlist = nlist

  if state.nlist NE 0 then $
    for i = 0,state.nlist-1 do state.lists[i] = fileandpath(a[i])

  listbase = widget_base(toolbar, /column, /align_center)
;  listlabel = widget_label(listbase, value='LISTS', /align_center)
  state.list_id = WIDGET_LIST(listbase, VALUE=state.lists[0:state.nlist-1], $
                              uvalue='LISTS', ysize = 6)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             '1. Choose a Standard field (list)', $ 
             '2. Identify the first star (LMB)', $
             '3. Modify the other stars (LMB)', $ 
             '4. Center (CENTER)', $ 
             '5. Proceed to the next (SAVE)', $
             'Save new positions: (SVLIST)', $
             'Start over on this field: (REDO)', $
             'Flip in x:  (FLIP)', $
             'Zoom Region: RMB+drag']
  help_text_id = widget_text(toolbar, value=strhelp, xsize=40, ysize=4,$
                             /scroll)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Drawing
  drawbase = WIDGET_BASE( state.base_id, /row, /base_align_center,/align_center)

  state.tv.winsize[0] = xsize
  state.tv.winsize[1] = ysize
  state.draw_id = widget_draw(drawbase, xsize=state.tv.winsize[0], $
                              ysize=state.tv.winsize[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Image Stuff
;  imgbase = WIDGET_BASE( drawbase, /column, /base_align_center,/align_center)

;        Contrast+Brightness
  conbase = WIDGET_BASE( toolbar, /column, /base_align_center,/align_center)
  contrtext = widget_label(conbase, value='Contr/Bright', /align_center)
  state.drag_contr_id = CW_FSLIDER(conbase, TITLE='Contrast', $
                                   VALUE=0.5, uvalue='CONTRAST', $
                                   maximum=1.0, minimum=0.0, ysize=15, $
                                   /frame, /suppress_value)
  state.drag_brght_id = CW_FSLIDER(conbase, VALUE=0.5, uvalue='BRIGHT', $
                                   maximum=1.0, minimum=0.0, ysize=15, $
                                   TITLE='Brightness', /frame, /suppress_value)
;      Buttons1
  butbase = widget_base(toolbar, /column, /align_center)
  flipx = WIDGET_BUTTON(butbase, value='FLIP-X',uvalue='FLIP')
  centerall = WIDGET_BUTTON(butbase, value='CENTER',uvalue='CENTER')
  full_img = WIDGET_BUTTON(butbase, value='FULLIMG',uvalue='FULLIMG')
;      Buttons2
  butbase2 = widget_base(toolbar, /column, /align_center)
  svlist = WIDGET_BUTTON(butbase2, value='SVLIST',uvalue='SVLIST')
  redo = WIDGET_BUTTON(butbase2, value='REDO',uvalue='REDO')
  done = WIDGET_BUTTON(butbase2, value='NEXT',uvalue='NEXT')
  
; Realize
  WIDGET_CONTROL, base, /realize

; Colors 
  device, get_decomposed=val_decomp
  device, decompose=0
  loadct, 0, /silent
  if (!d.table_size LT 12) then begin
      message, 'Too few colors available for color table'
      stop
  endif
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)

  ximgd_getct, state, 0, /CLR

; INVERT
  r_vector = reverse(r_vector)
  g_vector = reverse(g_vector)
  b_vector = reverse(b_vector)
  ximgd_stretchct, state
  x_starid_Reset, state

; TV
  tv_image = bytarr(xsize, ysize)

  
; Set Orientation and arcpix

  x_ccdinf, ccd, tel, arc, orient
  state.arcpix = arc
  state.orient = orient

if keyword_set(INORIENT) then state.orient = inorient

; Update
  x_starid_ReadImg, state

;  x_starid_UpdateImg, base
  x_starid_PlotImg, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'x_starid', base

; Finish
  delvarx, main_image, tv_image, display_image, header
  device, decomposed=val_decomp

  return
end

