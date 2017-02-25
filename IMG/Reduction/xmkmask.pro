;+ 
; NAME:
; xmkmask   
;   Version 1.1
;
; PURPOSE:
;    Allows the user to interactively choose pixels to mask in an
;    image.  The resulting images are saved as gzipped fits images.
;
; CALLING SEQUENCE:
;   
;   xmkmask, img, IMSK=, OUTDIR=, XSIZE=, YSIZE=, /SKYMSK
;
; INPUTS:
;   img        - Image(s) for Masking
;
; RETURNS:
;
; OUTPUTS:
;   mask       -  Creates a fits table of pixel values to mask
;
; OPTIONAL KEYWORDS:
;   outdir     - Output directory: Default = './Masks/'
;   IMSK       - Initial mask to include; Convenient for dealing with
;                  bad columns and the like in the CCD
;   XSIZE      - Size of gui in screen x-pixels (default = 700)
;   YSIZE      - Size of gui in screen y-pixels (default = 700)
;   SKYMSK     - Create sky mask (default into Masks/Sky/)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xmkmask,'OV/ccd001.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-July-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
; Common
;;;;

pro xmkmask_initcommon

;

common xmkmask_color, r_vector, g_vector, b_vector
common xmkmask_images, $
  main_image, $
  img_size, $ 
  display_image, $
  mask_img, $
  init_msk, $
  cmmn_regions, $
  cmmn_nreg, $
  cmmn_regtype, $
  header
	
end

;;;;
; Events
;;;;

pro xmkmask_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
;      'ERRORB' : widget_control, state.error_msg_id, set_value=''
      'CONTRAST': begin
          state.contrast = ev.value
          xmkmask_stretchct, state
          xmkmask_PlotImg, state
      end
      'BRIGHT': begin
          state.brightness = ev.value
          xmkmask_stretchct, state
          xmkmask_PlotImg, state
      end
      'REGLIST' : if( state.flgline NE 0) then print, 'Still in RAY' else $
            state.regtyp = ev.index
      'DRAW' : begin
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Start the Region
                          if state.flgdith EQ -1 then begin
                              print, 'You cant modify regions until dithering is done'
                              replot = 0
                              break
                          endif
                          xmkmask_SetReg, state, 0, replot=replot
                          state.press = 1
                      end 
                      4 : begin    ; Delete/Undelete a region
                          xmkmask_DelReg, state
                          if state.flgdith EQ -1 then xmkmask_Dither, state
                          replot = 1
                      end
                      else :
                  endcase
                  if replot EQ 1 then xmkmask_PlotImg, state
              end
              1 : begin ; Button release
                  if( state.regtyp NE 2 AND state.press EQ 1) then begin
                      xmkmask_SetReg, state, 1, replot=replot
                      if replot EQ 1 then xmkmask_PlotImg, state
                      state.press = 0
                  endif
              end
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
              end
          endcase
      end
      'DXTEXT' : begin
          state.dithx = ev.value
          xmkmask_Dither, state
          xmkmask_PlotImg, state
      end
      'DYTEXT' : begin
          state.dithy = ev.value
          xmkmask_Dither, state
          xmkmask_PlotImg, state
      end
      'SVDITH' : begin
          xmkmask_Dither, state, /SVDITH
          widget_control, state.dx_text_id, sensitive=0
          widget_control, state.dy_text_id, sensitive=0
          widget_control, state.dithsv_id, sensitive=0
          widget_control, state.nodith_id, sensitive=0
      end
      'NODITH' : begin
          xmkmask_Dither, state, /NODITH
          widget_control, state.dx_text_id, sensitive=0
          widget_control, state.dy_text_id, sensitive=0
          widget_control, state.dithsv_id, sensitive=0
          widget_control, state.nodith_id, sensitive=0
      end
      'NEXT' : begin
          if state.flgdith EQ -1 then begin
              print, 'You cant move on until settling the dither image'
              break
          endif
          xmkmask_UpdDisplay, state, /MASK
          xmkmask_Output, state
          if(state.curimg EQ state.nimg-1) then begin 
              cmmn_regions = state.reg
              cmmn_regtype = state.reg_typ
              cmmn_nreg = state.nreg
              widget_control, ev.top, /destroy 
              return
          endif else begin      ; Next image
              state.curimg = state.curimg + 1
              xmkmask_ReadImg, state
              xmkmask_PlotImg, state
          endelse
      end
      'DONE' : begin
          cmmn_regions = state.reg
          cmmn_regtype = state.reg_typ
          cmmn_nreg = state.nreg
          widget_control, ev.top, /destroy 
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xmkmask_PlotImg, state
  
common xmkmask_images

; Plot Data

  widget_control, state.draw_id, get_value=wind
  wset, wind

  tv, congrid(display_image, state.gridsize[0], state.gridsize[1])

end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro xmkmask_Reset, state

; Plotting

  state.brightness = 0.5
  widget_control, state.drag_contr_id, set_value=state.brightness
  state.contrast = 0.5
  widget_control, state.drag_brght_id, set_value=state.contrast
  xmkmask_stretchct, state


end

;;;;;;;;;;;;;;;;;;;;
;  ReadImage
;;;;;;;;;;;;;;;;;;;;

pro xmkmask_ReadImg, state

common xmkmask_images

  widget_control, /hourglass

; Read Fits
  delvarx, main_image
  main_image = xmrdfits( state.img[state.curimg], 0, header, /fscale, /silent)
  sz = size(main_image)
  img_size = lonarr(2)
  img_size[0] = sz[1]
  img_size[1] = sz[2]

; Set tv size

  if float(state.winsize[0])/float(sz[1]) LT $
    float(state.winsize[1])/float(sz[2]) then begin
      state.gridsize[0] = state.winsize[0]
      state.gridsize[1] = $
        round( float(sz[2])*float(state.winsize[0])/float(sz[1]))
      state.xymnx[0] = 0.
      state.xymnx[1] = 0.
      state.xymnx[2] = float(img_size[0]-1)
      state.xymnx[3] = float(state.winsize[1])*float(img_size[1]-1)/ $
        float(state.gridsize[1])
  endif else begin
      state.gridsize[1] = state.winsize[1]
      state.gridsize[0] = $
        round( float(sz[1])*float(state.winsize[1])/float(sz[2]))
      state.xymnx[0] = 0.
      state.xymnx[1] = 0.
      state.xymnx[2] = float(state.winsize[0])*float(img_size[0]-1)/ $
        float(state.gridsize[0])
      state.xymnx[3] = float(img_size[1]-1)
  endelse

; Set min max

  state.imgmin = min(main_image, max=fmax)
  state.imgmax = fmax
  med = median(main_image, /even)
  state.pltmax = med + 0.2*abs(med)
  state.pltmin = med - 0.2*abs(med)

  state.svxymnx = state.xymnx

; Display Image

  delvarx, display_image
  display_image = bytscl(main_image, min=state.pltmin, max=state.pltmax)

; Initial Mask

  if (state.flgimsk EQ 1 ) then begin
      if( state.curimg EQ 0) then begin ; Read in the Init Mask
          imgmsk = mrdfits( state.imsk, /silent )
                                ; Check size
          szi = size( imgmsk )
          if(sz[1] NE szi[1] OR sz[2] NE szi[2] ) then begin
              print, 'Init mask is wrong size!'
              state.flgimsk = 0
              return
          endif 
          init_msk = where( imgmsk NE 0 )
      endif 
      display_image[init_msk] = 0
  endif

; Name
  widget_control, state.name_id, set_value=strmid(state.img[state.curimg],0,18)

; Dither Image
  if(state.curimg NE 0) then begin
      state.flgdith = -1
      widget_control, state.dx_text_id, sensitive=1
      widget_control, state.dy_text_id, sensitive=1
      widget_control, state.dithsv_id, sensitive=1
      widget_control, state.nodith_id, sensitive=1
  endif else begin
      state.flgdith = 0
      widget_control, state.dx_text_id, sensitive=0
      widget_control, state.dy_text_id, sensitive=0
      widget_control, state.dithsv_id, sensitive=0
      widget_control, state.nodith_id, sensitive=0
  endelse


end



;;;;;;;;;
;  Update Display
;;;;;;;;;

pro xmkmask_UpdDisplay, state, MASK=mask

common xmkmask_images

  widget_control, /hourglass

; Mask
  if keyword_set( MASK ) then begin
      delvarx, mask_img
      mask_img = bytarr(img_size[0],img_size[1])
  endif

; Reset Display

  display_image = bytscl(main_image, min=state.pltmin, max=state.pltmax)

; Initial Mask

  if state.flgimsk EQ 1 then begin
      display_image[init_msk] = 0
      if keyword_set( MASK ) then mask_img[init_msk] = 1B
  endif

; Loop on Regions

  for q=1,state.nreg do begin
      case state.reg_typ[q] of
          0 : begin  ; Rectangle
              x0 = state.reg[q,0] > 0
              y0 = state.reg[q,1] > 0
              x1 = state.reg[q,2] < (img_size[0] - 1)
              y1 = state.reg[q,3] < (img_size[1] - 1)
              if(x0 LE x1 AND y0 LE y1) then begin
                  display_image[x0:x1,y0:y1] = 0
                  if keyword_set( MASK ) then  mask_img[x0:x1,y0:y1] = 1B
              endif
          end
          1 : begin ; Circle
              radius = sqrt( float(state.reg[q,0]-state.reg[q,2])^2 + $
                             float(state.reg[q,1]-state.reg[q,3])^2 )
              mskpix = xpix_circ(state.reg[q,0], $
                                 state.reg[q,1], radius, COUNT=count)
              if count NE 0 then begin
                  display_image[mskpix[0,*],mskpix[1,*]] = 0
                  if keyword_set( MASK ) then $
                    mask_img[mskpix[0,*], mskpix[1,*]] = 1B
              endif
          end
          2 : begin ; Line
              mskpix = xpix_line(state.reg[q,0],state.reg[q,1], $
                                state.reg[q,3],state.reg[q,4], $
                                state.reg[q,2], img_size[0]-1, img_size[1]-1, $
                               COUNT=count)
              if count NE 0 then begin
                  display_image[mskpix[0,*],mskpix[1,*]] = 0
                  if keyword_set( MASK ) then $
                    mask_img[mskpix[0,*], mskpix[1,*]] = 1B
              endif
          end
          else :
      endcase
  endfor

; Dither Image
;  if state.flgdith EQ 1 then begin
;      dith = where(dith_img EQ 1) 
;      display_image[dith] = 0
;      if keyword_set( MASK ) then mask_img[dith] = 1
;  endif

end

;;;;;;;;;
;  Dither
;;;;;;;;;

pro xmkmask_Dither, state, NODITH=nodith, SVDITH=svdith

common xmkmask_images


  widget_control, /hourglass

; Reset Display

  if not keyword_set( SVDITH ) then $
    display_image = bytscl(main_image, min=state.pltmin, max=state.pltmax)

; Initial Mask

  if state.flgimsk EQ 1 then begin
      display_image[init_msk] = 0
  endif

; NO DITH

  if keyword_set( NODITH ) then begin
      state.flgdith = 0
      state.nreg = 0
      return
  endif


; Loop on Regions

  if not keyword_set( SVDITH ) then begin
      for q=1,state.nreg do begin
          case state.reg_typ[q] of
              0 : begin         ; Rectangle
                  x0 = state.reg[q,0] + state.dithx > 0
                  y0 = state.reg[q,1] + state.dithy > 0
                  x1 = state.reg[q,2] + state.dithx < (img_size[0] - 1)
                  y1 = state.reg[q,3] + state.dithy < (img_size[1] - 1)
                                ; Check if it is even in the image
                  if(x0 LE x1 AND y0 LE y1) then begin
                      display_image[x0:x1,y0:y1] = 0
                  endif
              end
              1 : begin         ; Circle
                  radius = sqrt( float(state.reg[q,0]-state.reg[q,2])^2 + $
                                 float(state.reg[q,1]-state.reg[q,3])^2 )
                  mskpix = xpix_circ(state.reg[q,0]+state.dithx, $
                                     state.reg[q,1]+state.dithy, MAXX=img_size[0]-1, $
                                     MAXY=img_size[1]-1, radius, count=count)
                  if count NE 0 then begin
                      display_image[mskpix[0,*],mskpix[1,*]] = 0
                  endif
              end
              2 : begin         ; Line
                  mskpix = xpix_line(state.reg[q,0],state.reg[q,1], $
                                    state.reg[q,3],state.reg[q,4], $
                                    state.reg[q,2], img_size[0], img_size[1], $
                                    count=count)
                  for i=0L,count-1 do begin
                      x0 = mskpix[0,i] + state.dithx
                      y0 = mskpix[1,i] + state.dithy
                      if( (x0 GE 0) AND (x0 LE img_size[0]-1) AND $
                          (y0 GE 0) AND (y0 LE img_size[1]-1) ) then begin
                          display_image[x0,y0] = 0
                      endif
                  endfor
              end
              else :
          endcase
      endfor
  endif else begin    ; SVDITH
      for q=1,state.nreg do begin
          case state.reg_typ[q] of
              0 : begin         ; Rectangle
                  state.reg[q,0] = state.reg[q,0] + state.dithx
                  state.reg[q,1] = state.reg[q,1] + state.dithy
                  state.reg[q,2] = state.reg[q,2] + state.dithx 
                  state.reg[q,3] = state.reg[q,3] + state.dithy 
              end
              1 : begin         ; Circle
                  state.reg[q,0] = state.reg[q,0] + state.dithx
                  state.reg[q,1] = state.reg[q,1] + state.dithy
                  state.reg[q,2] = state.reg[q,2] + state.dithx 
                  state.reg[q,3] = state.reg[q,3] + state.dithy 
              end
              2 : begin         ; Line
                  state.reg[q,0] = state.reg[q,0] + state.dithx
                  state.reg[q,1] = state.reg[q,1] + state.dithy
                  state.reg[q,3] = state.reg[q,3] + state.dithx 
                  state.reg[q,4] = state.reg[q,4] + state.dithy 
              end
              else :
          endcase
      endfor
;
      state.flgdith = 1
      state.dithx = 0 
      widget_control, state.dx_text_id, set_value=0
      state.dithy = 0 
      widget_control, state.dy_text_id, set_value=0
  endelse

end
;;;;;;;;;
;  Set Regions
;;;;;;;;;

pro xmkmask_SetReg, state, flgb, REPLOT=replot

common xmkmask_images

  widget_control, /hourglass
  replot = 0

; flgb = 0 -> Button press: flgb=1 -> Release

  case state.regtyp of 
      0 : begin    ; Rectangle
          if flgb EQ 0 then begin   
              state.nreg = state.nreg + 1
              state.reg_typ[state.nreg] = 0
              state.reg[state.nreg,0] = nint(xgetx_plt(state.xcurs,state.pos, $
                                                 state.xymnx, state.winsize))
              state.reg[state.nreg,1] = nint(xgety_plt(state.ycurs,state.pos, $
                                                 state.xymnx, state.winsize))
          endif else begin
              xnew = nint(xgetx_plt(state.xcurs,state.pos, $
                               state.xymnx, state.winsize))
              ynew = nint(xgety_plt(state.ycurs,state.pos, $
                               state.xymnx, state.winsize))
              ; Require x1 < x2 AND y1 < y2
              if xnew LT state.reg[state.nreg,0] then begin
                  state.reg[state.nreg,2] = state.reg[state.nreg,0]
                  state.reg[state.nreg,0] = xnew
              endif else state.reg[state.nreg,2] = xnew
              if ynew LT state.reg[state.nreg,1] then begin
                  state.reg[state.nreg,3] = state.reg[state.nreg,1]
                  state.reg[state.nreg,1] = ynew
              endif else state.reg[state.nreg,3] = ynew
              ; Set some limits
              state.reg[state.nreg,0] = (0 > state.reg[state.nreg,0])
              state.reg[state.nreg,1] = (0 > state.reg[state.nreg,1])
              state.reg[state.nreg,2] = (img_size[0]-1 < $
                                         state.reg[state.nreg,2])
              state.reg[state.nreg,3] = (img_size[1]-1 < $
                                         state.reg[state.nreg,3])
              ; Set Image to 255
              display_image[state.reg[state.nreg,0]:state.reg[state.nreg,2], $
                            state.reg[state.nreg,1]:state.reg[state.nreg,3]] $
                = 0
              replot = 1
          end
      end
      1 : begin ; Circle
          if flgb EQ 0 then begin   
              state.nreg = state.nreg + 1
              state.reg_typ[state.nreg] = 1
              state.reg[state.nreg,0] = nint(xgetx_plt(state.xcurs,state.pos, $
                                                   state.xymnx, state.winsize))
              state.reg[state.nreg,1] = nint(xgety_plt(state.ycurs,state.pos, $
                                                   state.xymnx, state.winsize))
          endif else begin
              xnew = nint(xgetx_plt(state.xcurs,state.pos, $
                               state.xymnx, state.winsize))
              ynew = nint(xgety_plt(state.ycurs,state.pos, $
                               state.xymnx, state.winsize))
              state.reg[state.nreg,2] = xnew  
              state.reg[state.nreg,3] = ynew  
              ; Find Radius
              radius = sqrt( float(state.reg[state.nreg,0]-xnew)^2 + $
                             float(state.reg[state.nreg,1]-ynew)^2 )
              if radius LT 1. then begin    ; Error catch
                  display_image[xnew,ynew] = 0
                  replot = 1
                  return
              endif
              ; Find pixels
              mskpix = xpix_circ(state.reg[state.nreg,0], $
                                 state.reg[state.nreg,1], MAXX=img_size[0]-1, $
                                 MAXY=img_size[1]-1, radius, COUNT=count)
              if count NE 0 then display_image[mskpix[0,*],mskpix[1,*]] = 0
              replot = 1
          endelse
      end
      2 : begin ; Line
          case state.flgline of 
              0 : begin         ; endpoint
                  state.nreg = state.nreg + 1
                  state.reg_typ[state.nreg] = 2
                  state.reg[state.nreg,0] = nint(xgetx_plt(state.xcurs,state.pos, $
                                                           state.xymnx, state.winsize))
                  state.reg[state.nreg,1] = nint(xgety_plt(state.ycurs,state.pos, $
                                                           state.xymnx, state.winsize))
                  state.flgline = 1
              end
              1 : begin         ; width
                  x = xgetx_plt(state.xcurs,state.pos, $
                                     state.xymnx, state.winsize)
                  y = xgety_plt(state.ycurs,state.pos, $
                                state.xymnx, state.winsize)
                  width = sqrt( (state.reg[state.nreg,0] - x)^2 + $
                                (state.reg[state.nreg,1] - y)^2 )
                  state.reg[state.nreg,2] = round(width)
                  state.flgline = 2
              end
              2 : begin
                  state.reg[state.nreg,3] = nint(xgetx_plt(state.xcurs,state.pos, $
                                                           state.xymnx, state.winsize))
                  state.reg[state.nreg,4] = nint(xgety_plt(state.ycurs,state.pos, $
                                                           state.xymnx, state.winsize))
                  mskpix = xpix_line(state.reg[state.nreg,0],state.reg[state.nreg,1], $
                                    state.reg[state.nreg,3],state.reg[state.nreg,4], $
                                    state.reg[state.nreg,2], $
                                    img_size[0], img_size[1])
                  display_image[mskpix[0,*],mskpix[1,*]] = 0
                  replot = 1
                  state.flgline = 0
              end
          end
      end
      else :
  endcase
end
  
;;;;
; Delete Region(s)
;;;;

pro xmkmask_DelReg, state

common xmkmask_images

  widget_control, /hourglass
  if(state.nreg EQ 0) then return

; Find x,y and allow for dithering mode

  if state.flgdith NE -1 then begin
      xpix = nint(xgetx_plt(state.xcurs,state.pos, $
                            state.xymnx, state.winsize))
      ypix = nint(xgety_plt(state.ycurs,state.pos, $
                            state.xymnx, state.winsize))
  endif else begin
      xpix = nint(xgetx_plt(state.xcurs,state.pos, $
                            state.xymnx, state.winsize)) $
        - state.dithx
      ypix = nint(xgety_plt(state.ycurs,state.pos, $
                            state.xymnx, state.winsize)) $
        - state.dithy
  endelse

  flgb = 0

; Loop

  for q=1,state.nreg do begin
      case state.reg_typ[q] of
          0 : if( (xpix GE state.reg[q,0]) AND (xpix LE state.reg[q,2]) AND $
                  (ypix GE state.reg[q,1]) AND (ypix LE state.reg[q,3]) ) then flgb=1
          1 : begin ; Circle
              radius = sqrt( float(state.reg[q,0]-state.reg[q,2])^2 + $
                             float(state.reg[q,1]-state.reg[q,3])^2 )
              mskpix = xpix_circ(state.reg[q,0], $
                                 state.reg[q,1], MAXX=img_size[0]-1, $
                                 MAXY=img_size[1]-1, radius, COUNT=count)
              if count NE 0 then begin
                  a = where(mskpix[0,*] EQ xpix AND mskpix[1,*] EQ ypix, countB)
                  if( countB NE 0) then flgb=1
              endif
          end
          2 : begin ; Line
              mskpix = xpix_line(state.reg[q,0],state.reg[q,1], $
                                state.reg[q,3],state.reg[q,4], $
                                state.reg[q,2], img_size[0], img_size[1], $
                               COUNT=countA)
              if countA NE 0 then begin
                  a = where(mskpix[0,*] EQ xpix AND mskpix[1,*] EQ ypix, count)
                  if( count NE 0) then flgb=1
              endif
          end
          else :
      endcase
      if( q EQ state.nreg AND flgb NE 1) then begin
          print, 'You missed the region'
          return
      endif
      if( flgb EQ 1) then break
  endfor


; Adjust state.reg

  if (state.nreg NE 1) then begin
      for j=q,state.nreg-1 do begin
          state.reg[j,*] = state.reg[j+1,*]
          state.reg_typ[j] = state.reg_typ[j+1]
      endfor
  endif
  state.nreg = state.nreg - 1
  xmkmask_UpdDisplay, state

          
end

;;;;
; Output the Mask!
;;;;

pro xmkmask_Output, state

common xmkmask_images
  
  widget_control, /hourglass

; Name
  lastslsh = strpos(state.img[state.curimg],'/',/reverse_search)
  if lastslsh NE -1 then nopth = strmid(state.img[state.curimg], lastslsh+1) $
  else nopth = state.img[state.curimg]

; Take off the ov_ or f_
  if strmid(nopth,0,3) EQ 'ov_' then nopth = strmid(nopth,3)
  if strmid(nopth,0,2) EQ 'f_' then nopth = strmid(nopth,2)

; Add the Directory and header info
  if state.flgskmsk EQ 1 then begin
      outim = strjoin( [state.outdir, 'sm_', nopth] ) 
      bpm = strjoin( [state.outdir, 'sm_', nopth, '.gz'] ) 
      sxaddpar, header, 'SKYMASK', bpm
  endif else begin
      outim = strjoin( [state.outdir, 'mk_', nopth] ) 
      bpm = strjoin( [state.outdir, 'mk_', nopth, '.gz'] ) 
      sxaddpar, header, 'MASK', bpm
      sxaddpar, header, 'BPM', bpm
  endelse


; Write
  mwrfits, mask_img, outim, header, /silent, /create
  mwrfits, main_image, state.img[state.curimg], header, /silent, /create

; Compress
  spawn, 'gzip '+outim
 
return

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xmkmask_initcolors, state

common xmkmask_color

; Load a simple color table with the basic 8 colors in the lowest 
; 8 entries of the color table.  Also set top color to white.

rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
tvlct, 255*rtiny, 255*gtiny, 255*btiny

tvlct, [255],[255],[255], !d.table_size-1

end

;--------------------------------------------------------------------

pro xmkmask_getct, state, tablenum

common xmkmask_color

; Read in a pre-defined color table, and invert if necessary.

loadct, tablenum, /silent,  bottom=8
tvlct, r, g, b, 8, /get

xmkmask_initcolors, state

r = r[0:state.ncolors-2]
g = g[0:state.ncolors-2]
b = b[0:state.ncolors-2]

;if (state.invert_colormap EQ 1) then begin
;r = reverse(r)
;g = reverse(g)
;b = reverse(b)
;endif

r_vector = r
g_vector = g
b_vector = b

xmkmask_stretchct, state
;if (state.bitdepth EQ 24 AND (n_elements(pan_image) GT 10) ) then $
;  atv_refresh

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xmkmask_stretchct, state

; routine to change color stretch for given values of 
; brightness and contrast.
; Complete rewrite 2000-Sep-21 - Doug Finkbeiner
; This routine is now shorter and easier to understand.  

; if GETMOUSE then assume mouse positoin passed; otherwise ignore
; inputs

common xmkmask_color

x = state.brightness*(state.ncolors-1)
y = state.contrast*(state.ncolors-1) > 2   ; Minor change by AJB 
high = x+y & low = x-y
diff = (high-low) > 1

slope = float(state.ncolors-1)/diff ;Scale to range of 0 : nc-1
intercept = -slope*low
p = long(findgen(state.ncolors)*slope+intercept) ;subscripts to select
tvlct, r_vector[p], g_vector[p], b_vector[p], 8

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

pro xmkmask, img, IMSK=imsk, OUTDIR=outdir, XSIZE=xsize, YSIZE=ysize, $
             SKYMSK=skymsk, REGSTR=regstr

common xmkmask_color
common xmkmask_images

resolve_routine, 'xpix_circ', /either, /no_recompile
resolve_routine, 'xpix_line', /either

;
if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'xmkmask, img, IMSK=, OUTDIR=, XSIZE=,' + $
      'YSIZE=, /SKYMSK'
    return
  endif 

;  Optional Keywords

  if not keyword_set( OUTDIR ) then begin
      if not keyword_set ( SKYMSK ) then begin 
          a = findfile('Masks/..', count=count)
          if count EQ 0 then file_mkdir, 'Masks'
          outdir = 'Masks/' 
      endif else begin
          a = findfile('Masks/..', count=count)
          if count EQ 0 then file_mkdir, 'Masks'
          a = findfile('Masks/Sky/..', count=count)
          if count EQ 0 then file_mkdir, 'Masks/Sky'
          outdir = 'Masks/Sky/' 
      endelse
  endif

  if not keyword_set( IMSK ) then    imsk = ''

;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; INTERACTIVE

  if not keyword_set( XSIZE ) then    xsize = 800
  if not keyword_set( YSIZE ) then    ysize = 800

;    STATE
  state = {             $
            img: img, $
            nimg: 0, $
            curimg: 0, $
            imgmin: 0.0, $
            imgmax: 0.0, $
            imsk: imsk, $
            flgskmsk: 0, $
            flgimsk: 0, $
            outdir: outdir, $
            nreg: 0, $
            regtyp: 0, $ 
            reg_typ: intarr(1000), $ ; Type: 0 = Rect, 1= Circ
            reg: intarr(1000,10), $ 
            flgline: 0, $
            press: 0, $
            pos: [0.00,0.00,1.00,1.00], $ ; Plotting
            svxymnx: fltarr(4), $
            xymnx: fltarr(4), $
            winsize: intarr(2), $
            gridsize: intarr(2), $
            ncolors: 0, $
            brightness: 0.0, $
            contrast: 0.0, $
            pltmax: 0.0, $
            pltmin: 0.0, $
            xcurs: 0.0, $
            ycurs: 0.0, $
            dithx: 0, $         ; Dithering
            dithy: 0, $
            flgdith: 0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            drag_contr_id: 0L, $
            drag_brght_id: 0L, $
            draw_id: 0L, $
            dx_text_id: 0L, $
            dy_text_id: 0L, $
            reglist_id: 0L, $
            name_id: 0L, $
            dithsv_id: 0L, $
            nodith_id: 0L, $
            error_msg_id: 0L, $
            help_text_id: 0L $
          }
;  Images

  state.nimg = n_elements(img)
  state.curimg = 0

; Sky Mask

  if keyword_set( SKYMSK ) then state.flgskmsk = 1

;    WIDGET
  base = WIDGET_BASE( title = 'xmkmask: Make Mask', /column)
  state.base_id = base
  
;      Toolbar
  
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
;        Version + Name
  labelbase = widget_base(toolbar, /column, /base_align_center, /frame, $
                          /align_center)
  namelabel = WIDGET_LABEL(labelbase, value='xmkmask', /align_center)
  verslabel = WIDGET_LABEL(labelbase, value='(v1.1)', /align_center)
  state.name_id = WIDGET_LABEL(labelbase, value=strmid(img[0],0,18), /align_center)

;        Contrast
  contrbase = widget_base(toolbar, /column, /align_center)
  contrtext = widget_label(contrbase, value='Contr/Bright', /align_center)
  state.drag_contr_id = CW_FSLIDER(contrbase, TITLE='Contrast', $
                                   VALUE=0.5, uvalue='CONTRAST', $
                                   maximum=1.0, minimum=0.0, ysize=15, $
                                   /frame, /suppress_value)
;        Brightness
  state.drag_brght_id = CW_FSLIDER(contrbase, VALUE=0.5, uvalue='BRIGHT', $
                                   maximum=1.0, minimum=0.0, ysize=15, $
                                   TITLE='Brightness', /frame, /suppress_value)
;        Region Type
  regbase = widget_base(toolbar, /column, /align_center)
  reglabel = widget_label(regbase, value='Reg Type', /align_center)
  strreg = ['Rectangle','Circle','Line']
  state.reglist_id = WIDGET_LIST(regbase, VALUE=strreg, uvalue='REGLIST', $
                                 ysize = 3)
  widget_control, state.reglist_id, set_list_select=0
  
;        Dither
  dithbase = widget_base(toolbar, /column, /align_center)
  dxbase = widget_base(dithbase, /row, /align_right)
  state.dx_text_id = cw_field(dxbase, title='dx', /integer, $
                              value=0, xsize=5, ysize=1, /return_events, $
                                uvalue='DXTEXT')
  dybase = widget_base(dithbase, /row, /align_right)
  state.dy_text_id = cw_field(dybase, title='dy', /integer, $
                              value=0, xsize=5, ysize=1, /return_events, $
                                uvalue='DYTEXT')
  dbutbase = widget_base(dithbase, /row, /align_right)
  state.dithsv_id = WIDGET_BUTTON(dbutbase, value='SVDITH',uvalue='SVDITH')
  state.nodith_id = WIDGET_BUTTON(dbutbase, value='NODITH',uvalue='NODITH')
;        Help
  strhelp = strarr(50)
  strhelp = ['   Help Menu   ',$
             'Set Region to mask',$
             '  Rect: LMB+drag -> corner,corner',$
             '  Circ: LMB+drag -> center,edge',$
             '  Line: LMB+LMB+LMB -> center,edge,other center',$
             'Delete Region -- RMB', $
             'Dither (need to have created a mask first)', $
             '  dx,dy: Set value+return', $
             '  Save: SVDITH',$
             '  Erase dither mask: NODITH',$
             'Skip Mask -- NOSAVE', $
             'Save Mask -- NEXT' $
            ]
  help_text_id = widget_text(toolbar, value=strhelp, xsize=30, ysize=6,$
                             /scroll)
;      Error base
;  error_base_id  = widget_base(toolbar, /column)
;  error_btn_id = widget_button(error_base_id, value='Erase Error Msg', $
;                               uvalue='ERRORB')
;          Error Messages
;  state.error_msg_id = widget_text(error_base_id, value='', xsize=40, $
;                                   ysize=2, /scroll)
  
;      Drawing
  state.winsize[0] = xsize
  state.winsize[1] = ysize
  state.draw_id = widget_draw(base, xsize=state.winsize[0], $
                              ysize=state.winsize[1], /frame, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Done
  donebase = widget_base(toolbar, /column, /align_right)
  done = WIDGET_BUTTON(donebase, value='NEXT',uvalue='NEXT')
  nosave = WIDGET_BUTTON(donebase, value='No Save',uvalue='NOSAVE')
  bigq = WIDGET_BUTTON(donebase, value='DONE!',uvalue='DONE')
  
; Realize
  WIDGET_CONTROL, base, /realize

; Initial Mask

  if keyword_set( IMSK ) then state.flgimsk = 1
  
; Colors 
  loadct, 0, /silent
  if (!d.table_size LT 12) then begin
      message, 'Too few colors available for color table'
      stop
  endif
  state.ncolors = !d.table_size - 9
  r_vector = bytarr(state.ncolors)
  g_vector = bytarr(state.ncolors)
  b_vector = bytarr(state.ncolors)

  xmkmask_getct, state, 0

  xmkmask_Reset, state
; Update
  xmkmask_ReadImg, state
;  xmkmask_UpdateImg, base
  xmkmask_PlotImg, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'xmkmask', base

; Output
;  ffit = *pnt_ffit
;  PTR_FREE, pnt_ffit
  REGSTR = { $
           nreg: cmmn_nreg, $
           regions: cmmn_regions, $
           reg_type: cmmn_regtype $
           }
            
  
  return
end

