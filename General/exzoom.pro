
PRO widzoom_event, event
   WIDGET_CONTROL, event.id, GET_UVALUE=uvalue
   CASE uvalue OF
      'ZOOM': HELP, /STRUCT, event
      'DONE': WIDGET_CONTROL, event.top, /DESTROY
   ENDCASE
END

;Next, create the widget program:

PRO widzoom
   OPENR, lun, FILEPATH('people.dat', $
      SUBDIR = ['examples','data']), /GET_LUN
   image=BYTARR(192,192)
   READU, lun, image
   FREE_LUN, lun
   sz = SIZE(image)
   base=WIDGET_BASE(/COLUMN)
   zoom=CW_ZOOM(base, XSIZE=sz[1], YSIZE=sz[2], UVALUE='ZOOM')
   done=WIDGET_BUTTON(base, VALUE='Done', UVALUE='DONE')
   WIDGET_CONTROL, base, /REALIZE
   WIDGET_CONTROL, zoom, SET_VALUE=image
   XMANAGER, 'widzoom', base
END
