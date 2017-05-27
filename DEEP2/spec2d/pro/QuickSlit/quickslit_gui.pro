;+
;
; NAME
;     quickslit_gui.pro
;
; PURPOSE
;     quickslit_gui manages the GUI for the DEIMOS QuickSlit reduction
;     package. the routine launches an IDL widget which allows the
;     user to specify the DEIMOS image file(s) and the DEIMOS slit
;     number which are then passed to quickslit.pro. 
;
; SYNTAX
;     quickslit_gui
;
; INPUTS
;     None.
;
; KEYWORDS
;     None.
;
; OUTPUTS
;     The routine launches the QuickSlit GUI. 
;
; PROCEDURES CALLED 
;     quickslit.pro
;     
; EXAMPLES
;     None.
;
; COMMENTS
;     quickslit_gui also handles the error message reporting for the
;     QuickSlit reduction package. Errors caught by quickslit pass
;     messages back to quicklsit_gui and are displayed via the GUI.
;
; HISTORY
;     Created May 19, 2004 by mcc.
;
;-


pro quickslit_getfile, event
  common quickslit_cblock, file, slitno
  widget_control, event.id, get_value=file
  widget_control, event.top, get_uvalue=wstate
  wstate.file = file
end

pro quickslit_getslitno, event
  common quickslit_cblock, file, slitno
  widget_control, event.id, get_value=slitno
  widget_control, event.top, get_uvalue=wstate
  wstate.slitno = slitno
end

pro quickslit_run, event
  common quickslit_cblock, file, slitno
; get the start sign.
  widget_control, event.top, get_uvalue=wstate

  print
  print, '-------------------------------------------'
  print, 'Reducing file(s): ' + file
  print, 'Slit Number: ' + string(slitno,  format='(I3.3)')
  print, '-------------------------------------------'
  print

  nbase = widget_base(/column,  title='STATUS',  $
                     group_leader=event.top,  /align_center,  $
                      scr_ysize=200,  scr_xsize=200)
  label =  widget_label(nbase,  value='Reduction Running....')
  widget_control,  nbase,  /realize

  quickslit,  file[0], slitno[0], mess=mess, errval=errval
  if errval eq 1 then begin
     print, '-------------------------------'
     print, 'ERROR'
     print, '-------------------------------'
     ebase0 =  widget_base(/col,  title='ERROR', /align_center,  $
                           group_leader=event.top, xsize=1000)
     ebase1 = widget_base(ebase0, /row, /align_center)
     label =  widget_label(ebase0, value=mess,  ysize=500)
     ebase2 =  widget_base(ebase0, /row,  /align_center)
     ok = widget_button(ebase0, value='OK', xsize=50, $
                        event_pro='quickslit_redo', /align_center)

     widget_control, ebase0, /realize
  endif

  if errval eq 0 then begin
     ebase0 =  widget_base(/col,  title='SUCCESS', /align_center,  $
                           group_leader=event.top, xsize=1000)

     mess1 = 'This reduction has been brought to you by'
     mess2 = 'the members of the UC-Berkeley DEEP2 team.'
     mess3 = 'Any monetary donations should be sent to '
     mess4 = '     Michael Cooper   '
     mess5 = '     601 Campbell Hall   '
     mess6 = '     University of California  '
     mess7 = '     Berkeley, CA 94720-3411   '
     ebase1 = widget_base(ebase0, /row, /align_center)     
     label1 =  widget_label(ebase0, value=mess1)
     ebase2 = widget_base(ebase0, /row, /align_center)     
     label2 =  widget_label(ebase0, value=mess2)
     ebase3 = widget_base(ebase0, /row, /align_center)     
     label3 =  widget_label(ebase0, value=mess3)
     ebase4 = widget_base(ebase0, /row, /align_center)     
     label4 =  widget_label(ebase0, value=mess4)
     ebase5 = widget_base(ebase0, /row, /align_center)     
     label5 =  widget_label(ebase0, value=mess5)
     ebase6 = widget_base(ebase0, /row, /align_center)     
     label6 =  widget_label(ebase0, value=mess6)
     ebase7 = widget_base(ebase0, /row, /align_center)     
     label7 =  widget_label(ebase0, value=mess7)

     

     ebaseB =  widget_base(ebase0, /row,  /align_center)
     ok = widget_button(ebase0, value="YOU'RE WELCOME", xsize=200, $
                        event_pro='quickslit_redo', /align_center)

     widget_control, ebase0, /realize
  endif


  widget_control, nbase, /destroy

end


pro quickslit_redo, event
  widget_control, event.top, /destroy
end


pro quickslit_exit, event
; get the start sign.
  widget_control, event.top, /destroy
  retall
end



pro quickslit_gui
; define our common block.
  common quickslit_cblock, file, slitno
; set variables to default values.
  file = 'NOFILE'
  slitno = 0L

; make the widget bases.
  base = widget_base(/col, title='QUICKSLIT')
  base0a = widget_base(base, /row,  /align_center)
  base0b = widget_base(base, /row,  /align_center)
  base0c = widget_base(base, /row,  /align_center)
  base1 = widget_base(base, /row)
  base2 = widget_base(base, /row)
  base3 = widget_base(base, /row)
  base4 = widget_base(base, /row)

; make the author list.
  label = widget_label(base0a, value='QUICKSLIT Reduction Package (v0.1)')
  label = widget_label(base0b, value='Created by Michael Cooper')
  label = widget_label(base0c, value='UC-Berkeley')

; make the frame entry.
  label = widget_label(base1, value='Enter Science Image(s):')
  text1  = widget_text(base1, xsize=50, uvalue=0L, $
                       event_PRO='quickslit_getfile', $
                       value='', /editable, /all_events)

; make the slitno entry.
  label = widget_label(base2, value='Enter Slit Number:')
  text2  = widget_text(base2, xsize=50, uvalue=0L, $
                       event_pro='quickslit_getslitno', $
                       value='', /editable, /all_events)


; add GO button.
  start = widget_button(base3, value='START', xsize=500, $
                        event_pro='quickslit_run', /align_center)

; add EXIT button.
  start = widget_button(base4, value='EXIT', xsize=500, $
                        event_pro='quickslit_exit', /align_center)




; define the state structure.
  wstate = {file:file, slitno:slitno, start:start}

  widget_control, base, /realize
  widget_control, base, set_uvalue=wstate

  xmanager, 'quickslit_start', base


end
