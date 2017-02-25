; $Id: x_widget_displayfile.pro,v 1.1 2008-04-04 18:36:11 xavier Exp $
;
; Copyright (c) 1991-2006, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.

PRO x_widget_displayfile_write, wText, filename

  COMPILE_OPT hidden

  WIDGET_CONTROL, /HOURGLASS
  OPENW, unit, FILENAME, /GET_LUN, ERROR=i		;open the file and then
  if i lt 0 then begin		;OK?
	a = [ !error_state.msg, filename + ' could not be opened for writing.']  ;No
	void = DIALOG_MESSAGE(a, /ERROR, DIALOG_PARENT=wText)
  endif else begin
	WIDGET_CONTROL, wText, get_value=txtArray
	ON_IOERROR, done_writing
	; print out each line separately in order to get desired line breaks
	for j = 0, N_ELEMENTS(txtArray)-1 do PRINTF, unit, txtArray[j]
done_writing:
	ON_IOERROR, null
	FREE_LUN, unit				;free the file unit.
  endelse
END

PRO x_widget_displayfile_event, event

  COMPILE_OPT hidden

  WIDGET_CONTROL, event.top, GET_UVALUE=state
  CASE TAG_NAMES(event, /STRUCTURE_NAME) OF
    'WIDGET_BASE': BEGIN
      w = (event.x gt state.x_reserve) $
	? (event.x - state.x_reserve) : state.x_reserve
      if (!version.os_family eq 'Windows') then begin
	h=event.y
      endif else begin
        h = (event.y gt state.y_reserve) $
	  ? (event.y - state.y_reserve) : state.y_reserve
      endelse
      WIDGET_CONTROL, state.filetext, scr_xsize=w, scr_ysize=h

      RETURN
      END
    'WIDGET_KILL_REQUEST': retval = "EXIT"
    ELSE: WIDGET_CONTROL, event.id, GET_UVALUE = retval
  ENDCASE

  CASE retval OF
  	"SAVE": BEGIN
               if (LMGR(/DEMO)) then begin
                  tmp = DIALOG_MESSAGE( /ERROR, $
                        'Save: Feature disabled for demo mode.')
                  return
                endif
		IF (STRLEN(state.filename) EQ 0) THEN BEGIN
			state.filename = DIALOG_PICKFILE(/WRITE)
		ENDIF
		IF (STRLEN(state.filename) GT 0) THEN BEGIN
			x_widget_displayfile_write, state.filetext, state.filename
			WIDGET_CONTROL, event.top, SET_UVALUE=state
			IF state.notitle THEN WIDGET_CONTROL, event.top, $
				TLB_SET_TITLE=state.filename
		ENDIF

		RETURN
  	END
	"SAVE_AS": BEGIN
               if (LMGR(/DEMO)) then begin
                  tmp = DIALOG_MESSAGE( /ERROR, $
                        'Save As: Feature disabled for demo mode.')
                  return
                endif
		state.filename = DIALOG_PICKFILE(/WRITE)
		IF (STRLEN(state.filename) GT 0) THEN BEGIN
			x_widget_displayfile_write, state.filetext, state.filename
			WIDGET_CONTROL, event.top, SET_UVALUE=state
			IF state.notitle THEN WIDGET_CONTROL, event.top, $
				TLB_SET_TITLE=state.filename
		ENDIF
		RETURN
	END
	"EXIT": BEGIN
		WIDGET_CONTROL, event.top, /DESTROY
		IF (WIDGET_INFO(state.ourGroup, /VALID)) THEN $
			WIDGET_CONTROL, state.ourGroup, /DESTROY
	END
	ELSE:
  ENDCASE
END






PRO XDisplayFileGrowToScreen, tlb, text, height, nlines
; Grow the text widget so that it displays all of the text or
; it is as large as the screen can hold.

  max_y = (get_screen_size())[1] - 100
  cur_y = (WIDGET_INFO(tlb, /geometry)).scr_ysize

  ; If the display is already long enough, then there's nothing to do.
  ;
  ; We are only filling to grow the display, not shrink it, so if its
  ; already too big, there's nothing to do. This can only happen if
  ; the caller sets a large HEIGHT keyword, which is operator error.
  if ((nlines le height) || (cur_y gt max_y)) then return

  ; The strategy I use is binary divide and conquer. Furthermore,
  ; if nlines is more than 150, I limit the search to that much.
  ; (Consider that a typical screen is 1024 pixels high, and that
  ; using 10pt type, this yields 102 lines). This number may need to
  ; be adjusted as screen resolution grows but that will almost certainly
  ; be a slowly moving and easy to track target.
  ;
  ; Note: The variable cnt should never hit its limit. It is there
  ; as a "deadman switch".
  low = height
  high = MIN([150, nlines+1])
  cnt=0
  while ((low lt high) && (cnt++ lt 100)) do begin
    old_low = low
    old_high = high
    mid = low + ((high - low + 1) / 2)
    WIDGET_CONTROL, text, ysize=mid
    cur_y = (WIDGET_INFO(tlb, /geometry)).scr_ysize
    if (cur_y lt max_y) then low = mid else high = mid
    if ((old_low eq low) && (old_high eq high)) then break
  endwhile

end







PRO x_widget_displayfile, FILENAME, TITLE = TITLE, $
  GROUP = GROUP, WIDTH = WIDTH, $
  HEIGHT = HEIGHT, TEXT = TEXT, FONT = font, $
  DONE_BUTTON=done_button, MODAL=MODAL, $
  EDITABLE=editable, GROW_TO_SCREEN=grow_to_screen, $
  WTEXT=filetext, BLOCK=block, RETURN_ID=return_id, WRAP=wrap
;+
; NAME:
;	x_widget_displayfile
;
; PURPOSE:
;	Display an ASCII text file using widgets and the widget manager.
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	x_widget_displayfile, Filename
;
; INPUTS:
;     Filename:	A scalar string that contains the filename of the file
;		to display.  The filename can include a path to that file.
;
; KEYWORD PARAMETERS:
;	BLOCK:  Set this keyword to have XMANAGER block when this
;		application is registered.  By default the Xmanager
;               keyword NO_BLOCK is set to 1 to provide access to the
;               command line if active command 	line processing is available.
;               Note that setting BLOCK for this application will cause
;		all widget applications to block, not only this
;		application.  For more information see the NO_BLOCK keyword
;		to XMANAGER.
;
;	DONE_BUTTON: the text to use for the Done button.  If omitted,
;		the text "Done with <filename>" is used.
;
;	EDITABLE: Set this keyword to allow modifications to the text
;		displayed in x_widget_displayfile.  Setting this keyword also
;		adds a "Save" button in addition to the Done button.
;
;	FONT:   The name of the font to use.  If omitted use the default
;		font.
;	GROUP:	The widget ID of the group leader of the widget.  If this
;		keyword is specified, the death of the group leader results in
;		the death of x_widget_displayfile.
;
;       GROW_TO_SCREEN: If TRUE, the length of the display area is grown
;		to show as much of the text as possible without being too
;		large to fit on the screen. In this case, HEIGHT sets the
;		lower bound on the size instead of setting the size itself.
;
;	HEIGHT:	The number of text lines that the widget should display at one
;		time.  If this keyword is not specified, 24 lines is the
;		default.
;
;       RETURN_ID : A variable to be set to the widget ID of the top level
;               base of the resulting help application.
;	TEXT:	A string or string array to be displayed in the widget
;		instead of the contents of a file.  This keyword supercedes
;		the FILENAME input parameter.
;
;	TITLE:	A string to use as the widget title rather than the file name
;		or "XDisplayFile".
;
;	WIDTH:	The number of characters wide the widget should be.  If this
;		keyword is not specified, 80 characters is the default.
;
;	WTEXT:	Output parameter, the id of the text widget.  This allows
;		setting text selections and cursor positions programmatically.
;
; OUTPUTS:
;	No explicit outputs.  A file viewing widget is created.
;
; SIDE EFFECTS:
;	Triggers the XMANAGER if it is not already in use.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Open a file and create a widget to display its contents.
;
; MODIFICATION HISTORY:
;	Written By Steve Richards, December 1990
;	Graceful error recovery, DMS, Feb, 1992.
;       12 Jan. 1994  - KDB
;               If file was empty, program would crash. Fixed.
;       4 Oct. 1994     MLR Fixed bug if /TEXT was present and /TITLE was not.
;	2 jan 1997	DMS Added DONE_BUTTON keyword, made Done
;			button align on left, removed padding.
;	19 Nov 2004, GROW_TO_SCREEN and RETURN_ID keywords. Allow for
;                       user to resize display. General updating.
;-

  ; Establish defaults if keywords not specified
  IF(NOT(KEYWORD_SET(EDITABLE))) THEN editable = 0
  IF(NOT(KEYWORD_SET(HEIGHT))) THEN HEIGHT = 24
  IF(NOT(KEYWORD_SET(WIDTH))) THEN WIDTH = 80
  IF N_ELEMENTS(block) EQ 0 THEN block=0
  noTitle = N_ELEMENTS(title) EQ 0

  IF(NOT(KEYWORD_SET(TEXT))) THEN BEGIN
    IF noTitle THEN TITLE = FILENAME

    unit = -1
    CATCH, err
    if (err ne 0) then begin
      CATCH,/CANCEL
      if (unit ne -1) then FREE_LUN, 1
      a = [ !error_state.msg, ' Unable to display ' + filename]
      nlines = n_elements(a)
    endif else begin
      nlines = MIN([FILE_LINES(filename), 10000])
      OPENR, unit, FILENAME, /GET_LUN
      a = strarr(nlines)
      readf, unit, a
      CATCH, /CANCEL
      FREE_LUN, unit
    endelse
  ENDIF ELSE BEGIN
    IF(N_ELEMENTS(FILENAME) EQ 0) THEN FILENAME=''
    IF noTitle THEN TITLE = 'XDisplayFile'
    a = TEXT
    nlines = n_elements(a)
  ENDELSE

  ourGroup = 0L
  if KEYWORD_SET(MODAL) then begin
    if N_ELEMENTS(GROUP) GT 0 then begin
      filebase = WIDGET_BASE(TITLE = TITLE, $
			/TLB_KILL_REQUEST_EVENTS, TLB_FRAME_ATTR=1, $
                        /BASE_ALIGN_LEFT, /COLUMN, $
                        /MODAL, GROUP_LEADER = GROUP)
    endif else begin
      ; modal requires a group leader
      ourGroup = WIDGET_BASE()
      filebase = WIDGET_BASE(TITLE = TITLE, $
        		/TLB_KILL_REQUEST_EVENTS, TLB_FRAME_ATTR=1, $
                        /BASE_ALIGN_LEFT, /COLUMN, $
                        /MODAL, GROUP_LEADER = ourGroup)
    endelse
    menu_bar = filebase
  endif else begin
    filebase = WIDGET_BASE(TITLE = TITLE, $
        		/TLB_KILL_REQUEST_EVENTS, /TLB_SIZE_EVENTS, $
                        /BASE_ALIGN_LEFT, /COLUMN, MBAR=menu_bar, $
			GROUP_LEADER = GROUP)
  endelse
  return_id = filebase


  extra = ''
  IF (menu_bar NE filebase) THEN BEGIN
    IF (!VERSION.OS_FAMILY EQ 'Windows') THEN extra = '&'
    menu_bar = WIDGET_BUTTON(menu_bar, VALUE=extra+'File', /MENU)
  ENDIF ELSE $
    menu_bar = WIDGET_BASE(filebase, /ROW)


  IF (editable) THEN BEGIN
    ; add 'Save', 'Save as...' buttons here
    saveButton = WIDGET_BUTTON(menu_bar, VALUE = extra+'Save', UVALUE = "SAVE")
    saveAsButton = WIDGET_BUTTON(menu_bar, $
		VALUE = 'Save '+extra+'As...', UVALUE = "SAVE_AS")

  ENDIF

  ; Done button
  if n_elements(done_button) eq 0 then done_button = "Done with " + TITLE
  filequit = WIDGET_BUTTON(menu_bar, SEPARATOR=editable, $
		VALUE = extra+done_button, UVALUE = "EXIT")

  ; Create a text widget to display the text
  IF n_elements(font) gt 0 then begin
   filetext = WIDGET_TEXT(filebase, XSIZE = WIDTH, YSIZE = HEIGHT, $
		EDITABLE = editable, UVALUE='TEXT', /SCROLL, VALUE = a, $
		FONT = font, WRAP=wrap)
  endif else begin
    filetext = WIDGET_TEXT(filebase, XSIZE = WIDTH, YSIZE = HEIGHT, $
		EDITABLE = editable, UVALUE='TEXT', /SCROLL, VALUE = a,$
                          WRAP=wrap)
  endelse

if (keyword_set(grow_to_screen)) then $
  XDisplayFileGrowToScreen, filebase, filetext, height, nlines

WIDGET_CONTROL, filebase, /REALIZE

geo_base = WIDGET_INFO(filebase, /geometry)
geo_text = WIDGET_INFO(filetext, /geometry)

state={ ourGroup:ourGroup, filename: filename, $
        filetext:filetext, notitle:noTitle, $
        x_reserve:geo_base.scr_xsize - geo_text.scr_xsize, $
        y_reserve:geo_base.scr_ysize - geo_text.scr_ysize }
WIDGET_CONTROL, filebase, SET_UVALUE = state


xmanager, "x_widget_displayfile", filebase, GROUP_LEADER = GROUP, $
	NO_BLOCK=(NOT(FLOAT(block)))

end
