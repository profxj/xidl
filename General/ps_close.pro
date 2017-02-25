;+
; Name:
;
;  PS_Close
;
; Purpose:
;
;    Closes the PostScript plotting device (opened with PS_Open),
;    adds an optional time-stamp in the lower left corner,  
;    sends the PostScript file to the default printer ('lw'),
;    and restore the plotting device saved by PS_Open.
;
; Calling sequence:
;
;    PS_Close
;
; Keywords:
;
;    NOID - If present, the time-stamp ID is not added to the plot.
;  
;    NOPRINT - If present, the PostScript plotting device is closed, and 
;      the plotting device saved by PS_Open is restore, 
;      but the PostScript file is NOT sent to the printer. 
;  
;    SAVEAS = <string>. Save the PostScript file by copying it to the
;      specified file. 
;  
;    PRINTER = <string>. Send the plot the specified printer.
;      This is keyword ignored for color plots.
;
;    COMMENT - A string to include with time-stamp
;
;    TRANSPARENCY - If /COLOR has been set in PS_Open queue the plot to be
;      printed on 8.5x11 transparency, unless /LEDGER has been specified in
;      PS_Open.
;
;    LLWR - for very large plots, use the flag /LLWR, to use my 'llwr' script
;
; See also:
;    PS_Open
;
; History:
;   17-Sep-93: added color flag to common block 
;              hence uses /h/sylvain/sbin/lwiplr to print CPS  
;    9-Apr-94: added stuff to info messages
;   13-Apr-93: cleaned up "transparency" bug
;              added LEDGER support
;   11-May-94: added /LLWR
;
PRO Ps_close, noprint = np, saveas = fn, printer = pnk, noid = noid, $
             transparency = trsp, llwr = use_llwr, comment = comm, silent=silent
;
  COMMON ps_common, old_dname, old_pfont, file_name, opened, color, ledger
;
  lwiplr = '/h/sylvain/sbin/lwiplr'
  llwr = '/usr4/sabbey/scripts/llwr'
;
  IF n_elements(opened)  EQ 0 THEN opened = 0
  IF opened EQ 0 THEN  BEGIN
    message, /info, $
      'ERROR: PS device not opened, you must open it first (using PS_Open)'
    return
  ENDIF
;
; add a time stamp first
  IF NOT keyword_set(comm) then comm=''
  IF NOT keyword_set(noid) THEN put_id,size=.75,comment=comm
  device, /close
  opened = 0
;
  IF NOT keyword_set(pnk) THEN pn = 'lw' ELSE pn = pnk
;
  IF keyword_set(fn) THEN BEGIN
    cmd = ['cp', file_name, fn]
    spawn, cmd, /noshell
    if ~keyword_set(silent) then print, 'PostScript plot saved as ', fn
  END
;
  IF NOT keyword_set(np) THEN BEGIN
    IF color THEN BEGIN
      IF keyword_set(pnk) THEN $
        if ~keyword_set(silent) then message, /info, 'PRINTER specification ('+pn+') ignored'
      IF ledger EQ 1 THEN BEGIN
        mode =  'ledger'
      ENDIF ELSE BEGIN
        IF keyword_set(trsp) THEN BEGIN
          mode = 'transparency'
        ENDIF ELSE BEGIN
          mode = 'letter'
        ENDELSE
      ENDELSE
      
;      message, /info, $
;        'Color PostScript plot printed on WIPL printer ('+mode+')'
;      message, /info, $
;        '(using '+lwiplr+' -'+mode+')'
;      spawn, [lwiplr, '-'+mode, file_name], /noshell 
    stop,'Ahem, we do not have a color printer, save it to a file'
      
    ENDIF ELSE BEGIN
      IF keyword_set(use_llwr) THEN BEGIN
;        spawn, [llwr, file_name], /noshell
;        message, /info, 'PostScript plot printed using llwr'
         stop,'Ahem, do not use this option, or ask Chris to fix it.'
      ENDIF ELSE BEGIN
        spawn, ['lpr', '-P'+pn, file_name], /noshell
        if ~keyword_set(silent) then message, /info, 'PostScript plot printed on printer '+pn
      ENDELSE
    ENDELSE
  END
;  
  set_plot, old_dname
  !p.font =  old_pfont
  if ~keyword_set(silent) then message, /info, 'plotting device restored as ' + old_dname
;  
  return
END
