;+ 
; NAME:
; windowtofile
;
;   WINDOW_TO_FILE.PRO   5-14-97  DSFoster   [Unix]
;
; Routine to save contents of specified window to graphics file.
; Specify format using keyword: GIF, PICT, TIFF, BMP, JPEG. 
; Defaults to TIFF.
;
; Also use WINDOW_TO_FILE() to print graphics file containing window
; contents, using /PRINT keyword. If the argument FILENAME is the null
; string '' then a temporary file will be created for printing, and
; will be deleted after printing.
;
; Creates new pixmap window and reads from that, since reading
; from scrollable draws using TVRD() gives incorrect results (IDL bug!).
;
; Return values are:
;    0 : No error
;       -1 : Error writing to file
;       -2 : Specified window is unavailable or not open
;       -3 : Unknown graphics file format specified
;
; Modifications:
;
;  2-19-96 DSF Bug fix: give appropriate X and Y dimensions in 
;          DEVICE, COPY=[] command.
;  7-19-96 DSF Initialize argument "message".
;  5-14-97 DSF Add JPEG format.
 

FUNCTION window_to_file, window, filename, message, FORMAT=format, PRINT=print

message = ''

if (keyword_set(FORMAT)) then begin
        format = strupcase(format)
endif else begin
        format = 'TIFF'
endelse

; Determine "quality" of JPEG image

if (keyword_set(QUALITY)) then begin
        if (quality lt 0 or quality gt 100) then begin
                message, 'Keyword QUALITY must be between 0 and 100', /continue
                return, -3
        endif else if (format ne 'JPEG' and format ne 'JPG') then begin
                message, 'Keyword QUALITY requires FORMAT=JPEG', /continue
                return, -3
        endif
        quality = fix( round(quality) )
endif else begin
        quality = 75                                ; Default to 75 (Very good)
endelse

device, window_state=windows
if (windows(window) le 0) then begin
        message, 'Window unavailable', /continue
        return, -2
endif

if (filename eq '') then begin                  ; Generate temporary name
        t = string(bin_date(systime()))                  ; Use date/time
        t = t(0) + t(1) + t(2) + t(3) + t(4) + t(5)
        fname = '/tmp/IDL_' + strcompress(t, /remove_all)
endif else begin
        fname = filename
        openw, unit, fname, /get_lun, error=err     ; Else is file writable?
        if (err ne 0) then begin
                message, ' Error writing to file: ' + fname, /continue
                return, -1
        endif
        free_lun, unit
endelse

old_window = !d.window
wset, window                                    ; Get size of specified window
xsize = !d.x_size
ysize = !d.y_size

window, xsize=xsize, ysize=ysize, /free, /pixmap
DEVICE, copy=[0,0, xsize-1,ysize-1, 0,0, window]
array = tvrd()                                  ; Window contents to array

; Reverse array if necessary
if (!order ne 0 and format ne 'TIFF' and format ne 'TIF') then $
        array = reverse(temporary(array),2)

clr = getcolor(/load)
TVLCT, r,g,b, /GET                              ; Get color table

case (format) of
;    'TIFF': tiff_write, fname, array, red=clr.red, green=clr.green, $
;      blue=clr.blue
    'TIFF': tiff_write, fname, array, red=r, green=g, blue=b
    'TIF': tiff_write, fname, array, red=r, green=g, blue=b
    'BMP': write_bmp, fname, array, r,g,b
    'PICT': write_pict, fname, array, r,g,b
    'PIC': write_pict, fname, array, r,g,b
    'GIF': write_gif, fname, array, r,g,b
    'JPEG': write_jpeg, fname, array, quality=quality, order=0
    'JPG': write_jpeg, fname, array, quality=quality, order=0
    else: begin
        message, 'Unknown file format specified with FORMAT: ' + format, $
          /continue
        if (old_window ne -1) then wset, old_window
        return, -3
    end
endcase

if (keyword_set(PRINT)) then begin
        if (filename eq '') then begin       ; Make copy of file for print (-c)
                command = 'lp -c -w ' + fname + ' 2>&1'
        endif else begin
                command = 'lp -w ' + fname + ' 2>&1'
        endelse

        SPAWN, ["/bin/sh", "-c", command], results, /NOSHELL   ; PRINT IT!

        if (n_params() ge 3) then message = results
        if (filename eq '') then begin       ; Delete temporary file
                ret = delete_files(fname)
                if (ret ne 0) then $
                        message, 'Error deleting file: ' + fname, /continue
        endif
endif
                                                                                                                          
if (old_window ne -1) then wset, old_window
return, 0
END

