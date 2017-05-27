;+
; NAME:
;    slitchooser
;
; PURPOSE:
;    initiates a widget used to prepare for quicklook analysis
;
; CALLING SEQUENCE:
;     slitchooser
; 
; INPUTS:
;     arguments are specified in the widget
;
; OUTPUTS:
;     creates a series of xxxx.ql files for quicklook reduction, xxxx
;     is mask number.  These are formatted text files
;
; COMMENTS:
;     numerous masks can be created at one call of this widget.
;     Allows user to select slitlets for subsequent analysis
;
; REVISION HISTORY:
;    created by jcdt -Jul02
;    comments by MD, added polyflag 6sep02     
;----------------------------------------------------------------------


function isPosInt,  instring
  
                                ; This function returns 0 if the
                                ; string is not a postive integer, and
                                ; returns the integer otherwise

  len = (strlen(instring))[0]
  if (len eq 0) then return,0

  for i = 0L, len-1 DO  $
     if (strmid(instring, i, 1) gt '9' OR strmid(instring, i, 1) lt '0') $
           then return, 0

  num = fix(instring)
  
  if (num eq 0) then return, 0
  
  return, num
END



PRO resetFields, idList
                                ; This resets all the data entry
                                ; fields, except the data directory,
                                ; which will most likely be the same
                                ; for each file generated.
  
  widget_control, idList.mask, set_value=''
  ;widget_control, idList.dataDir, set_value=''
  widget_control, idList.flatname, set_value=''
  widget_control, idList.arcname, set_value=''
  widget_control, idList.slitNums, set_value=''

  for i=0, 3 do begin
     widget_control, idList.nArray[i, 0], set_value=0
     widget_control, idList.nArray[i, 1], set_value=''
  endfor

END


PRO setCalibFiles, maskID, dataDir, idList
                                ; This procedure find the calibration
                                ; files for the given mask, in the
                                ; given raw data directory.
  
  ; This code is not yet available.
  go = dialog_message ("This feature is not ready yet.")
  return

END

PRO slitChooser_event,  ev
 
  WIDGET_CONTROL, ev.top,  GET_UVALUE=idList
  WIDGET_CONTROL, ev.id, GET_UVALUE=uval
  
  

  CASE uval OF
     'GO': begin 
                                ; when the 'go' button is pressed, it
                                ; checks to make sure all the data is
                                ; in the correct format, then prints
                                ; it out (to a file?) to be parsed
                                ; later, then quits.


      ; Input the mask ID:
        
        Widget_control,  idList.mask,  get_value=tempStr  
        maskID= isPosInt (tempStr) 
        
        if (maskID le 0) then begin
           ok = dialog_message ("Enter a positive integer for the mask ID.")
           ; fail if the mask ID is not a pos. integer.
           return
        endif  

        ; Input the data directory and calibration files

        widget_control, idList.dataDir,  get_value=dataDir
        widget_control, idlist.flatname, get_value=flatname
        widget_control, idlist.arcName,  get_value=arcname

        ; Input the list of slits. The text may only contain spaces and digits.
        
        Widget_control,  idList.slitNums,  get_value=tempStr

                                ; Split the list, delimited by
                                ; spaces, into an array.

        strArray = STRSPLIT(tempStr, ' ', /EXTRACT) 
           
        
           
                                ; make sure there are only digits and spaces
           
        s =  size(strArray)
          
        if (s[0] gt 0) then begin
           
                                ; s[1] is the size of the array, if
                                ; s[0], the # of dimensions, is nonzero.
           for i=0, s[1]-1 do begin
              
              if (isPosInt(strArray[i]) eq 0) then begin
                 ok = dialog_message("Enter only a list of integers " + $
                                     "separated by spaces for the list " + $
                                     "of slits, or leave it blank.")
                 return
              endif
           endfor
        endif
        
        
        ; Input the nArray values:
        
        nArray = intarr(4)
        
        for i=0, 3 do begin
           widget_control,  idList.nArray[i, 0], get_value=selected
           
                                ; for each element, find out if it has
                                ; been selected:

           if (selected eq 1) then begin
              
                                ; Make sure it's a positive integer:

              Widget_control,  idList.nArray[i, 1], get_value=tempStr
              num = isPosInt (tempStr) 
              
              if (num le 0) then begin
                 ok = dialog_message ('Enter a positive integer for the ' $ 
                                      + idList.valueList[i])
                 return

              endif else begin 
                 nArray[i] = num ; Add the number to the array.
              endelse
              
           endif
        END
        
                                ; print out the desired output to the
                                ; file ####.ql where #### is the maskID
        
        close,  1
        file = strtrim(string(maskId), 2)+ '.ql'
        openw, 1, file

        ; print mask ID
        printf, 1, format= '("MASK: ", I0)', maskID
 
        ; print datadir and calib. files:
        ;if (dataDir ne '')  then printf, 1, "RAWDATADIR: ", dataDir
        printf, 1, 'POLYFLAG  --use polyflag for lambda analysis'
        if (flatname ne '') then printf, 1, "FLATNAME: ", flatname
        if (arcname ne '')  then printf, 1, "ARCNAME: ", arcname
        
        ; print list  of slits
        if (s[0] gt 0) then begin
           printf, 1, "SLITS: ",  strArray
        endif
        
        ; print nArray
        for i=0, 3 do begin
           if (nArray[i] gt 0) then begin
              printf, 1, format='(A, ": ", I0)', $
                (strsplit(idlist.valuelist[i], ' ', /extract))[1], nArray[i] 
           endif
        endfor

        close, 1

        ; print out success message:
        G0 = dialog_message(file + " was generated successfully", $
                            title="Success", /information)

        resetFields,  idList

        ;Widget_Control,  ev.top,  /destroy
        return
     endcase
                                ; this executes the function of the
                                ; 'cancel' button.
     
     'EXIT' : WIDGET_CONTROL, ev.top,  /DESTROY

     'SETDIR':begin

        dir =  DIALOG_PICKFILE(/DIRECTORY)
        if (dir eq '') then return
        WIDGET_CONTROL, idList.dataDir, set_value=dir

     endcase

     'SETARC': begin

        widget_control, idlist.dataDir, get_value=dataDir

        file= dialog_pickfile(/must_exist, path=dataDir,get_path=newPath, $
                              title= 'PLEASE SELECT AN ARCLAMP FILE')
       
        if (newPath ne '') then $
          widget_control, idList.dataDir, set_value=newPath  
        
        widget_control, idList.arcName, set_value=file
     
     endcase

     'SETFLAT': begin

        widget_control, idlist.dataDir, get_value=dataDir
        
        file= dialog_pickfile(/must_exist, path=dataDir, get_path=newpath, $
                              title= 'PLEASE SELECT A FLATLAMP FILE')
        
        if (newPath ne '') then $ 
          widget_control, idList.dataDir, set_value=newPath
       
        widget_control, idList.flatName, set_value=file
     
     endcase
 
     Else:
  ENDCASE
  
END
  
  
PRO slitChooser
 
  values = ['N brightest objects', 'N dimmest objects', $ 
            'N simplest objects', 'N longest slits']

  ; top-level base
  base =  widget_base(/ROW,  Title="Slit Chooser")
  

                                ; create structure to keep track of
                                ; all the widget ID's of the data
                                ; entry boxes

  idList = { mask:0L, dataDir:0L, flatName:0L,  arcName:0L,  slitNums:0L, $ 
             nArray:LONARR(4, 2), valueList:values}

                                ; this sets the uservalue of the base
                                ; to the structure, thereby passing it
                                ; to the event procedure.

  widget_control, base,  set_uvalue=idList
  
  ; create left and right columns

  leftBase =  widget_base(base, /column, /frame)
  rightBase = widget_base(base, /column, /frame)

  ; Mask ID input:
  idList.mask =  CW_FIELD(leftbase, TITLE='Mask: (required) ', xsize=5, $ 
                          uvalue = "?")
  tmp =  WIDGET_BASE(leftbase, ysize=10) ; spacer 

  ; Calibration file input.  Pressing the button fills attempts to
  ; locate the calibration files and puts them in the boxes below.

  tmp =  WIDGET_BUTTON(leftbase, value="Set flatlamp file", $
                       uvalue="SETFLAT")
  idList.flatName = widget_text(leftbase, xsize=30,  uvalue = "?")
  tmp =  WIDGET_BASE(leftbase, ysize=10) ; spacer 

  tmp =  WIDGET_BUTTON(leftbase, value="Set arclamp file", $
                       uvalue="SETARC")
  idList.arcName  = widget_text(leftbase, xsize=30, uvalue = "?")
  tmp =  WIDGET_BASE(leftbase, ysize=10) ; spacer 

  ; Raw data Directory input: pressing the button calls the file
  ; selection dialog, which puts the file name in the text box 
  ; below it.  The box is also editable.
 ; tmp =  WIDGET_BUTTON(leftbase, value="Set raw data directory",  $
   ;                    uvalue="SETDIR")
  tmp =  widget_label(leftbase, value="Current raw data directory:")
  idList.dataDir = WIDGET_TEXT(leftbase,  xsize= 30, $
                               value='', uvalue = "?")


  ; Input the data that specifies which slits to be analyzed:


  ; Optional list of slits
  tmp =  WIDGET_LABEL(rightbase, value="Enter a list of slits, separated by spaces:")
  idList.slitNums =  widget_text(rightbase, /editable,  xsize=30, uvalue = "?")


  ; Create the buttons for "N __est objects"
 
  base0 = widget_base(rightbase, /ROW,  /FRAME)
  idList.nArray[0, 0] = CW_BGROUP(base0,  values[0],  /COLUMN, $
                                  /NONEXCLUSIVE, UVALUE='button0')
  idList.nArray[0, 1] = widget_text(base0, xsize=5, $
                                    UVALUE='text0',  /EDITABLE) 

  base1   = widget_base(rightbase, /ROW,  /FRAME)
  idList.nArray[1, 0] = CW_BGROUP(base1,  values[1],  /COLUMN,  $
                                 /NONEXCLUSIVE,  UVALUE='button1')
  idList.nArray[1, 1] = widget_text(base1, xsize=5,  $
                                    UVALUE='text1',  /EDITABLE)
 
  base2   = widget_base(rightbase, /ROW,  /FRAME)
  idList.nArray[2, 0] = CW_BGROUP(base2,  values[2],  /COLUMN, $ 
                                  /NONEXCLUSIVE,  UVALUE='button2')
  idList.nArray[2, 1] = widget_text(base2, xsize=5, $
                                    UVALUE='text2',  /EDITABLE)
 
  base3   = widget_base(rightbase, /ROW,  /FRAME)
  idList.nArray[3, 0] = CW_BGROUP(base3,  values[3],  /COLUMN, $
                                  /NONEXCLUSIVE,  UVALUE='button3')
  idList.nArray[3, 1] = widget_text(base3, xsize=5, $
                                    UVALUE='text3',  /EDITABLE)


  tmp =  WIDGET_BUTTON(rightbase,  value= 'Generate file',  UVALUE='GO',  xsize=20)
  tmp =  WIDGET_BUTTON(rightbase,  value='Exit', UVALUE='EXIT', xsize=20)

  Widget_control, base,  set_uvalue=idList

  WIDGET_CONTROL,  base,  /REALIZE
  WIDGET_CONTROL,  leftbase,  /REALIZE
  WIDGET_CONTROL,  rightbase,  /REALIZE  
  WIDGET_CONTROL,  base0,  /REALIZE
  WIDGET_CONTROL,  base1,  /REALIZE
  WIDGET_CONTROL,  base2,  /REALIZE
  WIDGET_CONTROL,  base3,  /REALIZE

  XMANAGER,  'slitChooser',  base
  
END








