;+ 
; NAME:
; xdimg_starid   
;   Version 1.1
;
; PURPOSE:
;    Helps identify standard stars in a field by driving the program
;     x_starid which is an interactive, GUI program.
;
; CALLING SEQUENCE:
;   
;   xdimg_starid, struct, LST_PATH=, INDX=
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest.  This
;             program focuses on the STD frames only.
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   LST_PATH   - Path name for Lists (default: Take x_starid default)
;   OUT_PATH   - Output name for id files (default: 'Photo/')
;
; OPTIONAL OUTPUTS:
;   INDX -- Subset of indices for standards to fit
;   INORIENT -- Choose a number in the range [-4,4] to change the 
;                   orientation of standards plotted on the image. 
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_starid, strct
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Aug-2001 Written by JXP
;   08-Jan-2002 Added paths (JXP)
;   14-Oct-2002 Added indx (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_starid, struct, LST_PATH=lst_path, INDX=indx, INORIENT=inorient

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'xdimg_starid, struct, LST_PATH=, INDX=, INORIENT= [v1.1]'
      return
  endif 

  close, /all

;  Optional Keywords

  if not keyword_set( OUT_PATH ) then out_path = 'Photo/'

;  Find the Standard Stars

  stds = where(struct.type EQ 'STD' AND struct.flg_anly NE 0 AND $
               struct.flg_final NE 0, nstds)
  if nstds EQ 0 then begin
      print, 'No standard stars to analyse!'
      return
  endif
  if keyword_set(INDX) then stds = stds[indx]

;  Run x_starid

  x_starid, struct[stds].img_final, struct[stds[0]].date, $
    struct[stds[0]].ccd, struct[stds[0]].tel, OBJ=struct[stds].obj, $
    LST_PATH=lst_path, OUT_PATH=out_path, INORIENT=inorient

  ; Resave the updated structure so that you can pick up where you left off easily.
   mwrfits, struct, 'struct.fits', /create


end
