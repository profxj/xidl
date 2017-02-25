;+ 
; NAME:
; mike_chktrc   
;     Version 1.0
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  mike_chktrc, mike, setup, obj_id, side, [exp]
;
; INPUTS:
;   mike     -  ESI structure
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp]   -  Exposure frames (e.g. [0L, 1L]); Default = 0L
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STD     - Find object for standard star
;  /CHK     - Show overall trace
;  FITFIL=  - Map of pinholes (default: Maps/hole_fit.idl )
;  REFWV=   - Reference wavelength region (default: [5600., 5610.] )
;  REFORDR= - Reference order  (default: 4L)
;  /INTER   - Interactively identify objects and set apertures
;             (recommended)
;  MAXOFF=  - Minimum offset between max obj and center of slit
;             (default: 20.)
;  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_fndobj, mike, 0L, [0L, 1L], /CHK, /INTER, 
;      REFWV=[6500., 6520.], REFORDR=5L
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_chktrc, mike, setup, obj_id, side, exp, NOSTOP=NOSTOP

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_chktrc, mike, setup, obj_id, side, [exp], /NOSTOP [v1.0]'
      return
  endif 
  
  indx = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
               mike.side EQ side  AND $
               mike.obj_id EQ obj_id AND strtrim(mike.type,2) EQ 'OBJ', $
               nindx)

  if keyword_set( EXP ) then indx = indx[exp] else indx = indx[0]

;  Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  nordr = n_elements(ordr_str)

  ;; Read IMG+VAR
  imgfil = 'Final/f_'+mike[indx].img_root
  if not keyword_set( NOSTOP ) then begin
      print, 'mike_chktrcflt:  You need to run this command before '+ $
        'running this program:  '
      print, '  IDL>  xatv, ''', imgfil, ''''
      print, 'If you havent done this yet, return out of this program'
      print, '  run the command and rerun this program with /NOSTOP'
      print, 'If you have done this command, just continue'
      stop
  endif

  objfil = mike[indx].obj_fil
  objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
  
  ny = 4096L/mike[indx].rowbin

  for qq=0L,nordr-1 do xatvplot, objstr[qq].trace[0:ny-1], findgen(ny)

;  DONE
;  print, 'mike_fndobj: All done! '
  return
end
