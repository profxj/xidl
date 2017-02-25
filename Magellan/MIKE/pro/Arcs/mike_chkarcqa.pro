;+ 
; NAME:
; mike_chkarcqa
;     Version 1.1
;
; PURPOSE:
;  Simple program to bring up the ps files of the QA for a given arc
;
; CALLING SEQUENCE:
;   
;   mike_chkarcqa, mike, index
;
; INPUTS:
;   mike     -  MIKE structure
;   index    -  Index number of the arc
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_chkarcqa, mike, 33L
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Apr-2004 Written by JXP
;
;  Usage:
;
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_chkarcqa, mike, index
;

  if  N_params() LT 2 then begin
      print,'Syntax - ' + $
        'mike_allarc, mike, index [v1.1]'
      return
  endif 

  ;; Arcfile
  setup = mike[index].setup
  arc_fil = mike_getfil('arc_fil', setup, subfil=mike[index].img_root,/name)

  ;; Fitarc
  qafil = mike_getfil('qa_arcfit', setup, SUBFIL=arc_fil) 
  af = findfile(qafil+'*',count=na)
  if na NE 0 then $
    spawn, 'gv -noantialias -seascape -media letter '+af[0] $
    else print, 'mike_chkarcqa:  No file '+qafil

  ;; Arc2dFit
  qafil = mike_getfil('qa_arc2dfit', setup, SUBFIL=arc_fil) 
  af = findfile(qafil+'*',count=na)
  if na NE 0 then $
    spawn, 'gv -noantialias -seascape -media letter '+af[0] $
    else print, 'mike_chkarcqa:  No file '+qafil

  ;; Tracearc
  qafil = mike_getfil('qa_tracearc', setup, SUBFIL=arc_fil) 
  af = findfile(qafil+'*',count=na)
  if na NE 0 then $
    spawn, 'gv -noantialias -media letter '+af[0] $
    else print, 'mike_chkarcqa:  No file '+qafil

  ;; Fit Tracearc
  qafil = mike_getfil('qa_fittrcarc', setup, SUBFIL=arc_fil) 
  af = findfile(qafil+'*',count=na)
  if na NE 0 then $
    spawn, 'gv -noantialias -seascape -media letter '+af[0] $
    else print, 'mike_chkarcqa:  No file '+qafil
                                                             
  print, 'mike_chkarcqa: All done!'

  return
end
