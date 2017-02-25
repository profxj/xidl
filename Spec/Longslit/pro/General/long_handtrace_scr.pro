pro long_handtrace_scr, tracelist=tracelist, planfil=planfil
;+
; Script to carry out 'rescue' reductions by feeding information to
; LONG_REDUCE regarding sources that we want to trace, that did not
; automatically get traced.
;
; First read in TRACELIST which has the extra sources to be
; traced, and then run LONG_REDUCE with the HAND_X and HAND_Y
; information alongside the planfile specified in planfil for just 1 exposure
;-

  
  if not keyword_set(tracelist) then begin 
     print, 'Error: Enter list of objects to be traced, ' + $
            'in columns of slit number, HAND_X, HAND_Y, HAND_FWHM'
     stop
  endif 

  
readcol, tracelist, slitnum, hand_x, hand_y, hand_fwhm, $
         f='l,l,l,f', comment='#',/silent
long_reduce, planfil, hand_x=hand_x, hand_y=hand_y, $
             hand_fwhm=hand_fwhm

end
