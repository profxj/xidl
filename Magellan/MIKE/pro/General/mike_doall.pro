;+
; NAME: mike_doall
;     Version 1.1
;
; PURPOSE:
;    Run all of the steps for a single night of data.  You had better
;    know what you are doing!!
;
;
; CALLING SEQUENCE:
;   
;    mike_doall, setup, MIKE=, /NOCR, /CLOBBER, EXTRA=
;
; INPUTS:
;   mike    -  MIKE structure
;   setup   -  Setup ID
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOCR    - Turn off CR rejection routine (mike_objcr)
;  /CLOBBER - Overwrite previous images (mainly the processed image)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_doall, mike, 1L
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Sep-2004 Written by JXP 
;
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_doall, NOCR=nocr, CLOBBER=clobber, _EXTRA=EXTRA, MIKE=mike, $
                SETCRD=setcrd, SKIPGAIN=skipgain

;
  if  N_params() GT 1  then begin 
      print,'Syntax - ' + $
        'mike_allobj, MIKE=, /NOCR, /CLOBBER, EXTRA=, /SKIPGAIN [v1.1]'
      return
  endif 

  if not keyword_set( MIKE ) then begin
      ;; Create structure
      mike_strct, mike

      ;; Setup to 1 for everything
      mike.setup = 1
      
      ;; Set cards
      if keyword_set( SETCRD ) then begin
          print, 'mike_doall: Set your cards now!!  Then hit continue'
          stop
      endif
      
  endif

  ;; Setup
  mike_setup, mike
      
  setup = mike[uniq(mike.setup,sort(mike.setup))].setup
  nset = n_elements(setup)
      
      
  for qq=0L,nset-1 do begin
      ;; Gain
      if not keyword_set( SKIPGAIN) then mike_setgain, mike, setup[qq]

      ;; Flats
      mike_allflat, mike, setup[qq], _EXTRA=EXTRA, CLOBBER=clobber

      ;; Arcs
      mike_allarc, mike, setup[qq], _EXTRA=EXTRA, CLOBBER=clobber

      ;; Slit flat
      mike_slitflat, mike, setup[qq], _EXTRA=EXTRA, CLOBBER=clobber

      gd = where(mike.setup EQ setup[qq] AND mike.type EQ 'OBJ', ngd)
      if ngd EQ 0 then continue

      ;; Extraction
      obj = mike[gd[uniq(mike[gd].obj_id,sort(mike[gd].obj_id))]].obj_id
      nobj = n_elements(obj)

      for jj=0L,nobj-1 do begin
          ;; Gidde up
          mike_allobj, mike, obj[jj], /procall, NOCR=nocr, CLOBBER=clobber
      endfor
  endfor

;  DONE
  print, 'mike_doall: All done! '
  return
end


