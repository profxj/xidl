;+ 
; NAME:
; mike_mkarc   
;     Version 1.1
;
; PURPOSE:
;    Process and combine arc files  
;
; CALLING SEQUENCE:
;   
;  mike_mkarc, mike, slit, /CLOBBER
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Integer defining setup
;
; RETURNS:
;
; OUTPUTS:
;  One processed, combined Arc Image (e.g. Arcs/Arc_ECH##.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_mkarc, mike, 0.5
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_mkarc, mike, setup, side, CLOBBER=clobber, FLATFIL=flatfil, $
                SVFIN=svfin

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_mkarc, mike, setup, [side], FLATFIL=, /CLOBBER [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]  ; default to red side for now
  
; Setup
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 

  ;; QA
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

; Loop on side
  nside = n_elements(side)
  for ii=0L,nside-1 do begin

      qq = side[ii]
      ;; SIDE
      if qq EQ 1 then begin
          print, 'mike_mkarc: Making BLUE arc'
          nm = 'B'
      endif else begin
          print, 'mike_mkarc: Making RED arc'
          nm = 'R' 
      endelse

      ;; Grab all Arc files
      arcs = where(mike.setup EQ setup AND mike.flg_anly NE 0 AND $
                   mike.side EQ side[ii] AND strtrim(mike.type,2) EQ 'ARC', narc)
      if narc EQ 0 then begin
          print, 'mike_mkarc: No Arcs found! Returning' 
          return
      endif

      ;; Check for prior image
      outfil = 'Arcs/Arc_'+nm+'_'+c_s+'.fits'
      a = findfile(outfil, count=na)
      if na NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'mike_mkarc: Arc ', outfil, ' exists.  Returning'
          return
      endif

      ;; Process
      mike_proc, mike, arcs, FLATFIL=flatfil, CLOBBER=clobber, $
        REDOOV=clobber, /ARC

      ;; Combine
      case narc of 
          1: begin              ; 1 Image
              img_fil = mike[arcs[0]].img_final
              img = xmrdfits(img_fil, 0, head, /silent)
              ivar = xmrdfits(img_fil, 1, head, /silent)
          end
          
          2: begin              ; 2 Images
              head = xheadfits(mike[arcs[0]].img_final, /silent)
              img = mike_addtwo(mike, arcs, /SCALE, /ARC, IVAR=ivar)
          end

          else: begin           ; Muliple means median scaled by exp
              ;; IMG
              xcombine, mike[arcs].img_final, img, FCOMB=2, $
                SCALE=mike[arcs].exp, $
                GAIN=mike[arcs[0]].gain, RN=mike[arcs[0]].readno
              ;; VAR
              xcombine, mike[arcs].img_final, ivar, FCOMB=2, $
                SCALE=mike[arcs].exp, $
                GAIN=mike[arcs[0]].gain, RN=mike[arcs[0]].readno,IMGINDX=1L
          end
      endcase
  
      ;; Output
      print, 'mike_mkarc:  Writing ', outfil
      mwrfits, img, outfil, head, /create, /silent
      mwrfits, ivar, outfil, /silent
      spawn, 'gzip -f '+outfil
      
      ;; Cards
      objstd = where(mike.side EQ side[ii] AND mike.flg_anly NE 0 AND $
                     mike.setup EQ setup AND $
                     (mike.type EQ 'STD' OR mike.type EQ 'OBJ'), nobj)
      if nobj NE 0 then mike[objstd].arc_fil = outfil

      ;; DEL FIN
      if not keyword_set( SVFIN ) then $
        mike_delfin, mike, arcs, /silent

  endfor
;  ALL DONE
  print, 'mike_mkarc: All Done! '
  return

end

