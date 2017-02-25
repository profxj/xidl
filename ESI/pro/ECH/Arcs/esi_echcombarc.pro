;+ 
; NAME:
; esi_echcombarc   
;     Version 1.0
;
; PURPOSE:
;    Process arc file
;
; CALLING SEQUENCE:
;   
;  esi_echcombarc, esi, slit
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echcombarc, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echcombarc, esi, indx, FLATFIL=flatfil, BIASFIL=biasfil, CLOBBER = CLOBBER

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echcombarc, esi, indx [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  slit = esi[indx[0]].slit
  flg_lmp = strtrim(esi[indx[0]].arclamp,2)
  c_s = esi_slitnm(slit)

  nimg = n_elements(indx)

  out_fil = 'Arcs/ATMP_'+c_s+'_'+flg_lmp+'.fits'
  a = findfile(out_fil, count=na)
  if na NE 0 AND not keyword_set( CLOBBER ) then begin
      print, 'esi_echcombarc: Arc ', out_fil, ' exists.  Returning'
      return
  endif
  
  ;; Process
  esi_echproc, esi, indx, FLATFIL=flatfil, CLOBBER = CLOBBER, BIASFIL=biasfil

  case nimg of 
      1: begin ; 1 Image
          img_fil = esi[indx[0]].img_final
          img = xmrdfits(img_fil, 0, head, /silent)
          var = xmrdfits(img_fil, 1, /silent)
      end

      2: begin ; 2 Images
          head = xheadfits(esi[indx[0]].img_final, /silent)
          img = esi_addtwo(esi, indx, /SCALE, /ARC, VAR=var)
      end

      else: begin ; Muliple means median scaled by exp
          ;; IMG
          xcombine, esi[indx].img_final, img, FCOMB = 2, $
            SCALE=esi[indx].exp, $
            GAIN=esi[indx[0]].gain, RN=esi[indx[0]].readno
          ;; VAR
          xcombine, esi[indx].img_final, var, FCOMB=2, SCALE=esi[indx].exp, $
            GAIN=esi[indx[0]].gain, RN=esi[indx[0]].readno,IMGINDX=1L
      end
  endcase
  
  ;; Output
  mwrfits, img, out_fil, head, /create, /silent
  mwrfits, var, out_fil, /silent
  delvarx, img, var

  ;; DEL Final
  esi_delfin, esi, indx
  

  return
end
