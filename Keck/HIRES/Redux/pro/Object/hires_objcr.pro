;+ 
; NAME:
; hires_objcr   
;   Version 1.1
;
; PURPOSE:
;    Flags CRs given 2 or more object images.  Calls x_specobjcr
;
; CALLING SEQUENCE:
;   
;  hires_objcr, hires, setup, obj_id, chip, [iexp]
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1), Green (2) OR Red (3) chip
;   [iexp]  -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  GROW=  -- Radius to grow flagged CRs  (deafult: 2)
;  RTIO=  -- Minimum ratio between 'fiducial' and individual exposure
;            for flagging a CR  [default: 9]
;  LTHRESH=  -- Minimum value of CR for flagging
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Algorithm for 3 or more images is not well tested.
;
; EXAMPLES:
;   hires_objcr, hires, setup, obj_id, chip, 
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_objcr, hires, setup, obj_id, chip, iexp, CHK=chk, RTIO=rtio, $
                 LTHRESH=lthresh, GROW=grow

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_objcr, hires, setup, obj_id, chip, [exp], /CHK, RTIO= [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set(RTIO) then rtio = 9.
  if not keyword_set(LTHRESH) then lthresh = 80.
  if not keyword_set(CHIP) then chip = [1L,2L,3L]
  if n_elements(GROW) EQ 0 then grow = 2L 

  print, 'hires_objcr: Growing by ', grow, ' pixels'
  
  for ii=0L,n_elements(chip)-1 do begin
      qq = chip[ii]
      case qq of
          -1: print, 'hires_objcr: Single chip'
          1: print, 'hires_objcr: BLUE chip' 
          2: print, 'hires_objcr: GREEN chip' 
          3: print, 'hires_objcr: RED chip' 
          else:stop
      endcase
      ;;  Find the frames
      exp = where(hires.type EQ 'OBJ' AND hires.flg_anly NE 0 AND $
                   hires.setup EQ setup AND hires.chip EQ qq AND $
                   hires.obj_id EQ obj_id, nindx)
      if keyword_set( IEXP ) then exp = exp[iexp]
 
      nindx = n_elements(exp)

      if nindx LT 2 then begin
          print, 'hires_objcr: Less than 2 images!'
          continue
      endif

      ;; Call x_specobjcr
      x_specobjcr, hires[exp].img_final, EXPT=hires[exp].exp, $
        RTIO=rtio, LTHRESH=lthresh, CHK=chk, GROW=grow

  endfor

  if not keyword_set( SILENT ) then print, 'hires_objcr: All done!'
  return
end
