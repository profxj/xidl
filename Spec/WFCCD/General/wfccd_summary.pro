;+ 
; NAME:
; wfccd_summary   
;     Version 1.0
;
; PURPOSE:
;    Creates and outputs a structure for a series of ESI
;    spectroscopic images
;
; CALLING SEQUENCE:
;   
;  wfccd_summary, struct
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  Creates an IDL structure for direct images 
;         -  ASCII file summarizing the structure
;
; OPTIONAL KEYWORDS:
;   MKDIR      - Make directories
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_summary, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_summary, wfccd, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'wfccd_summary, struct, LIST=, NOFILE=, NOLIST= (v1.0)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTFIL ) then outfil = 'wfccd_summ.txt'

  
; Grab all ECH files

  gd = where(wfccd.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

; Open File
  print, 'wfccd_summary: Check the file ', outfil, ' for info'
  close, 17
  openw, 17, outfil

; Loop on Obj

  obj = uniq(wfccd[gd].mask_id, sort(wfccd[gd].mask_id))
  nobj = n_elements(obj)

  for ii=0,nobj-1 do begin
      ;; Slit
      tobj = wfccd[gd[obj[ii]]].mask_id
;      printf, 17, 'wfccd_summary: Checking for calibs for obj = ', tobj

      ;; Flat
      flg_flt = 0
      iflt = where(wfccd.mask_id EQ tobj AND wfccd.type EQ 'FLT' AND $  ; Flat
                   wfccd.flg_anly NE 0, nflt)
      if nflt EQ 0 then begin
          flg_flt = -1
          printf, 17, 'wfccd_summary:  WARNING! No Flats found for obj ', tobj
      endif

      ;; Arcs
      flg_arc = 0
      arcs = where(wfccd.mask_id EQ tobj AND wfccd.type EQ 'ARC' AND $
                   wfccd.flg_anly NE 0, narc)
      if narc EQ 0 then begin
          printf, 17, 'wfccd_summary:  WARNING! No Arcs found for obj ', tobj
          flg_arc = -1
      endif

      ;; Set Obj_id by Object name
      gdobj = where(wfccd.mask_id EQ tobj AND wfccd.flg_anly NE 0 AND $
                    wfccd.type EQ 'OBJ', ngdobj)
      for jj=0L, ngdobj-1 do begin
          ;; Print
          printf, 17, FORMAT='(i2,1x,a12,1x,a2,1x,a6,1x,i5,1x,2i2,1x,a30)', $
            tobj, $
            wfccd[gdobj[jj]].Obj, $
            wfccd[gdobj[jj]].masknm, $
            wfccd[gdobj[jj]].img_root, $
            long(wfccd[gdobj[jj]].exp), $
            flg_flt, flg_arc, $
            wfccd[gdobj[jj]].msk_fil
      endfor
  endfor
  printf, 17, '---------------------------------------------------'

  ;; ALL DONE
  print, 'wfccd_summary: All done!'
  close, 17

  return
end
