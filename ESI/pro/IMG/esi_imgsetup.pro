;+ 
; NAME:
; esi_imgsetup   
;     Version 1.0
;
; PURPOSE:
;    Creates and outputs a structure for a series of ESI
;    spectroscopic images
;
; CALLING SEQUENCE:
;   
;  esi_imgsetup, struct, LIST=list, MKDIR=mkdir, 
;               NOFILE=nofile, NOLIST=nolist
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
;   LIST       - Image list
;   MKDIR      - Make directories
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_imgsetup, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_imgsetup, esi, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_imgsetup, esi, OUTFIL= (v1.0)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTFIL ) then outfil = 'esiIMG_summ.txt'

  
; Grab all IMG files

  gd = where(esi.mode EQ 0 AND esi.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

; Open File
  print, 'esi_imgsetup: Check the file ', outfil, ' for info'
  close, 17
  openw, 17, outfil

; Loop on Filter width

  filter = uniq(esi[gd].imfilt, sort(esi[gd].imfilt))
  nfilt = n_elements(filter)

  nobj_id = 0L

  for ii=0,nfilt-1 do begin
      ;; Filter
      tfilt = esi[gd[filter[ii]]].imfilt
      printf, 17, 'esi_imgsetup: Checking for calibs for Filter = ', tfilt

      ;; Flat
      flg_flt = 0
      dflt = where(esi.imfilt EQ tfilt AND esi.type EQ 'DFLT' AND $  ; Internal
                   esi.flg_anly NE 0 AND esi.mode EQ 0, ndflt)
      if ndflt NE 0 then begin
          printf, 17, 'esi_imgsetup:   Found Dome Flats -- ', esi[dflt].img_root
          flg_flt = 1
      endif
      twflt = where(esi.imfilt EQ tfilt AND esi.type EQ 'TWI' AND $  ; Dome
                   esi.flg_anly NE 0 AND esi.mode EQ 0, ntwflt)
      if ntwflt NE 0 then begin
          printf, 17, 'esi_imgsetup:   Found Twilight Flats -- ', esi[twflt].img_root
          flg_flt = 2
      endif
      if flg_flt EQ 0 then $
        printf, 17, 'esi_imgsetup:  WARNING! No Flats found for filter ', tfilt

      ;; Standards
      flg_std = 0
      std = where(esi.imfilt EQ tfilt AND esi.type EQ 'STD' AND $  
                  esi.flg_anly NE 0 AND esi.mode EQ 0, nstd)
      if nstd NE 0 then begin
          printf, 17, 'esi_imgsetup:   Found Standard Stars -- ', esi[std].img_root
          flg_std = 1
      endif
      if flg_std EQ 0 then $
        printf, 17, 'esi_imgsetup:  WARNING! No Standards found for filter ', tfilt

      ;; Set Obj_id by Object name
      obj = where(esi.imfilt EQ tfilt AND esi.flg_anly NE 0 AND $
                  (esi.type EQ 'OBJ' OR esi.type EQ 'STD'), nobj)
      if nobj NE 0 then begin
          uobj = uniq(esi[obj].Obj, sort(esi[obj].Obj))
          for jj=0L, n_elements(uobj)-1 do begin
              a = where( esi[obj].Obj EQ (esi[obj].Obj)[uobj[jj]], na )
              esi[obj[a]].obj_id = nobj_id
              ;; Set Flat name (default to Dome)
              if flg_flt EQ 2 then $
                esi[obj[a]].flat_fil = 'Flats/FlatIMG_T'+tfilt+'.fits'
              ;; Print
              for kk=0,na-1 do begin
                  printf, 17, FORMAT='(i4,1x,a7,1x,a12,1x,i3,a4,1x,a2,1x,i5,1x,a16)', $
                    obj[a[kk]], $
                    esi[obj[a[kk]]].img_root, $
                    esi[obj[a[kk]]].Obj, $
                    esi[obj[a[kk]]].obj_id, $
                    esi[obj[a[kk]]].type, $
                    esi[obj[a[kk]]].imfilt, $
                    long(esi[obj[a[kk]]].exp), $
                    esi[obj[a[kk]]].flat_fil
              endfor
              ;; Increment
              nobj_id = nobj_id + 1
          endfor
      endif
      printf, 17, '---------------------------------------------------'

  endfor

  ;; ALL DONE
  print, 'esi_imgsetup: All done!'
  close, 17

  return
end
