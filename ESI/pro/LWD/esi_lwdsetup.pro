;+ 
; NAME:
; esi_lwdsetup   
;     Version 1.0
;
; PURPOSE:
;    Creates and outputs a structure for a series of ESI
;    spectroscopic images
;
; CALLING SEQUENCE:
;   
;  esi_lwdsetup, struct, LIST=list, MKDIR=mkdir, 
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
;   esi_lwdsetup, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdsetup, esi, LIST=list, MKDIR=mkdir, NOFILE=nofile, $
                   NOLIST=nolist, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_lwdsetup, struct, LIST=, MKDIR=, NOFILE=, NOLIST= (v1.0)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTFIL ) then outfil = 'esiLWD_summ.txt'

  
; Grab all LWD files

  gd = where(esi.mode EQ 1 AND esi.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

; Open File
  print, 'esi_lwdsetup: Check the file ', outfil, ' for info'
  close, 17
  openw, 17, outfil

; Loop on Slit width

  slits = uniq(esi[gd].slit, sort(esi[gd].slit))
  nslit = n_elements(slits)

  nobj_id = 0L

  for ii=0,nslit-1 do begin
      ;; Slit
      tslit = esi[gd[slits[ii]]].slit
      printf, 17, 'esi_lwdsetup: Checking for calibs for slit = ', tslit

      ;; Slit name
      c_s = esi_slitnm( tslit )

      ;; Flat
      flg_flt = 0
      iflt = where(esi.slit EQ tslit AND esi.type EQ 'IFLT' AND $  ; Internal
                   esi.flg_anly NE 0 AND esi.mode EQ 1, niflt)
      if niflt NE 0 then begin
          printf, 17, 'esi_lwdsetup:   Found Internal Flats -- ', esi[iflt].img_root
          flg_flt = 1
      endif
      dflt = where(esi.slit EQ tslit AND esi.type EQ 'DFLT' AND $  ; Dome
                   esi.flg_anly NE 0 AND esi.mode EQ 1, ndflt)
      if ndflt NE 0 then begin
          printf, 17, 'esi_lwdsetup:   Found Dome Flats -- ', esi[dflt].img_root
          flg_flt = 2
      endif
      if flg_flt EQ 0 then $
        printf, 17, 'esi_lwdsetup:  WARNING! No Flats found for slit ', tslit

      ;; Arcs
      arcs = where(esi.slit EQ tslit AND esi.arclamp GT 0 AND $
                   esi.flg_anly NE 0 AND esi.mode EQ 1, narc)
      if narc NE 0 then $
        printf, 17, 'esi_lwdsetup:   Found Arc Frames-- ', esi[arcs].img_root $
      else printf, 17, 'esi_lwdsetup:  WARNING! No Arcs found for slit ', tslit

      ;; STD
      std = where(esi.slit EQ tslit AND esi.type EQ 'STD' AND $
                   esi.flg_anly NE 0 AND esi.mode EQ 1, nstd)
      if nstd NE 0 then begin
          printf, 17, 'esi_echsetup:   Found Standard stars-- '
          for i=0L,n_elements(std)-1 do $
            printf, 17, '            ', esi[std[i]].img_root, ' ', esi[std[i]].Obj
          esi[std].arc_fil = 'Arcs/AIMG_LWD'+c_s+'.fits'
          esi[std].flat_fil = 'Flats/FlatLWD_'+c_s+'.fits'
      endif else printf, 17, 'esi_echsetup:  WARNING! No Standards found for slit '

      printf, 17, ' '
      ;; Set Obj_id by Object name
      obj = where(esi.slit EQ tslit AND esi.flg_anly NE 0 AND $
                  esi.mode EQ 1 AND esi.type EQ 'OBJ', nobj)
      if nobj NE 0 then begin
          uobj = uniq(esi[obj].Obj, sort(esi[obj].Obj))
          for jj=0L, n_elements(uobj)-1 do begin
              a = where( esi[obj].Obj EQ (esi[obj].Obj)[uobj[jj]], na )
              esi[obj[a]].obj_id = nobj_id
              ;; Set Arc Image
              esi[obj[a]].arc_fil = 'Arcs/AIMG_LWD'+c_s+'.fits'
              ;; Set Flat name (default to Dome)
              if flg_flt EQ 2 then $
                esi[obj[a]].flat_fil = 'Flats/FlatLWD_'+c_s+'.fits'
              ;; Print
              for kk=0,na-1 do begin
                  printf, 17, FORMAT='(i4,1x,a7,1x,a12,1x,i3,i5,1x,a15,1x,a16)', $
                    obj[a[kk]], $
                    esi[obj[a[kk]]].img_root, $
                    esi[obj[a[kk]]].Obj, $
                    esi[obj[a[kk]]].obj_id, $
                    long(esi[obj[a[kk]]].exp), $
                    esi[obj[a[kk]]].arc_fil, $
                    esi[obj[a[kk]]].flat_fil
              endfor
              ;; Increment
              nobj_id = nobj_id + 1
          endfor
      endif
      printf, 17, '---------------------------------------------------'

  endfor

  ;; ALL DONE
  print, 'esi_lwdsetup: All done!'
  close, 17

  return
end
