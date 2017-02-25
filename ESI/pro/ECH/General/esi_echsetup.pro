;+ 
; NAME:
; esi_echsetup   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of ESI
;    spectroscopic images
;
; CALLING SEQUENCE:
;   
;  esi_echsetup, esi
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   esi --  ESI structure
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echsetup, esi, OUTFIL=
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;   29-Jan-2003 Polished by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echsetup, esi, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echsetup, esi, OUTFIL= (v1.1)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTFIL ) then outfil = 'esiECH_summ.txt'
  

; Grab all ECH files

  allg = where(esi.mode EQ 2 AND esi.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

; Open File
  print, 'esi_echsetup: Check the file ', outfil, ' for info'
  close, 17
  openw, 17, outfil

; Loop on binning
  bintot = esi[allg].cbin + 10*esi[allg].rbin
  bintyp = bintot[uniq(bintot, sort(bintot))]
  nbin = n_elements(bintyp)

  for qq=0L,nbin-1 do begin
      rbin = bintyp[qq]/10
      cbin = bintyp[qq] - rbin*10

      gd = where(esi.mode EQ 2 AND esi.flg_anly NE 0 AND $
                 esi.rbin EQ rbin AND esi.cbin EQ cbin, ngd)

      ;; Loop on Slit width

      slits = uniq(esi[gd].slit, sort(esi[gd].slit))
      nslit = n_elements(slits)
      
      nobj_id = 0L

      for ii=0,nslit-1 do begin
          ;; Slit
          tslit = esi[gd[slits[ii]]].slit
          if tslit EQ 9.99 then continue
          c_s = esi_slitnm(tslit)
          printf, 17, 'esi_echsetup: Checking for calibs for slit = ', tslit
          
          ;; Flat
          flg_flt = 0
          iflt = where(esi.slit EQ tslit AND esi.type EQ 'IFLT' AND $ ; Internal
                       esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                       esi.flg_anly NE 0 AND esi.mode EQ 2, niflt)
          if niflt NE 0 then begin
              printf, 17, 'esi_echsetup:   Found Internal Flats -- ', $
                esi[iflt].img_root
              flg_flt = 1
          endif
          dflt = where(esi.slit EQ tslit AND esi.type EQ 'DFLT' AND $ ; Dome
                       esi.flg_anly NE 0 AND esi.mode EQ 2, ndflt)
          if ndflt NE 0 then begin
              printf, 17, 'esi_echsetup:   Found Dome Flats -- ', esi[dflt].img_root
              flg_flt = 2
           endif
          tflt =  where(esi.slit EQ tslit AND esi.type EQ 'TWI' AND $ ;TWI
                        esi.flg_anly NE 0 AND esi.mode EQ 2, ntwi)
          if ntwi NE 0 then begin
             printf, 17, 'esi_echsetup:   Found Twilight Flats -- ' $
                     , esi[tflt].img_root
          endif
          if flg_flt EQ 0 then $
            printf, 17, 'esi_echsetup:  WARNING! No Flats found for slit ' $
                    , tslit
          
          ;; Arcs
          arcs = where(esi.slit EQ tslit AND esi.arclamp GT 0 AND $
                       esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                       esi.flg_anly NE 0 AND esi.mode EQ 2, narc)
          if narc NE 0 then begin
              printf, 17, 'esi_echsetup:   Found Arc Frames-- (file, arclamp)'
              writecol, 'dum', esi[arcs].img_root, esi[arcs].arclamp, FILNUM=17
          endif else $
            printf, 17, 'esi_echsetup:  WARNING! No Arcs found for slit ', tslit

          ;; STD
          std = where(esi.slit EQ tslit AND esi.type EQ 'STD' AND $
                      esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                      esi.flg_anly NE 0 AND esi.mode EQ 2, nstd)
          if nstd NE 0 then begin
              printf, 17, 'esi_echsetup:   Found Standard stars-- ', $
                esi[std].img_root, ' ', esi[std].Obj 
              esi[std].arc_fil = esi_getfil('arc_img', SLIT=tslit, $
                                            cbin=cbin, rbin=rbin, /name)
              esi[std].flat_fil = esi_getfil('finflat_fil', SLIT=tslit, $
                                             cbin=cbin, rbin=rbin, /name)
              IF ntwi GT 0 THEN $
                 esi[std].twiflat_fil =  esi_getfil('twiflat_fil' $
                                                    , SLIT = tslit $
                                                    , cbin = cbin $
                                                    , rbin = rbin, /name) $
              ELSE  esi[std].twiflat_fil =  esi[std].flat_fil
              esi[std].obj_fil = 'Extract/Obj_'+esi[std].img_root
              esi[std].img_final = 'Final/f_'+esi[std].img_root
          endif else $
            printf, 17, 'esi_echsetup:  WARNING! No Standards found for slit '

          ;; Set Obj_id by Object name
          obj = where(esi.slit EQ tslit AND esi.flg_anly NE 0 AND $
                      esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                      esi.mode EQ 2 AND esi.type EQ 'OBJ', nobj)
          if nobj NE 0 then begin
              uobj = uniq(esi[obj].Obj, sort(esi[obj].Obj))
              for jj=0L, n_elements(uobj)-1 do begin
                  a = where( esi[obj].Obj EQ (esi[obj].Obj)[uobj[jj]], na )
                  esi[obj[a]].obj_id = nobj_id
                  ;; Slit
                  c_s = esi_slitnm( tslit )
                  ;; Set Arc Image
                  esi[obj[a]].arc_fil = esi_getfil('arc_img', SLIT=tslit, $
                                                   cbin=cbin, rbin=rbin, /name)
                  ;; Set Flat Image
                  esi[obj[a]].flat_fil = esi_getfil('finflat_fil' $
                                                    , SLIT = tslit $
                                                    , cbin = cbin, rbin = rbin $
                                                    , /name)
                  ;; If we have twilights use them, otherwise use illum 
                  ;; function from domeflat
                  IF ntwi GT 0 THEN $
                    esi[obj[a]].twiflat_fil =  esi_getfil('twiflat_fil' $
                                                          , SLIT = tslit $
                                                          , cbin = cbin $
                                                          , rbin = rbin $
                                                          , /name) $
                  ELSE  esi[obj[a]].twiflat_fil =  esi[obj[a]].flat_fil
                  ;; Set Object file
                  esi[obj[a]].obj_fil = 'Extract/Obj_'+esi[obj[a]].img_root
                  ;; Set Final Image
                  esi[obj[a]].img_final = 'Final/f_'+esi[obj[a]].img_root
                  ;; Print
                  for kk = 0, na-1 do begin
                      printf, 17, $
                        FORMAT='(i4,1x,a12,1x,a12,1x,i3,i5,1x,a15,1x,a17)', $
                        obj[a[kk]], $
                        esi[obj[a[kk]]].img_root, $
                        esi[obj[a[kk]]].Obj, $
                        esi[obj[a[kk]]].obj_id, $
                        long(esi[obj[a[kk]]].exp), $
                        esi[obj[a[kk]]].arc_fil
                  endfor
                  ;; Increment
                  nobj_id = nobj_id + 1
              endfor
          endif
          printf, 17, '---------------------------------------------------'
      endfor
  endfor

  ;; ALL DONE
  print, 'esi_echsetup: All done!'
  close, 17

  return
end
