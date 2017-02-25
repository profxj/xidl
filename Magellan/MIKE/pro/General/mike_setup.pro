;+ 
; NAME:
; mike_setup   
;     Version 1.1
;
; PURPOSE:
;   Examine the mike structure to determine if a given setup has the
;   requisite calibration files.  It then sets a number of tag names
;   accordingly.   Finally, it gives the object frames a running index
;   based on the object name.
;
;   It is recommended to rerun this program if: (1) you have bombed out of
;   IDL without saving the MIKE structure, (2) you have changed the
;   setup values or calibration files.  I doubt it will ever hurt to
;   rerun.
;
;
; CALLING SEQUENCE:
;   
;  mike_setup, mike
;
; INPUTS:
;   mike --  MIKE structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   OUTFIL -- Name for the ASCII file output (default:
;             'mike_summ.txt')
;
; COMMENTS:
;
; EXAMPLES:
;   mike_setup, mike, OUTFIL=
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   mike_getfil
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;   29-Jan-2003 Polished by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_setup, mike, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_setup, mike, OUTFIL= (v1.1)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTFIL ) then outfil = 'mike_summ.txt'
  
; Grab all ECH files

  gd = where(mike.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

; Open File
  print, 'mike_setup: Check the file ', outfil, ' for info'
  close, 17
  openw, 17, outfil

; Loop on Slit width

  setup = uniq(mike[gd].setup, sort(mike[gd].setup))
  nsetup = n_elements(setup)

  nobj_id = 0L

; Loop on side
  printf, 17, '--------------------------------------------------'
  for qq=1,2L do begin
      nobj_id = 1L
      ;; SIDE
      if qq EQ 1 then begin
          printf, 17, 'mike_setup: Checking BLUE images'
          nm = 'B'
      endif else begin
          printf, 17, 'mike_setup: Checking RED images'
          nm = 'R' 
      endelse

      for ii=0,nsetup-1 do begin
          ;; Slit
          tsetup = mike[gd[setup[ii]]].setup
          if tsetup GE 10 then c_s = strtrim(tsetup,2) $
          else c_s = '0'+strtrim(tsetup,2) 
          printf, 17, 'mike_setup: Checking for calibs for setup = ', tsetup

          ;; Milky Flat
          flg_flt = 0
          mflt = where(mike.setup EQ tsetup AND mike.side EQ qq AND $
                       mike.type EQ 'MFLT' AND $ 
                       mike.flg_anly NE 0, nmflt)
          if nmflt NE 0 then begin
              printf, 17, 'mike_setup:   Found Milky Flats -- ', $
                mike[mflt].img_root
              flg_flt = 1
          endif else begin
              printf, 17, 'mike_setup:  WARNING! No Milky flats ' + $
                'found for setup ', tsetup
          endelse

          ;; Trace Flat
          flg_flt = 0
          tflt = where(mike.setup EQ tsetup AND mike.side EQ qq AND $
                       mike.type EQ 'TFLT' AND $ 
                       mike.flg_anly NE 0, nmflt)
          if nmflt NE 0 then begin
              printf, 17, 'mike_setup:   Found Trace Flats -- ', $
                mike[tflt].img_root
              flg_flt = 1
          endif else begin
              printf, 17, 'mike_setup:  WARNING! No Trace flats ' + $
                'found for setup ', tsetup
              print, 'mike_setup:  WARNING! No Trace flats ' + $
                'found for setup ', tsetup
          endelse

          ;; Arcs
          arcs = where(mike.setup EQ tsetup AND mike.side EQ qq AND $
                       strtrim(mike.type,2) EQ 'ARC' AND $
                       mike.flg_anly NE 0, narc)
          if narc NE 0 then begin
              printf, 17, 'mike_setup:   Found Arc Frames-- (file, arclamp)'
              writecol, 'dum', mike[arcs].img_root, mike[arcs].arclamp, FILNUM=17
          endif else printf, 17, 'mike_setup:  WARNING! No Arcs ' + $
            'found for setup ', tsetup

          ;; STD
          std = where(mike.setup EQ tsetup AND mike.side EQ qq AND $
                      mike.type EQ 'STD' AND mike.flg_anly NE 0, nstd)
          if nstd NE 0 then begin
              printf, 17, 'mike_setup:   Found Standard stars-- ', $
                mike[std].img_root, ' ', mike[std].Obj 
;              mike[std].arc_fil = 'Arcs/ArcECH_'+c_s+'IMG.fits'
          endif else printf, 17, 'mike_setup:  WARNING! No Standards found for setup '

          ;;   Set Arc names for STD stars
          for jj=0L,nstd-1 do begin
              ;; Identify arc image in time
              diff = abs( (mike[std[jj]].date + $
                           mike[std[jj]].exp/2./3600./24.) - $
                          (mike[arcs].date + mike[arcs].exp/2./3600./24.) )
              mn = min(diff, gdarc)
              arc_fil = mike_getfil('arc_fil', $
                                    subfil=mike[arcs[gdarc]].img_root,$
                                    /name)
              mike[std[jj]].arc_fil = arc_fil
              mike[std[jj]].arc_img = $
                mike_getfil('arc_img', subfil=arc_fil,$
                            /name)
              ;; xyoff
              mike[std[jj]].arc_xyoff = mike[arcs[gdarc]].arc_xyoff
          endfor


          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Set Obj_id by Object name
          obj = where(mike.setup EQ tsetup AND mike.flg_anly NE 0 AND $
                      mike.side EQ qq AND $
                      (mike.type EQ 'OBJ' OR mike.type EQ 'STD'), nobj)
          if nobj NE 0 then begin
              ;; Trim
              mike[obj].Obj = strtrim(mike[obj].Obj,2)
              ;; Unique
              uobj = uniq(mike[obj].Obj, sort(mike[obj].Obj))
              for jj=0L, n_elements(uobj)-1 do begin
                  a = where( mike[obj].Obj EQ (mike[obj].Obj)[uobj[jj]], na )
                  mike[obj[a]].obj_id = nobj_id
                  ;; Set Flat Image
                  mike[obj[a]].flat_fil = 'Flats/Flat_'+nm+'_'+c_s+'_M.fits' 
                  ;; Set Object file
;                  mike[obj[a]].obj_fil = 'Extract/Obj_'+mike[obj[a]].img_root
                  mike[obj[a]].obj_fil = $
                    mike_getfil('obj_fil', tsetup, SUBFIL=mike[obj[a]].img_root, /name)
                  ;; Set Final Image
                  mike[obj[a]].img_final = $
                    mike_getfil('fin_fil', tsetup, SUBFIL=mike[obj[a]].img_root, /name)
                  ;; Print
                  for kk=0,na-1 do begin
                      idx = obj[a[kk]]

                      ;; Identify arc image in time
                      diff = abs( (mike[idx].date + $
                                   mike[idx].exp/2./3600./24.) - $
                                  (mike[arcs].date + mike[arcs].exp/2./3600./24.) )
                      mn = min(diff, gdarc)
;                      diff = abs( (mike[idx].ut + $
;                                   mike[idx].exp/2.) - $
;                                  (mike[arcs].ut + mike[arcs].exp/2.) )
;                      b = where( diff GT 21.*3600, nb)
;                      if nb NE 0 then diff[b] = abs(diff[b] - 24.*3600.)
;                      mn = min(diff, gdarc)
                      arc_fil = mike_getfil('arc_fil', $
                                            subfil=mike[arcs[gdarc]].img_root,$
                                            /name)
                      mike[obj[a[kk]]].arc_fil = arc_fil
                      mike[obj[a[kk]]].arc_img = $
                        mike_getfil('arc_img', subfil=arc_fil,$
                                    /name)
                      ;; Set xyoff
                      mike[obj[a[kk]]].arc_xyoff = mike[arcs[gdarc]].arc_xyoff

                      ;; Print
                      printf, 17, FORMAT='(i4,1x,a7,1x,a12,1x,i3,i5,1x,a15,1x,a28)', $
                        obj[a[kk]], $
                        mike[obj[a[kk]]].img_root, $
                        mike[obj[a[kk]]].Obj, $
                        mike[obj[a[kk]]].obj_id, $
                        long(mike[obj[a[kk]]].exp), $
                        mike[obj[a[kk]]].arc_img
                  endfor
                  ;; Increment
                  nobj_id = nobj_id + 1
              endfor
          endif
          printf, 17, '---------------------------------------------------'
      endfor

  endfor

  ;; ALL DONE
  print, 'mike_setup: All done!'
  close, 17

  return
end
