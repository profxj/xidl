;+ 
; NAME:
; kast_specsetup   
;     Version 1.1
;
; PURPOSE:
;    Examines the Kast IDL structure, checks for appropriate
;  Flats and Arcs and then outputs a summary ASCII file.
;
; CALLING SEQUENCE:
;  kast_specsetup, kast, OUTFIL=
;
; INPUTS:
;   kast --  Kast IDL structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  OUTFIL=  -- Name of ASCII file (default: kast_specsumm.txt')
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_specsetup, kast
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_specsetup, kast, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_specsetup, kast, OUTFIL= [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTFIL ) then outfil = 'kast_specsumm.txt'
  qfil = ' '
  afil = ' '
  
; Grab all spec files

  gd = where(kast.mode EQ 1 AND kast.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

; Open File
  print, 'kast_specsetup: Check the file ', outfil, ' for info'
  close, 17
  openw, 17, outfil

; Loop on side

  for qq=1,2 do begin
      setup = -1L
      nobj_id = 0L
      ;; Label
      case qq of 
          1: begin
              printf, 17, 'BLUE Setups: -------------------------' 
              lbl = 'B'
          end
          2: begin
              printf, 17, 'RED Setups: -------------------------' 
              lbl = 'R'
          end
          else:
      endcase

      clr = where(kast[gd].side EQ qq, nclr)
      if nclr EQ 0 then begin
          printf, 17, 'No frames.  Continuing...'
          continue
      endif

; Loop on Slit width
      slits = uniq(kast[gd[clr]].slit, sort(kast[gd[clr]].slit))
      nslit = n_elements(slits)

      for ii=0,nslit-1 do begin
          ;; Slit
          tslit = kast[gd[clr[slits[ii]]]].slit
          
          all = where(kast[gd[clr]].slit EQ tslit)
          all = gd[clr[all]]

; Loop on Grating/Grism
          grsng = uniq(kast[all].grising, sort(kast[all].grising))
          ngrsng = n_elements(grsng)

          for jj=0,ngrsng-1 do begin
              tgrsng = kast[all[grsng[jj]]].grising
              mall = where(kast[all].grising EQ tgrsng)
              mall = all[mall]

; Beam Splitter
              bsplit = uniq(kast[mall].splitter, sort(kast[mall].splitter))
              nbsplit = n_elements(bsplit)
              for ll=0,nbsplit - 1 do begin
                  ;; SETUP
                  tbsplit = kast[mall[bsplit[ll]]].splitter
                  tot = where(kast[mall].splitter EQ tbsplit)
                  tot = mall[tot]
                  setup = setup + 1L
                  kast[tot].setup = setup

                  printf, 17, '----------------------------'
                  printf, 17, 'Setup: ', setup
                  printf, 17, 'kast_specsetup: Slit = ', tslit
                  printf, 17, 'kast_specsetup: Grating/Grism = ', tgrsng
                  printf, 17, 'kast_specsetup: Splitter = ', tbsplit

                  ;; QTZ
                  qtz = where(kast[tot].type EQ 'QTZ',nqtz)
                  if nqtz NE 0 then begin
                      printf, 17, 'kast_specsetup:   Found Quartz -- ', $
                        kast[tot[qtz]].img_root
                      qfil = 'Flats/Flat'+lbl+strtrim(setup,2)+'.fits'
                  endif else $
                    printf, 17, $
                    'kast_specsetup:  WARNING! No Flats found for grising ', tgrsng

                  ;; Arcs
                  arcs = where(kast[tot].type EQ 'ARC',narc)
                  if narc NE 0 then begin
                      printf, 17, 'kast_specsetup:   Found Arc Frames-- (file, arclamp)'
                      writecol, 'dum', kast[tot[arcs]].img_root, $
                        kast[tot[arcs]].arclamp, FILNUM=17
                      afil = 'Arcs/Arc'+lbl+strtrim(setup,2)+'.fits'
                  endif else printf, 17, $
                'kast_specsetup:  WARNING! No Arcs found '

                  ;; STD
                  std = where(kast[tot].type EQ 'STD',nstd)
                  if nstd NE 0 then begin
                      printf, 17, 'kast_specsetup:   Found Standard stars-- ', $
                        kast[tot[std]].img_root, ' ', kast[tot[std]].Obj 
                      kast[tot[std]].arc_fil = afil
                      kast[tot[std]].flat_fil = qfil
                  endif else $
                printf, 17, 'kast_specsetup:  WARNING! No Standards found for slit '

                  ;; Set Obj_id by Object name
                  obj = where(kast[tot].type EQ 'OBJ',nobj)
                  if nobj NE 0 then begin
                      printf, 17, '  Objects: ============'
                      uobj = uniq(kast[tot[obj]].Obj, sort(kast[tot[obj]].Obj))
                      for mm=0L, n_elements(uobj)-1 do begin
                          a = where( kast[tot[obj]].Obj EQ $
                                     (kast[tot[obj]].Obj)[uobj[mm]], na )
                          kast[tot[obj[a]]].obj_id = nobj_id
                          ;; Set Arc Image
                          kast[tot[obj[a]]].arc_fil = afil
                          ;; Set Flat Image
                          kast[tot[obj[a]]].flat_fil = qfil
                          ;; Set Object file
                          kast[tot[obj[a]]].obj_fil = $
                            strtrim('Extract/Obj_'+kast[tot[obj[a]]].img_root,2)
                          ;; Set Final Image
                          kast[tot[obj[a]]].img_final = $
                            'Final/f_'+kast[tot[obj[a]]].img_root
                          ;; Print
                          for kk=0,na-1 do begin
                              printf, 17, $
                                FORMAT='(i4,1x,a8,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
                                tot[obj[a[kk]]], $
                                kast[tot[obj[a[kk]]]].img_root, $
                                kast[tot[obj[a[kk]]]].Obj, $
                                kast[tot[obj[a[kk]]]].obj_id, $
                                long(kast[tot[obj[a[kk]]]].exp), $
                                kast[tot[obj[a[kk]]]].arc_fil, $
                                kast[tot[obj[a[kk]]]].flat_fil
                          endfor
                          ;; Increment
                          nobj_id = nobj_id + 1
                      endfor
                      printf, 17, '           ============'
                  endif
                  printf, 17, '---------------------------------------------------'
              endfor
          endfor
      endfor
  endfor

  ;; ALL DONE
  print, 'kast_specsetup: All done!'
  close, 17

  return
end
