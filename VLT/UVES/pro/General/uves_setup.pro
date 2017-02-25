;+ 
; NAME:
; uves_setup   
;     Version 1.1
;
; PURPOSE:
;    Examines the Kast IDL structure, checks for appropriate
;  Flats and Arcs and then outputs a summary ASCII file.
;
; CALLING SEQUENCE:
;  uves_setup, uves, OUTFIL=
;
; INPUTS:
;   uves --  Kast IDL structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  OUTFIL=  -- Name of ASCII file (default: uves_summ.txt')
;  /KEEPSET -- Dont modify the setup numbers (useful to reset the
;              arcs)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_setup, uves
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro uves_setup_val, uves, setups, ECHANG=echang, XDANG=xdang, ETOLER=etoler,$
                     XTOLER=xtoler

  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.1

  nhir = n_elements(uves)

  ;; ECHANGL
  echang = uniq(uves.echangl, sort(uves.echangl))
  nech = n_elements(echang)
  msk = bytarr(nech)
  msk[*] = 1B
  lidx = lindgen(nech)
  for ll=0L,nech-1 do begin     ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(uves[echang[ll]].echangl-$
                   uves[echang].echangl) LT ETOLER AND $
                   lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  echang = uves[echang[where(msk)]].echangl
  nech = n_elements(echang)

  ;; XDANGL
  xdang = uniq(uves.xdangl, sort(uves.xdangl))
  nxd = n_elements(xdang)
  msk = bytarr(nxd)
  msk[*] = 1B
  lidx = lindgen(nxd)
  for ll=0L,nxd-1 do begin      ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(uves[xdang[ll]].xdangl-$
                       uves[xdang].xdangl) LT XTOLER $
                   AND lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  xdang = uves[xdang[where(msk)]].xdangl
  nxd = n_elements(xdang)

  ;; Setup
  setups = strarr(nhir)
  for qq=0L,nhir-1 do begin
      ;; Angles
      mtch = min(abs(uves[qq].echangl-echang),imn)
      echs = strtrim(imn,2)
      mtch = min(abs(uves[qq].xdangl-xdang),imn)
      xds = strtrim(imn,2)
      setups[qq] = string(round(100*uves[qq].slitwid),format='(i3)')+ $
                   '-'+echs+'-'+xds+'-'+$
                   strtrim(uves[qq].rowbin,2)+'-'+strtrim(uves[qq].colbin,2)
;                   '-'+uves[qq].block+'-'+echs+'-'+xds+'-'+$
  endfor
  setups = setups[uniq(setups, sort(setups))]

  return
end
      

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uves_setup, uves, OUTFIL=outfil, XTOLER=xtoler, ETOLER=etoler, $
                 KEEPSET=keepset

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'uves_setup, uves, OUTFIL=, ETOLER=, XTOLER=, /KEEPSET [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.1
  if not keyword_set( OUTFIL ) then outfil = 'uves_summ.txt'
  if not keyword_set( SUMMFIL ) then summfil = 'setup_summ.txt'

  qtz_fil = ' '
  arc_fil = ' '
  
  ;; Grab all  files
  gd = where(uves.type NE 'ZRO' and uves.type NE 'UNK' and uves.flg_anly NE 0, $
             ngd)
  if ngd EQ 0 then return

  ;; Open File
  print, 'uves_setup: Check the files ', outfil+' '+summfil, ' for info'
  close, 17
  openw, 17, outfil
  close, 19
  openw, 19, summfil


  setup = 0L
  nobj_id = 0L

  ;; Loop on Decker
  slitwid = uniq(uves[gd].slitwid, sort(uves[gd].slitwid))
  nslitwid = n_elements(slitwid)
;  bd = where(uves[gd[slitwid]].slitwid EQ 'D5',nbd)
;  if nbd NE 0 then stop

  uves_setup_val, uves[gd], setups, ECHANG=echang, XDANG=xdang, $
                   ETOLER=etoler, XTOLER=xtoler 
  nset = n_elements(setups)

  if keyword_set(KEEPSET) then begin
      numset = uves[uniq(uves[gd].setup, sort(uves[gd].setup))].setup
      nset = n_elements(numset)
  endif

  ;; Summary
  printf, 19, 'Setup SlitW BLKF  ECHANG   XDANG  RBIN CBIN'

  objset = 1L
  calibset = 50L
  svset = lonarr(nset)
  for bb=0L,nset-1 do begin
      ;; Keep setups
      if keyword_set(KEEPSET) then begin
          stop
      endif else kk=bb

      ;; Find files
      prs = strsplit(setups[kk], '-', /extract)
      tslitwid = prs[0]
;      tblock = prs[1]
      sech  = prs[1]
      sxd   = prs[2]
      rbin = prs[3]
      cbin = prs[4]

      mtch = where( round(100*uves[gd].slitwid) EQ round(long(tslitwid)) AND $
;                    uves[gd].block EQ tblock AND $
                    abs(uves[gd].echangl-echang[long(sech)]) LT etoler AND $
                    abs(uves[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                    uves[gd].rowbin EQ long(rbin) AND $
                    uves[gd].colbin EQ long(cbin),nm)
      if nm EQ 0 then begin
          if xdang[long(sxd)] LT 1 then continue else stop
      endif
      xall = gd[mtch]

      ;; Set Obj_id by Object name
      allobj = where(uves[xall].type EQ 'OBJ',naobj)
      if naobj NE 0 then begin
          uaobj = uniq(uves[xall[allobj]].Obj, $
                       sort(uves[xall[allobj]].Obj))
          valuobj = uves[xall[allobj[uaobj]]].Obj
          setup = objset
          objset = objset + 1
      endif else begin
          setup = calibset
          calibset = calibset + 1
      endelse
      if uves[xall[0]].setup NE setup AND uves[xall[0]].setup GT 0 then begin
          print, 'uves_setup:  About to change setup!!!  Think carefully...'
          stop
      endif

      uves[xall].setup = setup
      svset[kk] = setup
      
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, 'Setup: ', setup
      printf, 17, 'uves_setup: SWidth = ', tslitwid
;      printf, 17, 'uves_setup: Blocking = ', tblock
      printf, 17, 'uves_setup: ECH Angle = ', echang[long(sech)]
      printf, 17, 'uves_setup: XDAngle = ', xdang[long(sxd)]
      printf, 17, 'uves_setup: Rowbin = ', rbin
      printf, 17, 'uves_setup: Colbin = ', cbin

;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, 'Setup: ', setup
;      printf, 19, 'uves_setup: Decker = ', tslitwid
;      printf, 19, 'uves_setup: Blocking = ', tblock
;      printf, 19, 'uves_setup: ECH Angle = ', echang[long(sech)]
;      printf, 19, 'uves_setup: XDAngle = ', xdang[long(sxd)]
;      printf, 19, 'uves_setup: Rowbin = ', rbin
;      printf, 19, 'uves_setup: Colbin = ', cbin

      
      ;; Loop on WCEN
      wcen = uniq(uves[xall].xdangl, sort(uves[xall].xdangl))
      nwcen = n_elements(wcen)
      for pp=0L,nwcen-1 do begin
          twcen = round(uves[xall[wcen[pp]]].xdangl)
          call = where(round(uves[xall].xdangl) EQ twcen)
          call = xall[call]
          
          printf, 17, '-----'
          printf, 17, 'uves_setup: Chip (BGR) = ', twcen
          
          ;; QTZ
          qtz = where(uves[call].type EQ 'TFLT',nqtz)
          if nqtz NE 0 then begin
              printf, 17, 'uves_setup:   Found Trace Flats -- ', $
                uves[call[qtz]].img_root
              qtz_fil = uves_getfil('qtz_fil', $
                                     setup, WCEN=twcen, /name)
          endif else $
            printf, 17, $
            'uves_setup:  WARNING! No Trace Flats found..' 
          
          ;; Arcs
          arcs = where(uves[call].type EQ 'ARC',narc)
          if narc NE 0 then begin
              printf, 17, 'uves_setup:   Found Arc Frames-- ' + $
                '(file, arclamp)'
              writecol, 'dum', uves[call[arcs]].img_root, $
                replicate(' ', narc), uves[call[arcs]].lamp, FILNUM=17
              arc_fil = uves_getfil('arc_fil', setup, WCEN=twcen, $
                                     FRAME=uves[call[arcs[0]].frame, /name)
              arc_img = uves_getfil('arc_img', setup, WCEN=twcen, $
                                     FRAME=uves[call[arcs[0]].frame, /name)
          endif else printf, 17, $
            'uves_setup:  WARNING! No Arcs found '
          
          ;; STD
          std = where(uves[call].type EQ 'STD',nstd)
          if nstd NE 0 then begin
              printf, 17, 'uves_setup:   Found Standard stars-- ', $
                uves[call[std]].img_root, ' ', uves[call[std]].Obj 
              uves[call[std]].arc_fil = arc_fil
              uves[call[std]].flat_fil = qtz_fil
          endif else $
            printf, 17, 'uves_setup:  WARNING! No Standards found for slit '
          
          ;; Object
          obj = where(uves[call].type EQ 'OBJ',nobj)
          if nobj NE 0 then begin
              printf, 17, '  Objects: ============'
              uobj = uniq(uves[call[obj]].Obj, $
                          sort(uves[call[obj]].Obj))
              for tt=0L, n_elements(uobj)-1 do begin
                  ;; All matching objects
                  a = where( uves[call[obj]].Obj EQ $
                             (uves[call[obj]].Obj)[uobj[tt]], na )
                  ;; Objid
                  id = where(valuobj EQ $
                             (uves[call[obj]].Obj)[uobj[tt]], nid)
                  if nid NE 1 then stop
                  uves[call[obj[a]]].obj_id = id[0] + 1
                  ;; Set Arc Image
                  for mm=0L,na-1 do begin
                      if narc EQ 0 then continue
                      ;; Find arcfil
;                      arc_fil = uves_getfil(uves, call[obj[a[mm]]], $
;                                                call[arcs], ARC_img=arc_img)
                      arc_fil = hires_getarcfil(hires, call[obj[a[mm]]], $
                                                call[arcs], ARC_img=arc_img)
                      uves[call[obj[a[mm]]]].arc_fil = arc_fil
                      uves[call[obj[a[mm]]]].arc_img = arc_img
                      ;; Set Object file
                      uves[call[obj[a[mm]]]].obj_fil = $
                        uves_getfil('obj_fil', $
                                     OBJN=uves[call[obj[a[mm]]]].Obj, $
                                     WCEN=uves[call[obj[a[mm]]]].xdangl, $
                                     FRAME=uves[call[obj[a[mm]]]].frame, $
                                     /name)
                  endfor
                  ;; Set Flat Image
                  uves[call[obj[a]]].flat_fil = qtz_fil
                  ;; Set Final Image
                  lroot = strlen(uves[call[obj[a[0]]]].img_root)
                  uves[call[obj[a]]].img_final = $
                    uves_getfil('fin_fil', $
                                OBJN=strmid(uves[call[obj[a[0]]]].img_root, $
                                            0,lroot-5), /name)
                  ;; Print
                  for ss=0,na-1 do begin
                      printf, 17, $
                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
                        call[obj[a[ss]]], $
                        uves[call[obj[a[ss]]]].img_root, $
                        uves[call[obj[a[ss]]]].Obj, $
                        uves[call[obj[a[ss]]]].obj_id, $
                        round(uves[call[obj[a[ss]]]].exp), $
                        uves[call[obj[a[ss]]]].arc_fil, $
                        uves[call[obj[a[ss]]]].flat_fil
                      ;; summ
;                      if ss EQ 0 AND tt EQ 0 then $
;                        printf, 19, $
;                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
;                        call[obj[a[ss]]], $
;                        uves[call[obj[a[ss]]]].img_root, $
;                        uves[call[obj[a[ss]]]].Obj, $
;                        uves[call[obj[a[ss]]]].obj_id, $
;                        round(uves[call[obj[a[ss]]]].exp), $
;                        uves[call[obj[a[ss]]]].arc_fil, $
;                        uves[call[obj[a[ss]]]].flat_fil
                  endfor
                  ;; Increment
                  nobj_id = nobj_id + 1
              endfor
          endif
      endfor
      printf, 17, '---------------------------------------------------'
      printf, 17, '---------------------------------------------------'
      printf, 17, '---------------------------------------------------'
;      printf, 19, '---------------------------------------------------'
  endfor

  srt = sort(svset)
  for kk=0L,nset-1 do begin
      ;; Sort
      ii = srt[kk]
      ;; SETUP
      setup = svset[ii]

      ;; Find files
      prs = strsplit(setups[ii], '-', /extract)
      tslitwid = prs[0]
;      tblock = prs[1]
      sech  = prs[1]
      sxd   = prs[2]
      rbin = prs[3]
      cbin = prs[4]

;      printf, 19, setup, tslitwid, tblock, echang[long(sech)], xdang[long(sxd)], $
;        rbin, cbin, FORMAT='(i2,5x,a2,1x,a6,1x,f7.4,1x,f7.4,2x,i2,3x,i2)'
      printf, 19, setup, tslitwid, '     ', echang[long(sech)], xdang[long(sxd)], $
        rbin, cbin, FORMAT='(i2,5x,a3,1x,a5,1x,f7.4,1x,f7.1,2x,i2,3x,i2)'
      ;; Obj
      obj = where(uves.setup EQ setup and uves.type EQ 'OBJ',nobj)
      for qq=0L,nobj-1 do begin
          printf, 19, uves[obj[qq]].obj_id, uves[obj[qq]].Obj, $
            round(uves[obj[qq]].xdangl), uves[obj[qq]].frame, uves[obj[qq]].exp, $
            FORMAT='(10x,i2,1x,a12,1x,i3,1x,i4,1x,f7.1)'
      endfor
  endfor
  ;; ALL DONE
  print, 'uves_setup: All done!'
  close, 17
  close, 19

  return
end
