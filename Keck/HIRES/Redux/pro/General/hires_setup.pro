;+ 
; NAME:
; hires_setup   
;     Version 1.1
;
; PURPOSE:
;    Examines the HIRES structure, separates files according to the
;    echelle+xdangl, blocking filter and binning.  The routine then
;    creates summary files that describe the calibration files that
;    exist for each setup.
;
; CALLING SEQUENCE:
;  hires_setup, hires, OUTFIL=
;
; INPUTS:
;   hires --  Kast IDL structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  OUTFIL=  -- Name of ASCII file (default: hires_summ.txt')
;  /KEEPSET -- Dont modify the setup numbers (useful to reset the
;              arcs)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_setup, hires
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro hires_setup_val, hires, setups, ECHANG=echang, XDANG=xdang, ETOLER=etoler,$
                     XTOLER=xtoler, SET_NOBLCK=set_noblck

  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.00101

  nhir = n_elements(hires)
  ;; ECHANGL
  echang = uniq(hires.echangl, sort(hires.echangl))
  nech = n_elements(echang)
  msk = bytarr(nech)
  msk[*] = 1B
  lidx = lindgen(nech)
  for ll=0L,nech-1 do begin     ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(hires[echang[ll]].echangl-$
                   hires[echang].echangl) LT ETOLER AND $
                   lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  echang = hires[echang[where(msk)]].echangl
  nech = n_elements(echang)

  ;; XDANGL
  xdang = uniq(hires.xdangl, sort(hires.xdangl))
  nxd = n_elements(xdang)
  msk = bytarr(nxd)
  msk[*] = 1B
  lidx = lindgen(nxd)
  for ll=0L,nxd-1 do begin      ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(hires[xdang[ll]].xdangl-$
                       hires[xdang].xdangl) LT XTOLER $
                   AND lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  xdang = hires[xdang[where(msk)]].xdangl
  nxd = n_elements(xdang)

  ;; Setup
  setups = strarr(nhir)
  set_noblck = strarr(nhir)
  for qq=0L,nhir-1 do begin
      ;; Angles
      mtch = min(abs(hires[qq].echangl-echang),imn)
      echs = strtrim(imn,2)
      mtch = min(abs(hires[qq].xdangl-xdang),imn)
      xds = strtrim(imn,2)
      setups[qq] = hires[qq].decker+'-'+hires[qq].block+'-'+echs+'-'+xds+'-'+$
        strtrim(hires[qq].rowbin,2)+'-'+strtrim(hires[qq].colbin,2)
      ;; Ignoring blocking filter
      set_noblck[qq] = hires[qq].decker+'-'+echs+'-'+xds+'-'+$
        strtrim(hires[qq].rowbin,2)+'-'+strtrim(hires[qq].colbin,2)
  endfor
  setups = setups[uniq(setups, sort(setups))]
  set_noblck = set_noblck[uniq(set_noblck, sort(set_noblck))]

  return
end
      

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_setup, hires, OUTFIL=outfil, XTOLER=xtoler, ETOLER=etoler, $
                 KEEPSET=keepset

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_setup, hires, OUTFIL=, ETOLER=, XTOLER=, /KEEPSET [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.00101
  if not keyword_set( OUTFIL ) then outfil = 'hires_summ.txt'
  if not keyword_set( SUMMFIL ) then summfil = 'setup_summ.txt'

  qtz_fil = ' '
  arc_fil = ' '
  
  ;; Grab all  files
  gd = where(hires.decker NE 'D5' AND hires.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

  ;; Set D5
  alld5 = where(hires.decker EQ 'D5', nd5)
  if nd5 NE 0 then hires[alld5].setup = 99

  ;; Open File
  print, 'hires_setup: Check the files ', outfil+' '+summfil, ' for info'
  close, 17
  openw, 17, outfil
  close, 19
  openw, 19, summfil


  setup = 0L
  nobj_id = 0L

  ;; Loop on Decker
  deck = uniq(hires[gd].decker, sort(hires[gd].decker))
  ndeck = n_elements(deck)
  bd = where(hires[gd[deck]].decker EQ 'D5',nbd)
  if nbd NE 0 then stop

  hires_setup_val, hires[gd], setups, ECHANG=echang, XDANG=xdang, $
                   ETOLER=etoler, XTOLER=xtoler
  nset = n_elements(setups)

  if keyword_set(KEEPSET) then begin
      numset = hires[uniq(hires[gd].setup, sort(hires[gd].setup))].setup
      nset = n_elements(numset)
  endif

  ;; Summary
  printf, 19, 'Setup Deck  BLKF  ECHANG   XDANG  RBIN CBIN'

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
      tdeck = prs[0]
      tblock = prs[1]
      sech  = prs[2]
      sxd   = prs[3]
      rbin = prs[4]
      cbin = prs[5]

      mtch = where( hires[gd].decker EQ tdeck AND $
                    hires[gd].block EQ tblock AND $
                    abs(hires[gd].echangl-echang[long(sech)]) LT etoler AND $
                    abs(hires[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                    hires[gd].rowbin EQ long(rbin) AND $
                    hires[gd].colbin EQ long(cbin),nm)
      if nm EQ 0 then stop
      xall = gd[mtch]

      ;; Add in Arcs with different blocker on red side only (2nd
                       ;;  order light)
      ex_arcs = where( hires[gd].decker EQ tdeck AND $
                       hires[gd].block NE tblock AND $
                       abs(hires[gd].echangl-echang[long(sech)]) LT etoler AND $
                       abs(hires[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                       hires[gd].rowbin EQ long(rbin) AND $
                       hires[gd].colbin EQ long(cbin) AND $
                       hires[gd].cross EQ 'RED' AND $
                       hires[gd].type EQ 'ARC', nma)
      if nma NE 0 then begin
          print, 'hires_setup: Adding in arcs with different blocking'
          xall = [xall, gd[ex_arcs]]
      endif

      ;; Set Obj_id by Object name
      allobj = where(hires[xall].type EQ 'OBJ',naobj)
      if naobj NE 0 then begin
          uaobj = uniq(hires[xall[allobj]].Obj, $
                       sort(hires[xall[allobj]].Obj))
          valuobj = hires[xall[allobj[uaobj]]].Obj
          setup = objset
          objset = objset + 1
      endif else begin
          setup = calibset
          calibset = calibset + 1
      endelse
      if hires[xall[0]].setup NE setup AND hires[xall[0]].setup GT 0 then begin
          print, 'hires_setup:  About to change setup!!!  Think carefully...'
          stop
      endif

      hires[xall].setup = setup
      svset[kk] = setup
      
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, 'Setup: ', setup
      printf, 17, 'hires_setup: Decker = ', tdeck
      printf, 17, 'hires_setup: Blocking = ', tblock
      printf, 17, 'hires_setup: ECH Angle = ', echang[long(sech)]
      printf, 17, 'hires_setup: XDAngle = ', xdang[long(sxd)]
      printf, 17, 'hires_setup: Rowbin = ', rbin
      printf, 17, 'hires_setup: Colbin = ', cbin

;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, 'Setup: ', setup
;      printf, 19, 'hires_setup: Decker = ', tdeck
;      printf, 19, 'hires_setup: Blocking = ', tblock
;      printf, 19, 'hires_setup: ECH Angle = ', echang[long(sech)]
;      printf, 19, 'hires_setup: XDAngle = ', xdang[long(sxd)]
;      printf, 19, 'hires_setup: Rowbin = ', rbin
;      printf, 19, 'hires_setup: Colbin = ', cbin

      
      ;; Loop on Chip
      chip = uniq(hires[xall].chip, sort(hires[xall].chip))
      nchip = n_elements(chip)
      for pp=0L,nchip-1 do begin
          tchip = hires[xall[chip[pp]]].chip
          call = where(hires[xall].chip EQ tchip)
          call = xall[call]
          
          printf, 17, '-----'
          printf, 17, 'hires_setup: Chip (BGR) = ', tchip
          
          ;; QTZ
          qtz = where(hires[call].type EQ 'TFLT',nqtz)
          if nqtz NE 0 then begin
              printf, 17, 'hires_setup:   Found Trace Flats -- ', $
                hires[call[qtz]].img_root
              qtz_fil = hires_getfil('qtz_fil', $
                                     setup, CHIP=tchip, /name)
          endif else $
            printf, 17, $
            'hires_setup:  WARNING! No Trace Flats found..' 
          
          ;; Arcs
          arcs = where(hires[call].type EQ 'ARC',narc)
          if narc NE 0 then begin
              printf, 17, 'hires_setup:   Found Arc Frames-- ' + $
                '(file, arclamp)'
              writecol, 'dum', hires[call[arcs]].img_root, $
                replicate(' ', narc), hires[call[arcs]].lamp, FILNUM=17
          endif else printf, 17, $
            'hires_setup:  WARNING! No Arcs found '
          
          ;; STD
          std = where(hires[call].type EQ 'STD',nstd)
          if nstd NE 0 then begin
              printf, 17, 'hires_setup:   Found Standard stars-- ', $
                hires[call[std]].img_root, ' ', hires[call[std]].Obj 
;              hires[call[std]].arc_fil = arc_fil
              hires[call[std]].flat_fil = qtz_fil
          endif else $
            printf, 17, 'hires_setup:  WARNING! No Standards found for slit '
          
          ;; Object
          obj = where(hires[call].type EQ 'OBJ',nobj)
          if nobj NE 0 then begin
              printf, 17, '  Objects: ============'
              uobj = uniq(hires[call[obj]].Obj, $
                          sort(hires[call[obj]].Obj))
              for tt=0L, n_elements(uobj)-1 do begin
                  ;; All matching objects
                  a = where( hires[call[obj]].Obj EQ $
                             (hires[call[obj]].Obj)[uobj[tt]], na )
                  ;; Objid
                  id = where(valuobj EQ $
                             (hires[call[obj]].Obj)[uobj[tt]], nid)
                  if nid NE 1 then stop
                  hires[call[obj[a]]].obj_id = id[0] + 1
                  ;; Set Arc Image
                  for mm=0L,na-1 do begin
                      if narc EQ 0 then continue
                      ;; Find arcfil
                      arc_fil = hires_getarcfil(hires, call[obj[a[mm]]], $
                                                call[arcs], ARC_img=arc_img)
                      hires[call[obj[a[mm]]]].arc_fil = arc_fil
                      hires[call[obj[a[mm]]]].arc_img = arc_img
                      ;; Set Object file
                      hires[call[obj[a[mm]]]].obj_fil = $
                        hires_getfil('obj_fil', $
                                     CHIP=hires[call[obj[a[mm]]]].chip, $
                                     FRAME=hires[call[obj[a[mm]]]].frame, $
                                     /name)
                  endfor
                  ;; Set Flat Image
                  hires[call[obj[a]]].flat_fil = qtz_fil
                  ;; Set Final Image
                  hires[call[obj[a]]].img_final = $
                    hires_getfil('fin_fil', $
                                 CHIP=hires[call[obj[a[0]]]].chip, $
                                 FRAME=hires[call[obj[a[0]]]].frame, $
                                 /name)
                  ;; Print
                  for ss=0,na-1 do begin
                      printf, 17, $
                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
                        call[obj[a[ss]]], $
                        hires[call[obj[a[ss]]]].img_root, $
                        hires[call[obj[a[ss]]]].Obj, $
                        hires[call[obj[a[ss]]]].obj_id, $
                        round(hires[call[obj[a[ss]]]].exp), $
                        hires[call[obj[a[ss]]]].arc_fil, $
                        hires[call[obj[a[ss]]]].flat_fil
                      ;; summ
;                      if ss EQ 0 AND tt EQ 0 then $
;                        printf, 19, $
;                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
;                        call[obj[a[ss]]], $
;                        hires[call[obj[a[ss]]]].img_root, $
;                        hires[call[obj[a[ss]]]].Obj, $
;                        hires[call[obj[a[ss]]]].obj_id, $
;                        round(hires[call[obj[a[ss]]]].exp), $
;                        hires[call[obj[a[ss]]]].arc_fil, $
;                        hires[call[obj[a[ss]]]].flat_fil
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
      tdeck = prs[0]
      tblock = prs[1]
      sech  = prs[2]
      sxd   = prs[3]
      rbin = prs[4]
      cbin = prs[5]

      printf, 19, setup, tdeck, tblock, echang[long(sech)], xdang[long(sxd)], $
        rbin, cbin, FORMAT='(i2,5x,a2,1x,a6,1x,f7.4,1x,f7.4,2x,i2,3x,i2)'
      ;; Obj
      obj = where(hires.setup EQ setup and hires.type EQ 'OBJ',nobj)
      for qq=0L,nobj-1 do begin
          printf, 19, hires[obj[qq]].obj_id, hires[obj[qq]].Obj, $
            hires[obj[qq]].chip, hires[obj[qq]].frame, hires[obj[qq]].exp, $
            FORMAT='(10x,i2,1x,a12,1x,i2,1x,i4,1x,f7.1)'
      endfor
  endfor
  ;; ALL DONE
  print, 'hires_setup: All done!'
  close, 17
  close, 19

  return
end
