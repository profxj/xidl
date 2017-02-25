;+ 
; NAME:
; apf_setup   
;     Version 1.1
;
; PURPOSE:
;    Examines the HIRES structure, separates files according to the
;    echelle+xdangl, blocking filter and binning.  The routine then
;    creates summary files that describe the calibration files that
;    exist for each setup.
;
; CALLING SEQUENCE:
;  apf_setup, apf, OUTFIL=
;
; INPUTS:
;   apf --  Kast IDL structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  OUTFIL=  -- Name of ASCII file (default: apf_summ.txt')
;  /KEEPSET -- Dont modify the setup numbers (useful to reset the
;              arcs)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   apf_setup, apf
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro apf_setup_val, apf, setups, ECHANG=echang, XDANG=xdang, ETOLER=etoler,$
                     XTOLER=xtoler, SET_NOBLCK=set_noblck

  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.00101

  nhir = n_elements(apf)
  ;; ECHANGL
  echang = uniq(apf.echangl, sort(apf.echangl))
  nech = n_elements(echang)
  msk = bytarr(nech)
  msk[*] = 1B
  lidx = lindgen(nech)
  for ll=0L,nech-1 do begin     ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(apf[echang[ll]].echangl-$
                   apf[echang].echangl) LT ETOLER AND $
                   lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  echang = apf[echang[where(msk)]].echangl
  nech = n_elements(echang)

  ;; XDANGL
  xdang = uniq(apf.xdangl, sort(apf.xdangl))
  nxd = n_elements(xdang)
  msk = bytarr(nxd)
  msk[*] = 1B
  lidx = lindgen(nxd)
  for ll=0L,nxd-1 do begin      ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(apf[xdang[ll]].xdangl-$
                       apf[xdang].xdangl) LT XTOLER $
                   AND lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  xdang = apf[xdang[where(msk)]].xdangl
  nxd = n_elements(xdang)

  ;; Setup
  setups = strarr(nhir)
  set_noblck = strarr(nhir)
  for qq=0L,nhir-1 do begin
      ;; Angles
      mtch = min(abs(apf[qq].echangl-echang),imn)
      echs = strtrim(imn,2)
      mtch = min(abs(apf[qq].xdangl-xdang),imn)
      xds = strtrim(imn,2)
      setups[qq] = apf[qq].decker+'-'+apf[qq].block+'-'+echs+'-'+xds+'-'+$
        strtrim(apf[qq].rowbin,2)+'-'+strtrim(apf[qq].colbin,2)
      ;; Ignoring blocking filter
      set_noblck[qq] = apf[qq].decker+'-'+echs+'-'+xds+'-'+$
        strtrim(apf[qq].rowbin,2)+'-'+strtrim(apf[qq].colbin,2)
  endfor
  setups = setups[uniq(setups, sort(setups))]
  set_noblck = set_noblck[uniq(set_noblck, sort(set_noblck))]

  return
end
      

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro apf_setup, apf, OUTFIL=outfil, XTOLER=xtoler, ETOLER=etoler, $
                 KEEPSET=keepset

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'apf_setup, apf, OUTFIL=, ETOLER=, XTOLER=, /KEEPSET [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.00101
  if not keyword_set( OUTFIL ) then outfil = 'apf_summ.txt'
  if not keyword_set( SUMMFIL ) then summfil = 'setup_summ.txt'

  qtz_fil = ' '
  arc_fil = ' '
  
  ;; Grab all  files
  gd = where(apf.decker NE 'D5' AND apf.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

  ;; Set D5
  alld5 = where(apf.decker EQ 'D5', nd5)
  if nd5 NE 0 then apf[alld5].setup = 99

  ;; Open File
  print, 'apf_setup: Check the files ', outfil+' '+summfil, ' for info'
  close, 17
  openw, 17, outfil
  close, 19
  openw, 19, summfil


  setup = 0L
  nobj_id = 0L

  ;; Loop on Decker
  deck = uniq(apf[gd].decker, sort(apf[gd].decker))
  ndeck = n_elements(deck)
  bd = where(apf[gd[deck]].decker EQ 'D5',nbd)
  if nbd NE 0 then stop

  apf_setup_val, apf[gd], setups, ECHANG=echang, XDANG=xdang, $
                   ETOLER=etoler, XTOLER=xtoler
  nset = n_elements(setups)

  if keyword_set(KEEPSET) then begin
      numset = apf[uniq(apf[gd].setup, sort(apf[gd].setup))].setup
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

      mtch = where( apf[gd].decker EQ tdeck AND $
                    apf[gd].block EQ tblock AND $
                    abs(apf[gd].echangl-echang[long(sech)]) LT etoler AND $
                    abs(apf[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                    apf[gd].rowbin EQ long(rbin) AND $
                    apf[gd].colbin EQ long(cbin),nm)
      if nm EQ 0 then stop
      xall = gd[mtch]

      ;; Add in Arcs with different blocker on red side only (2nd
                       ;;  order light)
      ex_arcs = where( apf[gd].decker EQ tdeck AND $
                       apf[gd].block NE tblock AND $
                       abs(apf[gd].echangl-echang[long(sech)]) LT etoler AND $
                       abs(apf[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                       apf[gd].rowbin EQ long(rbin) AND $
                       apf[gd].colbin EQ long(cbin) AND $
                       apf[gd].cross EQ 'RED' AND $
                       apf[gd].type EQ 'ARC', nma)
      if nma NE 0 then begin
          print, 'apf_setup: Adding in arcs with different blocking'
          xall = [xall, gd[ex_arcs]]
      endif

      ;; Set Obj_id by Object name
      allobj = where(apf[xall].type EQ 'OBJ',naobj)
      if naobj NE 0 then begin
          uaobj = uniq(apf[xall[allobj]].Obj, $
                       sort(apf[xall[allobj]].Obj))
          valuobj = apf[xall[allobj[uaobj]]].Obj
          setup = objset
          objset = objset + 1
      endif else begin
          setup = calibset
          calibset = calibset + 1
      endelse
      if apf[xall[0]].setup NE setup AND apf[xall[0]].setup GT 0 then begin
          print, 'apf_setup:  About to change setup!!!  Think carefully...'
          stop
      endif

      apf[xall].setup = setup
      svset[kk] = setup
      
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, 'Setup: ', setup
      printf, 17, 'apf_setup: Decker = ', tdeck
      printf, 17, 'apf_setup: Blocking = ', tblock
      printf, 17, 'apf_setup: ECH Angle = ', echang[long(sech)]
      printf, 17, 'apf_setup: XDAngle = ', xdang[long(sxd)]
      printf, 17, 'apf_setup: Rowbin = ', rbin
      printf, 17, 'apf_setup: Colbin = ', cbin

;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, 'Setup: ', setup
;      printf, 19, 'apf_setup: Decker = ', tdeck
;      printf, 19, 'apf_setup: Blocking = ', tblock
;      printf, 19, 'apf_setup: ECH Angle = ', echang[long(sech)]
;      printf, 19, 'apf_setup: XDAngle = ', xdang[long(sxd)]
;      printf, 19, 'apf_setup: Rowbin = ', rbin
;      printf, 19, 'apf_setup: Colbin = ', cbin

      
      ;; Loop on Chip
      chip = uniq(apf[xall].chip, sort(apf[xall].chip))
      nchip = n_elements(chip)
      for pp=0L,nchip-1 do begin
          tchip = apf[xall[chip[pp]]].chip
          call = where(apf[xall].chip EQ tchip)
          call = xall[call]
          
          printf, 17, '-----'
          printf, 17, 'apf_setup: Chip (BGR) = ', tchip
          
          ;; QTZ
          qtz = where(apf[call].type EQ 'TFLT',nqtz)
          if nqtz NE 0 then begin
              printf, 17, 'apf_setup:   Found Flats -- ', $
                apf[call[qtz]].img_root
              qtz_fil = apf_getfil('qtz_fil', $
                                     setup, CHIP=tchip, /name)
          endif else $
            printf, 17, $
            'apf_setup:  WARNING! No Trace Flats found..' 
          
          ;; Arcs
          arcs = where(apf[call].type EQ 'ARC',narc)
          if narc NE 0 then begin
              printf, 17, 'apf_setup:   Found Arc Frames-- ' + $
                '(file, arclamp)'
              writecol, 'dum', apf[call[arcs]].img_root, $
                replicate(' ', narc), apf[call[arcs]].lamp, FILNUM=17
          endif else printf, 17, $
            'apf_setup:  WARNING! No Arcs found '
          
          ;; STD
          std = where(apf[call].type EQ 'STD',nstd)
          if nstd NE 0 then begin
              printf, 17, 'apf_setup:   Found Standard stars-- ', $
                apf[call[std]].img_root, ' ', apf[call[std]].Obj 
;              apf[call[std]].arc_fil = arc_fil
              apf[call[std]].flat_fil = qtz_fil
          endif else $
            printf, 17, 'apf_setup:  WARNING! No Standards found for slit '
          
          ;; Object
          obj = where(apf[call].type EQ 'OBJ',nobj)
          if nobj NE 0 then begin
              printf, 17, '  Objects: ============'
              uobj = uniq(apf[call[obj]].Obj, $
                          sort(apf[call[obj]].Obj))
              for tt=0L, n_elements(uobj)-1 do begin
                  ;; All matching objects
                  a = where( apf[call[obj]].Obj EQ $
                             (apf[call[obj]].Obj)[uobj[tt]], na )
                  ;; Objid
                  id = where(valuobj EQ $
                             (apf[call[obj]].Obj)[uobj[tt]], nid)
                  if nid NE 1 then stop
                  apf[call[obj[a]]].obj_id = id[0] + 1
                  ;; Set Arc Image
                  for mm=0L,na-1 do begin
                      if narc EQ 0 then continue
                      ;; Find arcfil
                      arc_fil = apf_getarcfil(apf, call[obj[a[mm]]], $
                                                call[arcs], ARC_img=arc_img)
                      apf[call[obj[a[mm]]]].arc_fil = arc_fil
                      apf[call[obj[a[mm]]]].arc_img = arc_img
                      ;; Set Object file
                      apf[call[obj[a[mm]]]].obj_fil = $
                        apf_getfil('obj_fil', $
                                     CHIP=apf[call[obj[a[mm]]]].chip, $
                                     FRAME=apf[call[obj[a[mm]]]].frame, $
                                     /name)
                  endfor
                  ;; Set Flat Image
                  apf[call[obj[a]]].flat_fil = qtz_fil
                  ;; Set Final Image
                  apf[call[obj[a]]].img_final = $
                    apf_getfil('fin_fil', $
                                 CHIP=apf[call[obj[a[0]]]].chip, $
                                 FRAME=apf[call[obj[a[0]]]].frame, $
                                 /name)
                  ;; Print
                  for ss=0,na-1 do begin
                      printf, 17, $
                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
                        call[obj[a[ss]]], $
                        apf[call[obj[a[ss]]]].img_root, $
                        apf[call[obj[a[ss]]]].Obj, $
                        apf[call[obj[a[ss]]]].obj_id, $
                        round(apf[call[obj[a[ss]]]].exp), $
                        apf[call[obj[a[ss]]]].arc_fil, $
                        apf[call[obj[a[ss]]]].flat_fil
                      ;; summ
;                      if ss EQ 0 AND tt EQ 0 then $
;                        printf, 19, $
;                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
;                        call[obj[a[ss]]], $
;                        apf[call[obj[a[ss]]]].img_root, $
;                        apf[call[obj[a[ss]]]].Obj, $
;                        apf[call[obj[a[ss]]]].obj_id, $
;                        round(apf[call[obj[a[ss]]]].exp), $
;                        apf[call[obj[a[ss]]]].arc_fil, $
;                        apf[call[obj[a[ss]]]].flat_fil
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
        rbin, cbin, FORMAT='(i2,5x,a7,1x,a6,1x,f7.4,1x,f7.4,2x,i2,3x,i2)'
      ;; Obj
      obj = where(apf.setup EQ setup and apf.type EQ 'OBJ',nobj)
      for qq=0L,nobj-1 do begin
          printf, 19, apf[obj[qq]].obj_id, apf[obj[qq]].Obj, $
            apf[obj[qq]].chip, apf[obj[qq]].frame, apf[obj[qq]].exp, $
            FORMAT='(10x,i2,1x,a12,1x,i2,1x,i5,1x,f7.1)'
      endfor
  endfor
  ;; ALL DONE
  print, 'apf_setup: All done!'
  close, 17
  close, 19

  return
end
