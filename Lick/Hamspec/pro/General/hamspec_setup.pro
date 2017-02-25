;+ 
; NAME:
; hamspec_setup   
;     Version 1.1
;
; PURPOSE:
;    Examines the HIRES structure, separates files according to the
;    echelle+xdangl, blocking filter and binning.  The routine then
;    creates summary files that describe the calibration files that
;    exist for each setup.
;
; CALLING SEQUENCE:
;  hamspec_setup, hamspec, OUTFIL=
;
; INPUTS:
;   hamspec --  Kast IDL structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  OUTFIL=  -- Name of ASCII file (default: hamspec_summ.txt')
;  /KEEPSET -- Dont modify the setup numbers (useful to reset the
;              arcs)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_setup, hamspec
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro hamspec_setup_val, hamspec, setups, ECHANG=echang, XDANG=xdang, ETOLER=etoler,$
                     XTOLER=xtoler, SET_NOBLCK=set_noblck

  if not keyword_set( ETOLER ) then etoler = 10
  if not keyword_set( XTOLER ) then xtoler = 10

  nhir = n_elements(hamspec)
  ;; ECHANGL
  echang = uniq(hamspec.echangl, sort(hamspec.echangl))
  nech = n_elements(echang)
  msk = bytarr(nech)
  msk[*] = 1B
  lidx = lindgen(nech)
  for ll=0L,nech-1 do begin     ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(hamspec[echang[ll]].echangl-$
                   hamspec[echang].echangl) LT ETOLER AND $
                   lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  echang = hamspec[echang[where(msk)]].echangl
  nech = n_elements(echang)

  ;; XDANGL
  xdang = uniq(hamspec.xdangl, sort(hamspec.xdangl))
  nxd = n_elements(xdang)
  msk = bytarr(nxd)
  msk[*] = 1B
  lidx = lindgen(nxd)
  for ll=0L,nxd-1 do begin      ; Cut on angles 
      if msk[ll] EQ 0B then continue
      mtch = where(abs(hamspec[xdang[ll]].xdangl-$
                       hamspec[xdang].xdangl) LT XTOLER $
                   AND lidx NE ll, nmtch)
      if nmtch NE 0 then msk[mtch] = 0B
  endfor
  xdang = hamspec[xdang[where(msk)]].xdangl
  nxd = n_elements(xdang)

  ;; Setup
  setups = strarr(nhir)
  set_noblck = strarr(nhir)
  for qq=0L,nhir-1 do begin
      ;; Angles
      mtch = min(abs(hamspec[qq].echangl-echang),imn)
      echs = strtrim(imn,2)
      mtch = min(abs(hamspec[qq].xdangl-xdang),imn)
      xds = strtrim(imn,2)
;      setups[qq] = hamspec[qq].decker+'-'+hamspec[qq].block+'-'+echs+'-'+xds+'-'+$
;        strtrim(hamspec[qq].rowbin,2)+'-'+strtrim(hamspec[qq].colbin,2)
      ;; Ignoring blocking filter
      setups[qq] = hamspec[qq].decker+'-'+echs+'-'+xds+'-'+$
        strtrim(hamspec[qq].rowbin,2)+'-'+strtrim(hamspec[qq].colbin,2)
  endfor
  setups = setups[uniq(setups, sort(setups))]
  set_noblck = set_noblck[uniq(set_noblck, sort(set_noblck))]

  return
end
      

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hamspec_setup, hamspec, OUTFIL=outfil, XTOLER=xtoler, ETOLER=etoler, $
                 KEEPSET=keepset

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hamspec_setup, hamspec, OUTFIL=, ETOLER=, XTOLER=, /KEEPSET [v1.1]'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.00101
  if not keyword_set( OUTFIL ) then outfil = 'hamspec_summ.txt'
  if not keyword_set( SUMMFIL ) then summfil = 'setup_summ.txt'

  qtz_fil = ' '
  arc_fil = ' '
  
  ;; Grab all  files
  gd = where(hamspec.decker NE 'D5' AND hamspec.flg_anly NE 0, ngd)
  if ngd EQ 0 then return

  ;; Set D5
  alld5 = where(hamspec.decker EQ 'D5', nd5)
  if nd5 NE 0 then hamspec[alld5].setup = 99

  ;; Open File
  print, 'hamspec_setup: Check the files ', outfil+' '+summfil, ' for info'
  close, 17
  openw, 17, outfil
  close, 19
  openw, 19, summfil


  setup = 0L
  nobj_id = 0L

  ;; Loop on Decker
  deck = uniq(hamspec[gd].decker, sort(hamspec[gd].decker))
  ndeck = n_elements(deck)
  bd = where(hamspec[gd[deck]].decker EQ 'D5',nbd)
  if nbd NE 0 then stop

  hamspec_setup_val, hamspec[gd], setups, ECHANG=echang, XDANG=xdang, $
                   ETOLER=etoler, XTOLER=xtoler, SET_NOBLCK=set_noblck
  nset = n_elements(setups)

  if keyword_set(KEEPSET) then begin
      numset = hamspec[uniq(hamspec[gd].setup, sort(hamspec[gd].setup))].setup
      nset = n_elements(numset)
  endif

  ;; Summary
;  printf, 19, 'Setup Deck  BLKF  ECHANG   XDANG  RBIN CBIN'
  printf, 19, 'Setup Deck   ECHANG   XDANG  RBIN CBIN'

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
;      tblock = prs[1]
      sech  = prs[1]
      sxd   = prs[2]
      rbin = prs[3]
      cbin = prs[4]


      decker = strsplit(strtrim(tdeck, 2),':',/extract)
      twidth = long(decker[0])
      tlength = float(decker[1])

      mtch = where( hamspec[gd].decker EQ tdeck AND $
;                    hamspec[gd].block EQ tblock AND $
                    abs(hamspec[gd].echangl-echang[long(sech)]) LT etoler AND $
                    abs(hamspec[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                    hamspec[gd].rowbin EQ long(rbin) AND $
                    hamspec[gd].type NE 'PFLT' and $
                    hamspec[gd].type NE 'MFLT' and $
                    hamspec[gd].colbin EQ long(cbin),nm)
      flg = 0
      if nm EQ 0 then flg= 1
      if flg then continue
         
      xall = gd[mtch]


      ;; Add in Flats with longer slit
      ex_flats = where( $;hamspec[gd].block NE tblock AND $
                       abs(hamspec[gd].echangl-echang[long(sech)]) LT etoler AND $
                       abs(hamspec[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                       hamspec[gd].length GT tlength AND $
                       hamspec[gd].rowbin EQ long(rbin) AND $
                       hamspec[gd].colbin EQ long(cbin) AND $
                       (hamspec[gd].type EQ 'TFLT' or hamspec[gd].type EQ 'PFLT'), nmf)
      if nmf NE 0 then begin
          print, 'hamspec_setup: Adding in flats with longer slits'
          print, 'hamspec_setup: And setting them to PFLT'
          xall = [xall, gd[ex_flats]]
          hamspec[gd[ex_flats]].type = 'PFLT'
      endif

      m_flats = where( $;hamspec[gd].block NE tblock AND $
                       abs(hamspec[gd].echangl-echang[long(sech)]) LT etoler AND $
                       abs(hamspec[gd].xdangl-xdang[long(sxd)]) LT xtoler AND $
                       hamspec[gd].length GT tlength AND $
                       hamspec[gd].rowbin EQ long(rbin) AND $
                       hamspec[gd].colbin EQ long(cbin) AND $
                       hamspec[gd].dewarfoc GT 1780 AND $
                       hamspec[gd].decker EQ '800:5.0' AND $ 
                       (hamspec[gd].type EQ 'TFLT' or hamspec[gd].type EQ 'PFLT'), nmf)
      if nmf NE 0 then begin
          print, 'hamspec_setup: Adding in flats with longer slits'
          print, 'hamspec_setup: And setting them to MFLT'
          xall = [xall, gd[m_flats]]
          hamspec[gd[m_flats]].type = 'MFLT'
      endif

      ;; Set Obj_id by Object name
      allobj = where(hamspec[xall].type EQ 'OBJ',naobj)
      if naobj NE 0 then begin
          uaobj = uniq(hamspec[xall[allobj]].Obj, $
                       sort(hamspec[xall[allobj]].Obj))
          valuobj = hamspec[xall[allobj[uaobj]]].Obj
          setup = objset
          objset = objset + 1
      endif else begin
          setup = calibset
          calibset = calibset + 1
      endelse
      if hamspec[xall[0]].setup NE setup AND hamspec[xall[0]].setup GT 0 then begin
          print, 'hamspec_setup:  About to change setup!!!  Think carefully...'
          stop
      endif

      hamspec[xall].setup = setup
      svset[kk] = setup
      
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, '------------------------------------------'
      printf, 17, 'Setup: ', setup
      printf, 17, 'hamspec_setup: Decker = ', tdeck
;      printf, 17, 'hamspec_setup: Blocking = ', tblock
      printf, 17, 'hamspec_setup: ECH Angle = ', echang[long(sech)]
      printf, 17, 'hamspec_setup: XDAngle = ', xdang[long(sxd)]
      printf, 17, 'hamspec_setup: Rowbin = ', rbin
      printf, 17, 'hamspec_setup: Colbin = ', cbin

;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, '------------------------------------------'
;      printf, 19, 'Setup: ', setup
;      printf, 19, 'hamspec_setup: Decker = ', tdeck
;      printf, 19, 'hamspec_setup: Blocking = ', tblock
;      printf, 19, 'hamspec_setup: ECH Angle = ', echang[long(sech)]
;      printf, 19, 'hamspec_setup: XDAngle = ', xdang[long(sxd)]
;      printf, 19, 'hamspec_setup: Rowbin = ', rbin
;      printf, 19, 'hamspec_setup: Colbin = ', cbin

      
      ;; Loop on Chip
;      chip = uniq(hamspec[xall].chip, sort(hamspec[xall].chip))
;      nchip = n_elements(chip)
      nchip = 1L
      for pp=0L,nchip-1 do begin
;          tchip = hamspec[xall[chip[pp]]].chip
;          call = where(hamspec[xall].chip EQ tchip)
          call = xall
          
          printf, 17, '-----'
;          printf, 17, 'hamspec_setup: Chip (BGR) = ', tchip
          
          ;; Trace flats
          qtz = where(hamspec[call].type EQ 'TFLT',nqtz)
          if nqtz NE 0 then begin
              printf, 17, 'hamspec_setup:   Found Trace Flats -- ', $
                hamspec[call[qtz]].img_root
              qtz_fil = hamspec_getfil('qtz_fil', setup, /name)
          endif else $
            printf, 17, $
            'hamspec_setup:  WARNING! No Trace Flats found..' 

          ;; Pixel flats
          pixf = where(hamspec[call].type EQ 'PFLT',npixf)
          if npixf NE 0 then begin
              printf, 17, 'hamspec_setup:   Found Pixel Flats -- ', $
                hamspec[call[pixf]].img_root
              pix_fil = hamspec_getfil('pixflt_fil', setup, /name)
          endif else $
            printf, 17, $
            'hamspec_setup:  WARNING! No Pixel Flats found..' 
          
          ;; Arcs
          arcs = where(hamspec[call].type EQ 'ARC',narc)
          if narc NE 0 then begin
              printf, 17, 'hamspec_setup:   Found Arc Frames-- ' + $
                '(file, arclamp)'
              writecol, 'dum', hamspec[call[arcs]].img_root, $
                replicate(' ', narc), hamspec[call[arcs]].lamp, FILNUM=17
          endif else printf, 17, $
            'hamspec_setup:  WARNING! No Arcs found '
          
          ;; STD
          std = where(hamspec[call].type EQ 'STD',nstd)
          if nstd NE 0 then begin
              printf, 17, 'hamspec_setup:   Found Standard stars-- ', $
                hamspec[call[std]].img_root, ' ', hamspec[call[std]].Obj 
;              hamspec[call[std]].arc_fil = arc_fil
              hamspec[call[std]].flat_fil = qtz_fil
          endif else $
            printf, 17, 'hamspec_setup:  WARNING! No Standards found for slit '
          
          ;; Object
          obj = where(hamspec[call].type EQ 'OBJ',nobj)
          if nobj NE 0 then begin
              printf, 17, '  Objects: ============'
              uobj = uniq(hamspec[call[obj]].Obj, $
                          sort(hamspec[call[obj]].Obj))
              for tt=0L, n_elements(uobj)-1 do begin
                  ;; All matching objects
                  a = where( hamspec[call[obj]].Obj EQ $
                             (hamspec[call[obj]].Obj)[uobj[tt]], na )
                  ;; Objid
                  id = where(valuobj EQ $
                             (hamspec[call[obj]].Obj)[uobj[tt]], nid)
                  if nid NE 1 then stop
                  hamspec[call[obj[a]]].obj_id = id[0] + 1
                  ;; Set Arc Image
                  for mm=0L,na-1 do begin
                      if narc EQ 0 then continue
                      ;; Find arcfil
                      arc_fil = hamspec_getarcfil(hamspec, call[obj[a[mm]]], $
                                                call[arcs], ARC_img=arc_img)
                      hamspec[call[obj[a[mm]]]].arc_fil = arc_fil
                      hamspec[call[obj[a[mm]]]].arc_img = arc_img
                      ;; Set Object file
                      hamspec[call[obj[a[mm]]]].obj_fil = $
                        hamspec_getfil('obj_fil', $
                                     FRAME=hamspec[call[obj[a[mm]]]].frame, $
                                     /name)
                  endfor
                  ;; Set Flat Image
                  hamspec[call[obj[a]]].flat_fil = qtz_fil
                  if npixf NE 0 then begin
                     hamspec[call[obj[a]]].pflat_fil = pix_fil
                  endif
                  ;; Set Final Image
                  hamspec[call[obj[a]]].img_final = $
                    hamspec_getfil('fin_fil', $
                                 FRAME=hamspec[call[obj[a[0]]]].frame, $
                                 /name)
                  ;; Print
                  for ss=0,na-1 do begin
                      printf, 17, $
                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
                        call[obj[a[ss]]], $
                        hamspec[call[obj[a[ss]]]].img_root, $
                        hamspec[call[obj[a[ss]]]].Obj, $
                        hamspec[call[obj[a[ss]]]].obj_id, $
                        round(hamspec[call[obj[a[ss]]]].exp), $
                        hamspec[call[obj[a[ss]]]].arc_fil, $
                        hamspec[call[obj[a[ss]]]].flat_fil
                      ;; summ
;                      if ss EQ 0 AND tt EQ 0 then $
;                        printf, 19, $
;                        FORMAT='(i4,1x,a9,1x,a12,1x,i3,i5,1x,a15,1x,2a17)', $
;                        call[obj[a[ss]]], $
;                        hamspec[call[obj[a[ss]]]].img_root, $
;                        hamspec[call[obj[a[ss]]]].Obj, $
;                        hamspec[call[obj[a[ss]]]].obj_id, $
;                        round(hamspec[call[obj[a[ss]]]].exp), $
;                        hamspec[call[obj[a[ss]]]].arc_fil, $
;                        hamspec[call[obj[a[ss]]]].flat_fil
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
;      tblock = prs[1]
      sech  = prs[1]
      sxd   = prs[2]
      rbin = prs[3]
      cbin = prs[4]

      printf, 19, setup, tdeck, echang[long(sech)], xdang[long(sxd)], $
        rbin, cbin, FORMAT='(i2,2x,a8,1x,i5,3x,i5,2x,i2,3x,i2)'
      ;; Obj
      obj = where(hamspec.setup EQ setup and hamspec.type EQ 'OBJ',nobj)
      for qq=0L,nobj-1 do begin
          printf, 19, hamspec[obj[qq]].obj_id, hamspec[obj[qq]].Obj, $
            hamspec[obj[qq]].frame, hamspec[obj[qq]].exp, $
            FORMAT='(10x,i2,1x,a12,1x,i4,1x,f7.1)'
      endfor
  endfor
  ;; ALL DONE
  print, 'hamspec_setup: All done!'
  close, 17
  close, 19

  return
end
