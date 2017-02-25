



PRO mase_event, ev

COMMON share, magestrct, flatinter, singlesinter
 
COMPILE_OPT hidden

  widget_control, ev.id, get_uvalue=uvalue ;event handler

IF N_ELEMENTS(uvalue) EQ 1 THEN CASE uvalue OF
   'list':
   'tab':
   'quit':  WIDGET_CONTROL, ev.TOP, /DESTROY
   'about': about = DIALOG_MESSAGE(['MAge Spectral Extractor v0.2 - Written by John Bochanski', '', 'For more information:  jjb@mit.edu', 'http://web.mit.edu/jjb/www/MASE.html'], title='About')


   'raw_pick': BEGIN
       ;get IDs from the main widget
      WIDGET_CONTROL, ev.top, get_uvalue=s
      WIDGET_CONTROL, s.setupdir, get_value=inipath
      ;pick working directory
      path=DIALOG_PICKFILE(/Directory, /must_exist, title='Browse', path=inipath)
     
      ;figure out which setup dir, send the path back to the top
      WIDGET_CONTROL, s.setupdir, set_value=path
      WIDGET_CONTROL, s.wtext, set_value='Raw path has been updated to '+path, /append
      ;;cd, path
      ;;cd, '..'
      spawn, 'pwd', cwd
      WIDGET_CONTROL, s.wtext, set_value='Current working directory is'+cwd, /append
   END

 'blue_trace_pick': BEGIN
      ;pick working directory
      WIDGET_CONTROL, ev.top, get_uvalue=s
      WIDGET_CONTROL, s.setupdir, get_value=rawpath
      file=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist)
      ;figure out which setup dir, send the path back to the top
      WIDGET_CONTROL, s.skyfile, set_value=file
      WIDGET_CONTROL, s.wtext, set_value='Sky trace flat added as '+file, /append
   END
 
 'green_trace_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    file=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.greenfile, set_value=file
    WIDGET_CONTROL, s.wtext, set_value='Green trace flat added as '+file, /append
 END

 'red_trace_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    file=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.redfile, set_value=file
    WIDGET_CONTROL, s.wtext, set_value='Red trace flat added as '+file, /append
 END
 'go_trace': BEGIN
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    WIDGET_CONTROL, s.redfile, get_value=redfile
    WIDGET_CONTROL, s.greenfile, get_value=greenfile
    WIDGET_CONTROL, s.skyfile, get_value=bluefile
    ;do simple coadd here
    ;;a = xmrdfits(redfile[0], 0, hdr)
    ;;b = xmrdfits(greenfile[0], 0)
    ;;c = xmrdfits(bluefile[0], 0)
    ;;d = a+b+c
    ;;mwrfits, d, 'temp_trace.fits', hdr, /create
                                ;While I don't pass tset_slits
                                ;back right now, there are two files
                                ;created "Orders.fits" and
                                ;"Ostr_mage.fits".  The 1st extension
                                ;of Orders.fits is the same as
                                ;tset_slits
    tracefiles =  [redfile[0], greenfile[0], bluefile[0]]
    WIDGET_CONTROL, s.wtext, set_value='Launching ATV;  Press q to continue after inspecting trace', /append
    tset_slits=mage_findslits(tracefiles)
    WIDGET_CONTROL, s.wtext, set_value='Orders have been traced.', /append
 END

 'pix_arc_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    file=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.pixarcfile, set_value=file
    WIDGET_CONTROL, s.wtext, set_value='Pixel arc for flats added as '+file, /append
 END

;;  'go_piximage': BEGIN
;;     WIDGET_CONTROL, ev.top, get_uvalue=s
;;     WIDGET_CONTROL, s.setupdir, get_value=rawpath
;;     WIDGET_CONTROL, s.pixarcfile, get_value=pixarcfile
;;     WIDGET_CONTROL, s.piximgfile, get_value=piximgfile
;;     tset_slits=mrdfits('Orders.fits', 1)
;;     IF pixinter THEN WIDGET_CONTROL, s.wtext, set_value='Launching ATV;  Press q to continue through each order', /append
;;     mage_makepix, pixarcfile, piximgfile, tset_slits, chk=pixinter
;;     WIDGET_CONTROL, s.wtext, set_value='Pixel Image written to '+ piximgfile, /append
;;  END

;;  'pix_arc_inter_tog':BEGIN
;;     WIDGET_CONTROL, ev.top, get_uvalue=s
;;     pixinter = ev.select
    
;;  END
 


 'red_flat_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    files=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist, /multiple_files)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.redflatfiles, set_value=files
    WIDGET_CONTROL, s.wtext, set_value='Red flat files added as '+files, /append
 END

 'green_flat_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    files=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist, /multiple_files)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.greenflatfiles, set_value=files
    WIDGET_CONTROL, s.wtext, set_value='Green flat files added as '+files, /append
 END

 'blue_flat_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    files=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist, /multiple_files)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.blueflatfiles, set_value=files
    WIDGET_CONTROL, s.wtext, set_value='Blue flat files added as '+files, /append
 END

 'illum_flat_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    files=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist, /multiple_files)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.illumflatfiles, set_value=files
    WIDGET_CONTROL, s.wtext, set_value='Illum flat files added as '+files, /append
 END

 'flat_inter_tog':BEGIN
    WIDGET_CONTROL, ev.top, get_uvalue=s
    flatinter = ev.select
    
 END

 'single_tog':BEGIN
    WIDGET_CONTROL, ev.top, get_uvalue=s
    singlesinter = ev.select
 END


 'go_flat': BEGIN
    WIDGET_CONTROL, ev.top, get_uvalue=s
    
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    WIDGET_CONTROL, s.pixarcfile, get_value=pixarcfile
    WIDGET_CONTROL, s.redflatfiles, get_value=redflatfiles
    WIDGET_CONTROL, s.greenflatfiles, get_value=greenflatfiles
    WIDGET_CONTROL, s.blueflatfiles, get_value=blueflatfiles
    WIDGET_CONTROL, s.illumflatfiles, get_value=illumflatfiles
    ;;WIDGET_CONTROL, s.flatfile, get_value=flatfile
    WIDGET_CONTROL, s.wtext, set_value='Generating Flat Field...this could take a while', /append    
    WIDGET_CONTROL, s.wtext, set_value = 'Press q after xatv displays the flat', /append 

    
    mage_makeflat, red=redflatfiles, green=greenflatfiles, blue=blueflatfiles, illum=illumflatfiles, orders='Orders.fits', arcfile = pixarcfile[0], chk=flatinter ;;, output=flatfile[0]  ;automatically puts flatfiles in cwd/Flat/
;    pixflatfile='./Flat/'+flatfile[0]
;    illumflatfile='./Flat/Illum.fits'

    WIDGET_CONTROL, s.wtext, set_value='Flat field creation finished properly', /append
 END

  'make_struct':BEGIN
                                ;This generates the mage structure,
                                ;but it needs to know how many files
                                ;are in it.  So here we go...make the
                                ;structure, write it to a file, have
                                ;another program read the file, count
                                ;the number of stars, then update the
                                ;main state structure with the proper
                                ;number of stars

     
     WIDGET_CONTROL, ev.top, get_uvalue=s
     WIDGET_CONTROL, s.setupdir, get_value=rawpath
     WIDGET_CONTROL, s.strctname, get_value=strctname
     ;make structure
     mage_mkstrct, mage, rawpath=rawpath[0]
     ;write it to a file
     mwrfits, mage, strctname[0], /create 
     magestrct = mage
     targets = where(magestrct.obj_id GT 0 AND $
                      strcompress(magestrct.exptype,/rem) EQ 'SCIENCE' OR strcompress(magestrct.exptype,/rem) EQ 'BRIGHT')
     WIDGET_CONTROL, s.w_spec_target, set_value=magestrct(targets).object
     WIDGET_CONTROL, s.w_spec_target, set_uvalue=magestrct(targets).object
     WIDGET_CONTROL, s.wtext, set_value='Mage structure created and written to '+strctname[0], /append
  END


  'write_struct':BEGIN
    
     WIDGET_CONTROL, ev.top, get_uvalue=s
     ;WIDGET_CONTROL, s.setupdir, get_value=rawpath
     WIDGET_CONTROL, s.strctname, get_value=strctname
     mage = magestrct
     mwrfits, mage, strctname[0], /create
      targets = where(magestrct.obj_id GT 0 AND $
                      strcompress(magestrct.exptype,/rem) EQ 'SCIENCE' OR strcompress(magestrct.exptype,/rem) EQ 'BRIGHT')
     WIDGET_CONTROL, s.w_spec_target, set_value=magestrct(targets).object
     WIDGET_CONTROL, s.w_spec_target, set_uvalue=magestrct(targets).object
     WIDGET_CONTROL, s.wtext, set_value='Mage structure written to '+strctname[0], /append
  END

  'load_struct':BEGIN
     WIDGET_CONTROL, ev.top, get_uvalue=s
     WIDGET_CONTROL, s.strctname, get_value=strctname
     mage = xmrdfits(strctname[0], 1)
     magestrct = mage
     targets = where(magestrct.obj_id GT 0 AND $
                     strcompress(magestrct.exptype,/rem) EQ 'SCIENCE' OR strcompress(magestrct.exptype,/rem) EQ 'BRIGHT')
     ;; If there are no science targrets, just list everything
     IF targets[0] EQ -1 THEN targets=lindgen(n_elements(magestrct))
     WIDGET_CONTROL, s.w_spec_target, set_value=magestrct(targets).object
     WIDGET_CONTROL, s.w_spec_target, set_uvalue=magestrct(targets).object
     WIDGET_CONTROL, s.wtext, set_value='Mage structure loaded from '+strctname[0], /append
  END

 'edit_struct':BEGIN
    WIDGET_CONTROL, ev.top, get_uvalue=s 
    WIDGET_CONTROL, s.strctname, get_value=strctname
    mage = xmrdfits(strctname[0], 1)
    mage_editstrct, strctname[0]
    mage = xmrdfits(strctname[0], 1)
    magestrct = mage
    targets = where(magestrct.obj_id GT 0 AND $
                 strcompress(magestrct.exptype,/rem) EQ 'SCIENCE' OR strcompress(magestrct.exptype,/rem) EQ 'BRIGHT', ntarg)
    IF ntarg GT 0 THEN BEGIN
    WIDGET_CONTROL, s.w_spec_target, set_value=magestrct(targets).object
    WIDGET_CONTROL, s.w_spec_target, set_uvalue=magestrct(targets).object
    ENDIF
    WIDGET_CONTROL, s.wtext, set_value='Mage structure has been edited and written to '+strctname[0], /append

  END

 'std_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    file=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.stdfile, set_value=file
    WIDGET_CONTROL, s.wtext, set_value='Standard star file added as '+file, /append
 END

;;  'std_arc_pick': BEGIN
;;                                 ;pick working directory
;;     WIDGET_CONTROL, ev.top, get_uvalue=s
;;     WIDGET_CONTROL, s.setupdir, get_value=rawpath
;;     file=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist)
;;                                 ;figure out which setup dir, send the path back to the top
;;     WIDGET_CONTROL, s.stdarcfile, set_value=file
;;     WIDGET_CONTROL, s.wtext, set_value='Standard star arc file added as '+file, /append
;;  END
 
 'std_flux_pick': BEGIN
                                ;pick working directory
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    file=DIALOG_PICKFILE(filter='*.dat', /read, path=getenv('XIDL_DIR')+'/Spec/Flux/', /must_exist)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.stdflxfile, set_value=file
    WIDGET_CONTROL, s.wtext, set_value='Standard star flux table added as '+file, /append
 END

 'go_sensfunc': BEGIN
    mage = magestrct
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.stdfile, get_value=stdfile
    ;;WIDGET_CONTROL, s.stdarcfile, get_value=stdarcfile
    WIDGET_CONTROL, s.stdflxfile, get_value=stdflxfile
    WIDGET_CONTROL, s.stdoutfile, get_value=stdoutfile
    ;;WIDGET_CONTROL, s.piximgfile, get_value=piximgfile
    ;;WIDGET_CONTROL, s.flatfile, get_value=flatfile
    WIDGET_CONTROL, s.wtext, set_value='Computing sensitivity function.  Press H for help on the interactive fitting ', /append
    stdfile = strsplit(stdfile,'\.gz',/extra,/regex)
    istd = WHERE(mage.FITSFILE EQ (fileandpath(stdfile))[0])
    mage[istd].pixflatfile = 'Flat/Pixflat.fits'
    mage[istd].illumflatfile  = 'Flat/Illumflat_' + $
                                strcompress(mage[istd].slit, /rem) + '.fits'
    mage.orderfile      = 'OStr_mage.fits'
    mage.slitfile       = 'Orders.fits'
    mage_sensfunc,mage,stdfile[0],stdflxfile[0],stdoutfile[0] 
    ;;mage_sensfunc, mage, stdfile[0], arcfile=stdarcfile[0], fluxtable=stdflxfile[0], sensfile=stdoutfile[0], pixflatfile='./Flat/'+flatfile[0], pixfile=piximgfile[0]
END 

 'go_pipe': BEGIN
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.stdoutfile, get_value=stdoutfile
    WIDGET_CONTROL, s.stdfile, get_value=stdfile
    ;;WIDGET_CONTROL, s.piximgfile, get_value=piximgfile
    ;;WIDGET_CONTROL, s.flatfile, get_value=flatfile
    WIDGET_CONTROL, s.w_spec_target, get_uvalue=targetfiles
    print, 'targ', targetfiles ;this actually works...now I need to figure out which ones are highlighted...
    mage = magestrct
    WIDGET_CONTROL, s.wtext, set_value='Running pipeline.  Go get lunch, this will take a while...', /append

    ;; Assign the required files to standard locations
    mage.pixflatfile    = 'Flat/Pixflat.fits'
    mage.illumflatfile  = 'Flat/Illumflat_' + $
                          strcompress(mage.SLIT, /rem) + '.fits'
    mage.orderfile      = 'OStr_mage.fits'
    mage.slitfile       = 'Orders.fits'
    objfile             = 'Object/ObjStr' + $
                          strcompress(strsplit(fileandpath(stdfile[0]) $
                                               , 'mage', /extract,/regex), /rem)
    mage.stdfile        = (strsplit(objfile[0],'\.gz',/extra,/regex))[0]
    mage.sensfunc = stdoutfile[0]
    mage_pipe, mage, singles = singlesinter
 END

 'ql_file_pick':BEGIN
    WIDGET_CONTROL, ev.top, get_uvalue=s
    WIDGET_CONTROL, s.setupdir, get_value=rawpath
    file=DIALOG_PICKFILE(filter='*.fit*', /read, path=rawpath, /must_exist)
                                ;figure out which setup dir, send the path back to the top
    WIDGET_CONTROL, s.qlfile, set_value=file
    WIDGET_CONTROL, s.wtext, set_value='Quicklook target is '+file, /append
 END


 'go_ql':BEGIN
     WIDGET_CONTROL, ev.top, get_uvalue=s
     WIDGET_CONTROL, s.qlfile, get_value=qlfile
     mage_quicklook, qlfile
     WIDGET_CONTROL, s.wtext, set_value='Running Quicklook on '+qlfile, /append
     END

ENDCASE



END




;Main program
PRO mase

COMMON share, magestrct, flatinter, singlesinter

   flatinter=1  ;Default setup
   singlesinter=1 ;Default setup
   
  ;Get current working directory
  spawn, 'pwd', cwd
  spawn, 'whoami', name
  tmp = { magestrct }
  
  
                                                   
  main = widget_base(/Col, Title='MASE', MBAR=bar)      ; main base

  
  menu1 = WIDGET_BUTTON(bar, VALUE='File', /MENU)
  bhelp = WIDGET_BUTTON(bar, Value='Help', /help)
  babout = WIDGET_BUTTON(bhelp, VALUE='About', uvalue='about')
  
  wlabel = WIDGET_LABEL(main, value='Output from MASE')
  wtext = WIDGET_TEXT(main, XSIZE=60, ysize=20, /frame, value='Welcome Back to MASE, ' + name, /wrap, /scroll)  
  wtab = widget_tab(main, Location=0, uvalue='tab')


  wt1 = Widget_base(wtab, Title='Setup', /row)
  setup_label = WIDGET_LABEL(wt1, Value='Raw Directory')
  setupdir  = WIDGET_TEXT(wt1, /editable, value=cwd)
  setupbt =  WIDGET_BUTTON(wt1, Value='Browse', uvalue='raw_pick')

  ;Here the user selects one (or more) files, then it calls mage_traceorders
  wt2 = WIDGET_BASE(wtab, Title='Trace', /col, /base_align_right)
  wsky = WIDGET_BASE(wt2, /row)
  sky_label = WIDGET_LABEL(wsky, Value='Sky Flat')
  skyfile  = WIDGET_TEXT(wsky, /editable, value='', xsize=40)
  skybt =  WIDGET_BUTTON(wsky, Value='Browse', uvalue='blue_trace_pick')

  wgreen = WIDGET_BASE(wt2, /row)
  green_label = WIDGET_LABEL(wgreen, Value='Xe Trace Flat')
  greenfile  = WIDGET_TEXT(wgreen, /editable, value='', xsize=40)
  greenbt =  WIDGET_BUTTON(wgreen, Value='Browse', uvalue='green_trace_pick')

  wred = WIDGET_BASE(wt2, /row)
  red_label = WIDGET_LABEL(wred, Value='Red Flat')
  redfile  = WIDGET_TEXT(wred, /editable, value='', xsize=40)
  redbt =  WIDGET_BUTTON(wred, Value='Browse', uvalue='red_trace_pick')
  
  go_trace_bt = WIDGET_BUTTON(wt2, value='Trace Orders', uvalue='go_trace')

  ;Next tab does the pixel image
  ;; wt3 = WIDGET_BASE(wtab, Title='Pixel Image', /col,  /base_align_right)
;;   wpixarc = WIDGET_BASE(wt3, /row)
;;   pixarc_label = WIDGET_LABEL(wpixarc, Value='Pixel Arc')
;;   pixarcfile  = WIDGET_TEXT(wpixarc, /editable, value='', xsize=40)
;;   pixarcbt =  WIDGET_BUTTON(wpixarc, Value='Browse', uvalue='pix_arc_pick')
 
 
;;   wpiximg = WIDGET_BASE(wt3, /row)
;;   wpixarcinter = WIDGET_BASE(wt3, /nonexclusive)
;;   pixarcinterbt = WIDGET_BUTTON(wpixarcinter, Value='Interactive', uvalue='pix_arc_inter_tog', TOOLTIP='Toggle interactive viewing of pixel tilt in each order')
;;   Widget_Control, pixarcinterbt, Set_Button=1   ;Keep interactive on as a default


;;   piximg_label = WIDGET_LABEL(wpiximg, Value='Pixel Image')
;;   piximgfile  = WIDGET_TEXT(wpiximg, /editable, value='Piximg.fits')

;;   go_piximg_bt = WIDGET_BUTTON(wt3, value='Make Pixel Image', uvalue='go_piximage')


  ;Now we deal with flat files
  wt4 = WIDGET_BASE(wtab, Title='Flats', /col, /base_align_left)

  wleft = WIDGET_BASE(wt4, /row)
  wright = WIDGET_BASE(wt4, /row)
  wred = WIDGET_BASE(wleft, /row)
  red_label = WIDGET_LABEL(wred, Value='Red Flat Files')
  redflatfiles  = WIDGET_TEXT(wred, /editable, value='', xsize=40, ysize=4, /scroll)
  redbt =  WIDGET_BUTTON(wred, Value='Browse', uvalue='red_flat_pick')

  wgreen = WIDGET_BASE(wleft, /row)
  green_label = WIDGET_LABEL(wgreen, Value='Green Flat Files')
  greenflatfiles  = WIDGET_TEXT(wgreen, /editable, value='', xsize=40, ysize=4, /scroll)
  greenbt =  WIDGET_BUTTON(wgreen, Value='Browse', uvalue='green_flat_pick')

  wblue = WIDGET_BASE(wright, /row)
  blue_label = WIDGET_LABEL(wblue, Value='Blue Flat Files')
  blueflatfiles  = WIDGET_TEXT(wblue, /editable, value='', xsize=40, ysize=4, /scroll)
  bluebt =  WIDGET_BUTTON(wblue, Value='Browse', uvalue='blue_flat_pick')


  willum = WIDGET_BASE(wright, /row)
  illum_label = WIDGET_LABEL(willum, Value='Illum Flat Files')
  illumflatfiles  = WIDGET_TEXT(willum, /editable, value='', xsize=40, ysize=4, /scroll)
  illumbt =  WIDGET_BUTTON(willum, Value='Browse', uvalue='illum_flat_pick')

  ;wout = WIDGET_BASE(wt4, /row)
  ;flat_label = WIDGET_LABEL(wout, Value='Output Filename')
  ;flatfile  = WIDGET_TEXT(wout, /editable, value='Flat.fits', xsize=40)

  wbottom = WIDGET_BASE(wt4, /row)
 
  
  wflatpixarc_label = WIDGET_LABEL(wbottom, Value='Arc File')
  pixarcfile = WIDGET_TEXT(wbottom, /editable, value='', xsize=60)
  pixarcbt  = WIDGET_BUTTON(wbottom, Value='Browse', uvalue='pix_arc_pick') 
 
  wflatinter = WIDGET_BASE(wbottom, /nonexclusive)
  flatinterbt = WIDGET_BUTTON(wflatinter, Value='Interactive', uvalue='flat_inter_tog', TOOLTIP='Toggle interactive viewing of flats in each order')
  Widget_Control, flatinterbt, Set_Button=1 ;Keep interactive on as a default
  
  go_flat_bt = WIDGET_BUTTON(wbottom, value='Make Flat Field', uvalue='go_flat')
 
  ;Now make the MAGE structure
  wt5 = WIDGET_BASE(wtab, Title='Structure', /row, /base_align_right)
  strctlabel = WIDGET_LABEL(wt5, Value='Structure Filename')
  strctname = WIDGET_TEXT(wt5, Value='mage.fits', /editable)
  strctbt = WIDGET_BUTTON(wt5, Value='Generate Structure', uvalue='make_struct')
  sbut2 = WIDGET_BUTTON(wt5, Value='Load Structure', uvalue='load_struct')
  sbut3 = WIDGET_BUTTON(wt5, Value='Edit Structure',  uvalue='edit_struct')
  sbut4 = WIDGET_BUTTON(wt5, Value='Save Structure', uvalue='write_struct')

                                ;Sensitivity Function ... takes in
                                ;mage structure, fits file of standard
                                ;star, arcfile, a fluxtable and an
                                ;output

  wt6 = WIDGET_BASE(wtab, Title='Sens Func', /col, /base_align_right)
  wstd = WIDGET_BASE(wt6, /row)
  std_label = WIDGET_LABEL(wstd, Value='Standard Star')
  stdfile  = WIDGET_TEXT(wstd, /editable, value='', xsize=40)
  stdbt =  WIDGET_BUTTON(wstd, Value='Browse', uvalue='std_pick')

  ;;wstdarc = WIDGET_BASE(wt6, /row)
  ;;stdarc_label = WIDGET_LABEL(wstdarc, Value='Standard Star Arc')
  ;;stdarcfile  = WIDGET_TEXT(wstdarc, /editable, value='', xsize=40)
  ;;;stdarcbt =  WIDGET_BUTTON(wstdarc, Value='Browse', uvalue='std_arc_pick')

  wstdflx = WIDGET_BASE(wt6, /row)
  stdflx_label = WIDGET_LABEL(wstdflx, Value='Standard Star Flux Table')
  stdflxfile  = WIDGET_TEXT(wstdflx, /editable, value='', xsize=40)
  stdflxbt =  WIDGET_BUTTON(wstdflx, Value='Browse', uvalue='std_flux_pick')

  wstdout = WIDGET_BASE(wt6, /row)
  stdout_label = WIDGET_LABEL(wstdout, Value='Sensitivity Function')
  stdoutfile  = WIDGET_TEXT(wstdout, /editable, value=cwd+'/std.sav', xsize=40)
  
  go_sensfunc_bt = WIDGET_BUTTON(wt6, value='Make Sensitivity Function', uvalue='go_sensfunc')


  ;PIPELINE
  wt7 = WIDGET_BASE(wtab, Title='Reduce Data',  /col)
  wt7a = WIDGET_BASE(wt7, /row)
  go_pipe_bt = WIDGET_BUTTON(wt7a, Value='Run Pipeline', uvalue='go_pipe', xsize=350, ysize=150)
  wt7b = WIDGET_BASE(wt7a, /row)
  w_spec_target_label = WIDGET_LABEL(wt7b, Value ='Reduce ONLY these objects')
  w_spec_target =  WIDGET_LIST(wt7b, Value='', uvalue='select_target', /multiple, xsize=50, ysize=4)
  wt7c = WIDGET_BASE(wt7, /nonexclusive)
  singleinterbt = WIDGET_BUTTON(wt7c, Value='Singles', uvalue='single_tog', TOOLTIP='Extraction of each spectrum of an object (before coaddition)', /align_left)
  Widget_Control, singleinterbt, Set_Button=1 ;Keep interactive on as a default
  


  wt8 = WIDGET_BASE(wtab, Title='Quicklook', /col)
  wql = WIDGET_BASE(wt8, /row)
  ql_label = WIDGET_LABEL(wql, Value='Quicklook Target')
  qlfile  = WIDGET_TEXT(wql, /editable, value='', xsize=40)
  qlfilebt =  WIDGET_BUTTON(wql, Value='Browse', uvalue='ql_file_pick')

  qlbt = WIDGET_BUTTON(wt8, Value='Quicklook', uvalue='go_ql')
  
  ;button9 = WIDGET_BUTTON(menu1, VALUE='Open')
  bquit = WIDGET_BUTTON(menu1, VALUE='Quit', uvalue='quit')


  
  widget_control, main, /realize ; create the widgets


  state = {cwd:cwd, wtext:wtext, setupdir:setupdir, skyfile:skyfile, greenfile:greenfile, redfile:redfile, pixarcfile:pixarcfile, flatinterbt:flatinterbt, singleinterbt:singleinterbt, redflatfiles:redflatfiles, greenflatfiles:greenflatfiles, blueflatfiles:blueflatfiles, illumflatfiles:illumflatfiles, strctname:strctname,  stdfile:stdfile, stdflxfile:stdflxfile, stdoutfile:stdoutfile, w_spec_target:w_spec_target, qlfile:qlfile, bquit:bquit}



  widget_control, main, set_uvalue = state
  ;I think this makes the draw window the spot for plotting
  ;widget_control, draw, get_value=window
  ;wset, window
  
  xmanager, 'mase', main , /no_block      ; wait for events
END
