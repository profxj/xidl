;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; mase_non_gui.pro
; Author: Alex Hedglen                      Date: 8 Feb 2016
; Co-Author: Kathy Cooksey
; Project: Non-GUI, robust reduction of MagE data with
;          Kathy Cooksey
; 
; Description: This program reads a specifically formatted 
;              MagE configuration file of information regarding
;              a night's customization, and uses it to reduce
;              data for one night of observation. See mage_example.config
;
; Main call:
;   mase_non_gui, config_fil, [/no_edit, /no_trace,
;               /no_flat,/chk_flat,/reduce_std,/no_sens,
;               /no_pipe,/no_singles,/chk_trc]
;
;
; Inputs:
;   config_fil -- pipe (|) delimited file of
;                 param/tag | exp. #s | description
;                 (see mage_example.config)   
;
;   Default Parameters  -- dir           
;                          illum     
;                          blue        
;                          green       
;                          red           
;                          trace_arc  
;                          std1         
;                          fluxtab1 
;   
;   Optional Parameters -- Anything you want, but for editing the mage
;                          structure, the exptype parameters in the
;                          config_fil must be:
;                        
;                          trash
;                          science
;                          arc
;                          xe_flash
;                          domeflt
;                          bright     
;
;   For changing object names, the parameters in the configuration
;   file should be nc_* where the * is the name you want the
;   object's name to be changed to
;
;   If there is more than one standard star, you should add more
;   parameters named std* where * is the number 2,3,etc up to however
;   many standard stars you have. Then, the corresponding flux table
;   parameter should be added (fluxtab2,fluxtab3,etc.)
;
; Optional Inputs:
;   /debug     -- if set, print or stop at suitable places.
;
;   no_edit    -- option to edit the structure. If you know the
;                 structure has already been edited, you may set
;                 no_edit = 1 or /no_edit. You will want to do this
;                 when you are turning off other functions, so that
;                 all the variables you need are defined.
; 
;   no_trace   -- option to turn off trace orders function. Default is
;                 on. If you want off, you must define
;                 no_trace =  1 or /no_trace in the main call with
;                 mage_reduce_all, config_fil.
;
;   no_flat    -- option to turn off flats. Same deal as
;                 no_trace. If you want off set it as an optional
;                 call (no_flat = 1 or /no_flat).
;   
;   chk_trc    -- option to turn ON trace orders pop up
;                 ('Orders.fits'). Set chk_trc = 1 or /chk_trc.
;
;   chk_flat   -- option to have flats images pop up when
;                 finished (Flat/Illumflat_0.70.fits and
;                 Flat/Pixflat.fits). This option may be useful to
;                 check if everything is working properly. Default is
;                 off. To turn on, define in main call (chk_flat = 1 or
;                 /chk_flat). 
;
;   no_sens    -- to turn off the sensitivity function
;                 (mage_run_sensfunc), set no_sens = 1 or /no_sens. 
;
;   reduce_std -- if you turn off all functions but mage_run_sensfunc,
;                 you may want to reduce just the spectrophotometric
;                 standard star(s) and exposures. This will allow you
;                 to do so if turned on (reduce_std = 1 or
;                 /reduce_std). Default is off since the standard
;                 star(s) will be reduced in the main pipeline anyway. 
;
;   no_pipe    -- to turn off the main reduction pipeline, set no_pipe
;                 = 1 or /no_pipe.
;
;   no_singles -- to only reduce the main exposures, and not any other
;                 exposures for an object, set no_singles = 1 or
;                 /no_singles.  
; 
;   
;
; mage_reduce_all outputs:
;   mage.fits  -- the file name that the mage structure is saved
;                 as. This may be changed, but currently this is the name
;                 the code relies on.
;   
;   mage_night.fits -- copied version of 'mage.fits' with night's date
;   
; MASE outputs:
;   FSpec is the directory that contains the fully reduced spectra of
;   each object for the given night
; 
;   Trace Orders -- Orders.fits and OStr_mage.fits
;
;   Flats -- piximg_flats.fits, Orders.fits (overwrite), Flat
;            (directory), Illumflat_0.70.fits, and Pixflat.fits
;
;   Sensitivity Function -- FSpec (directory), Arcs (directory) and
;                           corresponding files, Final (directory) and
;                           standard star file, Object (directory)
;                           which contains std star structure,
;                           tmp.idl, and std.sav
;  
;   Pipeline -- reduced files are placed in FSpec, while other files
;               are placed in Arcs, Final, and Object directories as
;               the reduction pipeline is running
;
; History:
;  27 Nov 2015 -- created by KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Paste mage_rdconfig.pro functions here later and then
;; remove mage_rdconfig from git repo

function mage_rdconfig_range, val, count=count
  ;; Evaluate val (string). If a single number, return the
  ;; integer. Given a hyphenated range, return expanded array of
  ;; integers.

  prs = strsplit(val,'-',/extract,count=count)

  case count of
     1: exp_arr = [fix(val)]    ; array type
     2: begin
        prs = fix(prs)
        count = prs[1] - prs[0] + 1 ; over-write strsplit() return
        exp_arr = fix(prs[0]) + indgen(count)
     end
     else: begin
        print, 'mage_rdconfig_range(): given val improperly formatted:'
        print, val
        print,'     returning -1'
        exp_arr = -1
        count = 0               ; flag
     end
  endcase                       ; count
  
  return, exp_arr
end                             ; mage_rdconfig_range()


function mage_rdconfig_parse, val, count=count
  ;; Evaluate val (string). Handle commas and call
  ;; mage_rdconfig_range()

  ;; Check for commas first
  prs = strsplit(val,',',/extract,count=nprs)
  for jj=0,nprs-1 do begin
     tmp = mage_rdconfig_range(prs[jj],count=num)
     if num eq 0 then $
        stop,'mage_rdconfig_parse() stop: invalid return',prs[jj]
     
     if jj eq 0 then arr = tmp $ ; instantiate
     else arr = [arr,tmp]        ; append
  endfor                         ; loop jj=nprs
  count = (size(arr,/dim))[0]

  return, arr
end                             ; mage_rdconfig_parse()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mage_rdconfig, config_fil, nexp=nexp, debug=debug
  
  ;; Define structure; may be appended with more std# and fluxtab#
  tags_req = ['dir',$ ; string; full path to night's directory (assume Raw/ underneath)
              'illum',$ ; int arrays of arbitrary length; exposure number(s) for flat files
              'blue',$
              'green',$
              'red',$
              'trace_arc',$     ; example arc file for trace
              'std1',$          ; exposure number(s) of standard star
              'fluxtab1'] ; string; $XIDL_DIR/Spec/Flux/ flux table; may have std2, fluxtab2, etc.
  ntags = n_elements(tags_req)
                   
  if not keyword_set(config_fil) then begin
     if not keyword_set(nexp) then nexp = 50 ; assume won't be more than 50 exposures
     exp_arr = replicate(-1,nexp)
     ;; Make an empty structure
     for tt=0,ntags-1 do begin
        case tags_req[tt] of
           'dir': begin
              ;; handle string
              if not keyword_set(cfg) then $
                 cfg = create_struct(tags_req[tt],'') $      ; create
              else cfg = create_struct(cfg, tags_req[tt],'') ; append
           end
           'fluxtab1': begin
              ;; handle string
              if not keyword_set(cfg) then $
                 cfg = create_struct(tags_req[tt],'') $      ; create
              else cfg = create_struct(cfg, tags_req[tt],'') ; append
           end
           else: begin
              ;; handle strings
              if not keyword_set(cfg) then $
                 cfg = create_struct(tags_req[tt],exp_arr) $      ; create
              else cfg = create_struct(cfg, tags_req[tt],exp_arr) ; append
           end
        endcase
     endfor                     ; loop tt=ntags
     return, cfg                ; empty
  endif

  ;; Read in configuration file 
  readcol,config_fil,descript,val,format='a,a',delimiter='|',/silent
  descript = strtrim(descript,2)
  val = strtrim(val,2)
  nval = n_elements(val)

  ;; Check for un-used config structure elements
  mask = intarr(nval) ; 0 = unused; 1 = used
  
  ;; Loop to instantiate config struct
  for tt=0,ntags-1 do begin
     mtch = where(tags_req[tt] eq strlowcase(descript),nmtch)

     ;; Basic checks
     if nmtch eq 0 then begin
        print,'mage_rdconfig(): required value not in config file ',tags_req[tt]
        ;; Must add and set to a default flag (-1)
        if not keyword_set(cfg) then $
           cfg = create_struct(tags_req[tt],-1) $     ; create 
        else cfg = create_struct(cfg,tags_req[tt],-1) ; append
        continue 
     endif 
     
     if nmtch gt 1 then $
        stop,'mage_rdconfig() stop: multiple values in config file ',tags_req[tt]
     
     mtch = mtch[0]
     mask[mtch] = 1 ; used

     ;; Parse the line:
     if strlowcase(tags_req[tt]) eq 'dir' or $
        stregex(tags_req[tt],'fluxtab1',/boolean,/fold_case) then $
           ;; string (other fluxtab* added later)
        param = val[mtch] $
     else $
        ;; Otherwise a single, list, and/or range of exposure numbers.
        param = mage_rdconfig_parse(val[mtch],count=num)

     if not keyword_set(cfg) then $
        cfg = create_struct(tags_req[tt],param) $     ; create
     else cfg = create_struct(cfg,tags_req[tt],param) ; append
     
  endfor                        ; loop tt=ntags
  
  ;; Figure out what's in config_fil that hasn't been used and add,
  ;; including std#, fluxta# (# > 1)
  bd = where(mask eq 0,nbd)
  for ii=0,nbd-1 do begin
     mtch = bd[ii]
     
     ;; Is it a string?
     if stregex(val[mtch],'[a-z]',/boolean,/fold_case) then begin ; any letter
        ;; Add verbatim
        cfg = create_struct(cfg,descript[mtch],val[mtch])
     endif else begin
        ;; Otherwise parse and add appropriate
        arr = mage_rdconfig_parse(val[mtch],count=num)
        cfg = create_struct(cfg,descript[mtch],arr)
     endelse
     
  endfor                        ; loop ii=nbd

  return, cfg
end                             ; mage_rdconfig()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mage_run_sensfunc,mage,table_path,rawpath,_extra=extra
;reduce_std=reduce_std,no_singles=no_singles, clobber=clobber ; KLC: added this

  if n_params() ne 3 then begin
     print,'Syntax -- mage_run_sensfunc, mage, table_path, rawpath, '
     print,'                  [/reduce_std, /no_singles, /clobber]'
     return 
  endif
  
  ;; Standard file frame
  istd = where(stregex(mage.exptype,'STD',/boolean)) ; may be multiples
  
  ;; Define flux table path
  stdoutfile = mage[0].sensfunc ; default first standard

  ;; Define the sensitivity function's file name
  stdoutfile_new = 'std_'+strtrim(mage[istd[0]].object,2)+'.sav'


;Run sens func
  mage_sensfunc, mage, rawpath+mage[istd[0]].FITSFILE, getenv('XIDL_DIR') + '/Spec/Flux/'+table_path, stdoutfile, CLOBBER = CLOBBER
;  mage_sensfunc, mage, stdfile[0], stdflxfile, stdoutfile, CLOBBER = CLOBBER

  ;; Check if std.sav has been copied to std_star.sav and if not, copy
  test = file_test(stdoutfile_new)
  IF NOT KEYWORD_SET(test) THEN BEGIN
     file_copy,stdoutfile,stdoutfile_new
  ENDIF

  ;; Reduce Standard Star and Exposures (optional)
  IF KEYWORD_SET(reduce_std) THEN BEGIN
     mage_pipe,mage[istd],singles = keyword_set(no_singles) ? 0:1
  ENDIF
END                             ; mage_run_sensfunc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mase_non_gui, config_fil, no_edit=no_edit, no_trace=no_trace, $
                     no_flat=no_flat, chk_flat=chk_flat, $
                     reduce_std=reduce_std, no_sens=no_sens, $
                     no_pipe=no_pipe, no_singles=no_singles, $
                     chk_trc=chk_trc, clobber=clobber 

  if n_params() ne 1 then begin
     print,'Syntax -- mage_reduce_all, config_fil, [/no_edit, /no_trace,'
     print,'          /no_flat,/chk_flat,/reduce_std,/no_sens,'
     print,'          /no_pipe,/no_singles,/chk_trc]'
     return 
  endif

  if size(config_fil,/type) eq 7 then $ ;type 7 is structure
     cfg = mage_rdconfig(config_fil) $ ; read config file
  else cfg = config_fil 
  rawpath = cfg.dir             ; define in config_fil
  
;; Make the mage structure
  magestrct_fil = 'mage.fits'   ; default name
  IF NOT KEYWORD_SET(no_edit) THEN BEGIN
     ;; Create default array of mage structures (one per exposure)
     ;; with best-guess values
     mage_mkstrct, mage, rawpath = rawpath 

;; Edit Structure
     ;; Change Object Name if config struct has nc_* param
     tags=tag_names(cfg)
     gd=where(stregex(tags,'nc_',/boolean,/fold_case),ngd)
     prs=strsplit(tags[gd],'_',/extract,count=nprs)
     if ngd GE 1 then begin
        ;;if object name will be two words, we will have to deal with
        ;;the space in between them, so that is what I am doing inside
        ;;the if statement below
        for i = 0, ngd - 1 do begin ;nprs may be an array
           
           ;;prs may be an array
           if ngd GT 1 then par = prs[i] $
           else par = prs

           if nprs[i] eq 2 then begin            ;ex: nc_ J015650g ThAr
              idx = where(mage.frame EQ cfg.(gd[i])[0])
              mage[idx].object = par[1]
           endif else begin
              name = strjoin(par[1:*],' ') ; par[1]+' '+par[2]
              idx = where(mage.frame EQ cfg.(gd[i])[0])
              mage[idx].object = name
           endelse
       endfor
     endif
     
     ;; Change EXPTYPE (TRASH,SCIENCE,ARC,XE-FLASH,DOMEFLT,BRIGHT,STD)
     ;; Parameters in config_fil must satisfy names in exptyp_par

     exptyp_par = ['trash','science','arc','xe_flash','domeflt','bright']
     ;; Include any standards
     istd = where(stregex(tags,'std',/boolean,/fold_case),nstd) ; id std* tag elements
     if nstd eq 0 then begin
        ;; Should not get here because as of this writing,
        ;; mage_rdconfig() adds STD tag even if not provided in config
        ;; file. It is good to have this message prepared anyway
        print, 'mage_reduce_all: STD tag not found in config structure. Exiting.'
        return
     endif else $
        exptyp_par = [exptyp_par,strlowcase(tags[istd])] ; add 'std1', 'std2', etc
     nexptyp_par = n_elements(exptyp_par)
     for pp=0,nexptyp_par-1 do begin
        mtch = where(stregex(tags,exptyp_par[pp],/boolean,/fold_case),nmtch)      
        if nmtch eq 0 then begin
           continue             ; b/c don't have to change mage.exptype for any
        endif 
        
        ;; nmtch = 1 for the remainder
        itag = mtch[0]
        nexp = n_elements(cfg.(itag)) ; variable length (may include [-1])
        for ee=0,nexp-1 do begin
           indx = where(mage.frame eq cfg.(itag)[ee],nindx)
           if nindx eq 0 then $ ;begin
              stop,'mage_reduce_all stop: No frame #'+strtrim(cfg.(itag)[ee],2)+$
                   ' in mage structure for parameter '+tags[itag]
           ;; Generically handle underscores-should-be-hyphens (e.g., XE-FLASH)
           if stregex(exptyp_par[pp],'_',/boolean) then begin
              prs = strsplit(exptyp_par[pp],'_',/extract)
              val = strjoin(prs,'-')
           endif else begin
              if stregex(exptyp_par[pp],'std',/boolean) then $
                 val = strmid(exptyp_par[pp],0,3) $ ; trim to STD
              else val = exptyp_par[pp]
           endelse              ; not 'xe_flash'
           mage[indx].exptype = strupcase(val)
        endfor                  ; loop ee=nexp
     endfor                     ; loop pp=nexptyp_par

     ;; define mage structure parameters 
     istd = where(stregex(mage.exptype,'STD',/boolean))
     mage.stdfile        = 'Object/ObjStr'+string(mage[istd[0]].frame,format='(i04)')+'.fits' ; KLC thinks 
     mage.pixflatfile    = 'Flat/Pixflat.fits'
     mage.illumflatfile  = 'Flat/Illumflat_' $
;                        + strcompress(mage[istd[0]].SLIT, /rem) + '.fits'
                           + strcompress(mage.SLIT, /rem) + '.fits'
     mage.orderfile      = 'OStr_mage.fits'
     mage.slitfile       = 'Orders.fits'
     mage.sensfunc       =  'std.sav'
     
     ;; Write Structure to magestrct_fil (default: 'mage.fits')
     mwrfits, mage, magestrct_fil, /create 

     ;;Copy structure to night's name
     prs = strsplit(rawpath,'/',/extract,count=nprs)
     dir = prs[nprs-2]                 ; prints 11sep12_adh
     raw_dir = strmid(prs[nprs-2],0,7) ; prints 11sep12
     prs = strsplit(magestrct_fil,'.',/extract,count=nprs) ; assumes root.fits
     newstrctname = prs[0]+'_'+raw_dir+'.'+prs[1]          ; prs[0]='mage' and prs[1]='fits'

     ;;test to see if 'mage.fits' is already copied to newstrctname
     ;;(KLC: clobber business)
     test = file_test(newstrctname)
     IF NOT KEYWORD_SET(test) THEN $
        file_copy,magestrct_fil,newstrctname
   
  ENDIF ELSE BEGIN
     ;; must read in magestrct_fil (default: 'mage.fits') to be mage structure
     test = file_test(magestrct_fil)
     IF NOT KEYWORD_SET(test) THEN $
        stop,'mage_reduce_all stop: default structure '+magestrct_fil+' DNE'
     mage = xmrdfits(magestrct_fil,1,/silent)
  ENDELSE    
 

;; Flats -- Define these file numbers in the config_fil (cfg.illum,
;;          cfg.green, etc)
  
  ;; Illumination          
  illumflatfiles = cfg.dir + 'mage' + string(cfg.illum,format='(i04)')+'.fits'

  ;; GREEN: Xe Flash 5.0" slit
  greenflatfiles = cfg.dir + 'mage' + string(cfg.green,format='(i04)')+'.fits'

  ;; RED: Quartz 5.0" slit
  redflatfiles = cfg.dir + 'mage' + string(cfg.red,format='(i04)')+'.fits'

  ;; BLUE: Twilight flats 5.0" slit
  blueflatfiles = cfg.dir + 'mage' + string(cfg.blue,format='(i04)')+'.fits'

  ;; Assign arc
  pixarcfile = cfg.dir + 'mage' + string(cfg.trace_arc,format='(i04)')+'.fits'


;; Trace slits -- for option to turn off trace, see header
  IF NOT keyword_set(no_trace) THEN BEGIN
     trace_files = [illumflatfiles[0],redflatfiles[0],greenflatfiles[0]]
     tset_slits = mage_findslits(trace_files, no_chk_trc = keyword_set(chk_trc) ? 0:1)
  ENDIF


;; Make Flat -- for option to turn off flats, see header
  IF NOT KEYWORD_SET(no_flat) THEN BEGIN
     mage_makeflat, red = redflatfiles, green = greenflatfiles $
                    , blue = blueflatfiles, illum = illumflatfiles $
                    , orders = mage[0].slitfile, arcfile = pixarcfile[0] $ ; 'Orders.fits' better be set the same for all
                    , chk = chk_flat ;, allchk = flatinter
  ENDIF

;; Sensitivity Function
  IF NOT KEYWORD_SET(no_sens) THEN BEGIN
     mage_run_sensfunc, mage,  cfg.fluxtab1, rawpath, clobber=clobber, $ ; KLC: added clobber= but will be linked to clobber used below
                        reduce_std=keyword_set(reduce_std),no_singles=no_singles
  ENDIF


;; Pipeline
  IF NOT KEYWORD_SET(no_pipe) THEN BEGIN
    mage_pipe, mage, clobber = clobber, singles = (keyword_set(no_singles) ? 0:1)
  ENDIF
END



