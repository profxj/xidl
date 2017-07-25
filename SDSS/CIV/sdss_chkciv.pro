;+ 
; NAME:
; sdss_chkciv
;    Version 2.0
;
; PURPOSE:
;   Visually check SDSS CIV absorption detected with sdss_civsearch
;
; CALLING SEQUENCE:
;   
;   sdss_chkciv, qalfil, civfil, dblt_name=, /sav_ech
;                XSIZE=, YSIZE=, /debug
;
; INPUTS: 
;   qalfil -- FITS file output from sdss_civsearch {sdsscivstrct}
;   civfil -- name of output non-rejected candidates
;
; RETURNS: 
;
; OUTPUTS: 
;   civfil -- FIT file of all CIV candidate systems flagged as
;             DEFINITE (4), GOOD (3), MAYBE (2), BAD (0), unrated
;             (-1); see sdss_getrating()
;
; OPTIONAL KEYWORDS:
;   dblt_name= -- string of doublet to focus on (default 'CIV'; see
;                what dblt_retrieve() accepts)
;   /multi -- input civfil is mix of doublets; assume current doublet
;             is defined by wrest[idx_dblt]
;   /hidealt -- Don't make bottom brackets indicating alternate (wide)
;               doublets in main GUI.
;   altdblt_name= -- string of doublet to be the "if not dblt then
;                    altdblt" (default: 'MgII' if dblt_name='CIV' and
;                    visa versa)
;   XSIZE= -- Size of gui in screen x-pixels (default = 1000)
;   YSIZE= -- Size of gui in screen y-pixels (default = 600)
;   dvqso= -- Buffer to mark before zqso
;   /label -- print rating to top screen of GUI
;   initials= -- string array prepended to corresponding numbers in
;                civstr.rating[] 
;   /sav_ech -- save every step (default: save every 100; faster not
;               to save every step)
;   /debug -- print lots of messages
;   /svrtg -- if set, then do not flip rating[0] to -1
;   mlsfil= -- sdsscivstrct of other lines to show in GUI
;   inspec= -- inputing simple structure inspec={wave:wv,flux:fx,error:er}
;   idblt= -- index in sdsscivstrct arrays of doublet of interest
;             (default: 0; passed via _extra to sdss_chkciv_icmmn)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: 
;    To be called after sdss_civsearch.pro.  Can input output civfil
;    into sdss_fixciv.pro 
;
; EXAMPLES:
;   sdss_chkciv, 'sdss_civsearch.fit', 'good_civ.fit'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Dec-2003 Written by GEP/SHF
;   09-May-2011  Exit GUI more gracefully, continually save, alter so
;                can be used prior to sdss_civsearch and before 
;                sdss_fixciv, MMK
;   31-May-2011  Revamp to new structures and generic for all
;                doublets,/ KLC
;   02-Jun-2011  Elminate dividing input b/c safer, KLC
;   07-Jun-2011  Enable toggling BALFLG, KLC
;   14-Jul-2011  Add BLEND flag to civstr.rating[9], KLC
;   21-Jul-2011  Enable external note capabilities
;   10-Nov-2011  Allow input of structures for viewing, KLC
;   17-Nov-2011  Change how blend in rating[9] flagged, KLC
;   22-Apr-2016  Enable /multi option, KLC
;   10-Aug-2016  Enable locking on single aritrary spectrum, KLC
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;
;; Initialize the common blcok
pro sdss_chkciv_icmmn, state, fail=fail, idblt=idblt, inspec=inspec
  ;; The reason this subroutine exists is because IDL needs the common
  ;; block _definition_ to occur before any reference call. So if this
  ;; was in sdss_chkciv ("main") then all the subroutins before it
  ;; (sdss_chkciv_*) would not be able to refer to the common block.
  ;; It seems like the common block can be accomplished with
  ;; everything just in the state structure (though then
  ;; sdss_chkciv_event couldn't work).


  ;; This replaces an @sdss_functions.  Those files are stored 
  ;; in a separate.pro file, so the following command tells this file
  ;; to treat the following functions as to be resolved later, rather 
  ;; than treating them (erroneously) as floating pt arrays.
  forward_function sdss_getsdssdir, sdss_getrating, sdss_getbalflg, $
     sdss_getcflg, sdss_getrateddblt, sdss_chkciv_ew

  common sdss_chkciv_cmm, $
     clr, $
     npix, $
     fx, $
     fxnorm, $
     wv, $
     sig, $
     signorm, $
     clrsig, $
     fx_lin, $                  ; contrived to show detected centroids
     wv_lin, $
     cindx, $
     conti, $                   ; array of [npix,3] like abslin structure
     clrconti, $                ; array of 3 for colors
     nameconti, $               ; array of 3 for name of conti
     fit, $
     qalstr, $                  ; input
     civstr, $                  ; output
     loscand, $                 ; indices of candidates in current LOS
     nloscand, $
     mlsstr, $                  ; metal lines (firmly identified)
     losmls, $                  ; indices of metal lines
     nlosmls, $
     sdssdir, $
     cinv, $
     idx_dblt                   ; sdsscivstrct index to make primary
  cinv = 1./2.998e5             ; km^-1 s
  clr = getcolor(/load)
  clrconti = [clr.magenta, clr.cyan, clr.limegreen]
  nameconti = ['Eigen','Spline','Hybrid']
  
  if keyword_set(idblt) then idx_dblt = idblt $
  else idx_dblt = 0             ; default

  ;; Since everything set up referenced to $SDSSDIR just store it
  sdssdir = sdss_getsdssdir()
  fail = 0

  ;; QALSTR
  qalstr = xmrdfits(state.qalfil, 1, /silent) ; candidates with sdsscivstrct format
  if idx_dblt ge n_elements(qalstr[0].wrest) then $
     stop,'sdss_chkciv_icmmn stop: idblt = ',idx_dblt,' not allowed; reset idx_dblt'


  ;; metal-line system structure
  if state.mlsfil ne '' then $
     mlsstr = xmrdfits(state.mlsfil, 1, /silent) $
  else mlsstr = 0               ; overwrite anything previous in common block

  if not keyword_set(state.save_ratings) then begin
     qalstr.rating[0] = sdss_getrating(/unrated) ; rating and blend flag storage locales (rating[9] = blend)
     qalstr.notes = ''
     qalstr.rating[9] = qalstr.rating[9] > sdss_getblendflg(/unblend)
  endif 
  state.nqal = n_elements(qalstr)

  state.curqal = 0L

  ;; CIV struct
  test = file_search(state.civfil+'*',count=ntest) ; x_chkfil() slow
  if ntest ge 1 then begin
     print,''
     print,'sdss_chkciv_icmmn: Start from existing '+state.civfil+'?: y(es), n(o).'
     done = 0
     while not done do begin
        ch = get_kbrd(1)        ; w/(1): wait for character
        case byte(strlowcase(ch)) of
           110: begin           ; n
              ;; Make backup and start afresh
              spawn, "mv "+ test + ' ' + test+'.bkp'
              print,'                   Saved to '+test+'.bkp'
              civstr = qalstr
              done = 1
           end
           121: begin           ; y
              ;; Read in civfil and do some checks
              civstr = xmrdfits(state.civfil, 1, /silent) 
              nciv = n_elements(civstr)
              if nciv ne state.nqal then begin
                 print,'sdss_chkciv_icmmn: '+state.civfil+$
                       ' must be same size as '+state.qalfil+'; exiting.'
                 fail = 1
                 return
              endif 
              bd = where(strtrim(qalstr.sdss_obs[0],2) ne $
                         strtrim(civstr.sdss_obs[0],2))
              if bd[0] ne -1 then begin
                 print,'sdss_chkciv_icmmn: '+state.civfil+$
                       ' must be in same order as '+state.qalfil+'; exiting.'
                 fail = 1
                 return
              endif 
              ;; Find first un-rated one and start there
              gd = where(civstr.rating[0] eq sdss_getrating(/unrated),ngd,complement=rtd)
              if ngd eq 0 then $
                 print,'sdss_chkciv_icmmn: Everything rated; starting over.' $
              else state.curqal = gd[0] 
              if keyword_set(state.debug) then $
                 print,'sdss_chkciv_icmmn debug: number of unrated in civfil:',ngd
              done = 1
           end 
           else: print,'sdss_chkciv_icmmn: invalid option: ',ch
        endcase                 ; ch = get_kbrd()
     endwhile                   ; wait for done
     print,''
  endif else begin
     ;; civfil DNE so create new structure array
     civstr = qalstr
  endelse

  ;; Now that state.curqal settled, can assess whether to change
  ;; state.dblt; still necessary for sdss_ewciv() call
  if keyword_set(state.multi) then begin
     state.dblt = dblt_retrieve(qalstr[state.curqal].wrest[idx_dblt])
     if keyword_set(state.debug) then $
        print,'sdss_chkciv_icmmn debug: doublet now ',state.dblt.ion
  endif
  
  tags = tag_names(civstr)
  chk = where(tags eq 'NOTES')
  if chk[0] eq -1 then civstr = sdss_cpstrct(civstr,{sdsscivstrct},excl_tag=['NOTES'])
  civstr.notes = strtrim(civstr.notes,2)

  ;; Set up first value
  if keyword_set(state.debug) then $
     print,'sdss_chkciv_icmmn debug: first state.curqal=',strtrim(state.curqal,2)

  if keyword_set(inspec) then begin
     ;; Handle input of single spectrum
     ;; (snglspec=1/state.flg_snglspec=1)
     ;; Assume about continuum
     cindx = 0
     wv = inspec.wave
     npix = n_elements(wv)
     if npix ne n_elements(inspec.flux) or $
        npix ne n_elements(inspec.error) then $
           stop,'hst_chkciv_icmmn stop: number of elements in inspec.wave/flux/error must match'
     fx = inspec.flux
     fxnorm = fx
     sig = inspec.error
     signorm = sig
     conti = replicate(1,npix,3)
     ;; Other common variables sdss_chkciv_setup would handle
     ;; otherwise
     state.filnam = strtrim(qalstr[state.curqal].sdss_obs[0],2)
     ;; Set common QSO info
     state.zqso = qalstr[state.curqal].z_qso
     state.ra = qalstr[state.curqal].ra
     state.dec = qalstr[state.curqal].dec
     state.mag = qalstr[state.curqal].rmag
     ;; Figure out the other candidates in sightline
     loscand = where(strtrim(qalstr[state.curqal].qso_name,2) eq $
                     strtrim(qalstr.qso_name,2),nloscand)
     ;; Figure out other metal-line systems in sightline
     if keyword_set(mlsstr) then $
        losmls = where(strtrim(qalstr[state.curqal].qso_name,2) eq $
                       strtrim(mlsstr.qso_name,2),nlosmls) $
     else nlosmls = 0
     filnam = strtrim(qalstr[state.curqal].abslin_fil,2)
     test = file_search(sdssdir+filnam+'*',count=ntest)
     if ntest eq 0 then $       ; try locally
        test = file_search(filnam+'*',count=ntest)
     if ntest gt 0 then begin
        cstrct = xmrdfits(test[0],1,/silent)
        cindx = fix(alog(qalstr[state.curqal].cflg)/alog(2))
        ;; Set contrived "spectrum" to plot all detected lines, making a
        ;; square waveform
        sdss_pltabslin, fx, cstrct, /noplt, wave=wv, fx_lin=fx_lin, $
                        wv_lin=wv_lin
     endif else begin
        print,'sdss_chkciv_icmmn: abslin file DNE.'
        wv_lin = [0.,0.]
        fx_lin = [0.,0.]
     endelse 
  endif


end                             ; sdss_chkciv_icmmn

;;;;;;;;;;;;;;;;;;;;;;
;; Events
pro sdss_chkciv_event, ev

  common sdss_chkciv_cmm

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  if state.curqal lt state.nqal then $
     civstr[state.curqal].notes = strtrim(civstr[state.curqal].notes,2)

  case uval of
     ;; Ratings
     'DEFINITE': begin
        sdss_chkciv_setciv, state, /definite
        sdss_chkciv_next, state, KILL=kill
        if kill eq 0 then begin
           sdss_chkciv_setup, state
           sdss_chkciv_update, state
        endif ; else pause
     end
     'GOOD': begin
        sdss_chkciv_setciv, state, /good
        sdss_chkciv_next, state, KILL=kill
        if kill eq 0 then begin
           sdss_chkciv_setup, state
           sdss_chkciv_update, state
        endif ; else pause
     end

     'MAYBE': begin
        sdss_chkciv_setciv, state, /maybe
        sdss_chkciv_next, state, KILL=kill
        if kill eq 0 then begin
           sdss_chkciv_setup, state
           sdss_chkciv_update, state
        endif ; else pause
     end
     'BAD': begin 
        sdss_chkciv_setciv, state, /bad
        sdss_chkciv_next, state, KILL=kill
        if kill eq 0 then begin
           sdss_chkciv_setup, state
           sdss_chkciv_update, state
        endif ; else pause
     end
     ;; Movement
     'DONE' : begin
        sdss_chkciv_svciv, state
        widget_control, ev.top, /destroy
        return
     end
     'BACK': begin 
        sdss_chkciv_next, state, /back
        sdss_chkciv_setup, state
        sdss_chkciv_update, state
     end
     'NEXT': begin              ; increment of one; no grading
        sdss_chkciv_next, state, kill=kill
        if kill eq 0 then begin
           sdss_chkciv_setup, state
           sdss_chkciv_update, state
        endif ; else pause
     end
     'SKIP': begin              ; Go to next unrated, no grading 
        sdss_chkciv_next, state, /skip, kill=kill
        if kill eq 0 then begin
           sdss_chkciv_setup, state
           sdss_chkciv_update, state
        endif ; else pause
     end

     ;; niceties
     'SAVE' : begin
        sdss_chkciv_svciv, state
     end
     'CHGBAL': begin            ; set BAL flag
        sdss_chkciv_chgbalflg, state
        sdss_chkciv_setup, state
        sdss_chkciv_update, state
     end
     'BLEND': begin             ; set blended flag (rating[9])
        sdss_chkciv_blend, state
        sdss_chkciv_setup, state
        sdss_chkciv_update, state        
     end
     'CURNUM': begin            ; jump to number
        widget_control, state.num_id, get_value=tmp
        if tmp lt state.nqal then begin
           state.curqal = tmp 
           ;; Now that state.curqal settled, can assess whether to change
           ;; state.dblt
           if keyword_set(state.multi) then begin
              state.dblt = dblt_retrieve(qalstr[state.curqal].wrest[idx_dblt])
              if keyword_set(state.debug) then $
                 print,'sdss_chkciv_event debug: doublet now ',state.dblt.ion
           endif              
        endif else print,'sdss_chkciv_event: Index beyond number of candidates.'
        sdss_chkciv_setup, state ; revert number
        sdss_chkciv_update, state        
     end 
     'NOTES': begin             ; save notes
        widget_control, state.note_id, get_value=tmp
        civstr[state.curqal].notes = tmp
     end 
     'SPLT': x_specplot, fx, sig, ytwo=conti[*,cindx], $
                         psym2=-3, wave=wv, zin=state.zabs, /lls, $
                         /block, inflg=4, xrange=qalstr[state.curqal].wrest[idx_dblt:idx_dblt+1]*$
                         (1.+state.zabs)+[-100,100],$
                         ythree=fx_lin, three_wave=wv_lin, $
                         title=qalstr[state.curqal].qso_name+$
                         string(qalstr[state.curqal].z_qso,format='(1x,"zqso = ",f7.5)')
     state.dbltalt.ion: begin   ; MgII overplot
;        zalt = state.dblt.wvI/state.dbltalt.wvI*(1.+state.zabs) ; 1 + z
        zalt = qalstr[state.curqal].wrest[idx_dblt]/state.dbltalt.wvI*(1.+state.zabs) ; 1 + z
        x_specplot, fx, sig, ytwo=conti[*,cindx], $
                    psym2=-3, wave=wv, zin=zalt, /lls, /block, $
                    inflg=4, xrange=[state.dbltalt.wvI,state.dbltalt.wvII]*$
                    zalt+[-100,100],$
                    ythree=fx_lin, three_wave=wv_lin, $
                    title=qalstr[state.curqal].qso_name+$
                    string(qalstr[state.curqal].z_qso,format='(1x,"zqso = ",f7.5)')
     end

     ;; Shortcuts
     'TEXT' : begin             ; almost all left handed keys
        eventch = string(ev.ch)
        case eventch of
           '3': begin                       ; DEFINITE
              if (ev.press eq 1) then begin ; this check prevents double click
                 sdss_chkciv_setciv, state, /definite
                 sdss_chkciv_next, state, KILL=kill
                 if kill eq 0 then begin
                    sdss_chkciv_setup, state
                    sdss_chkciv_update, state
                 endif ; else pause
              endif 
           end
           '2': begin           ; GOOD
              if (ev.press eq 1) then begin
                 sdss_chkciv_setciv, state, /good
                 sdss_chkciv_next, state, KILL=kill
                 if kill eq 0 then begin
                    sdss_chkciv_setup, state
                    sdss_chkciv_update, state
                 endif ; else pause
              endif 
           end
           '1': begin           ; MAYBE
              if (ev.press eq 1) then begin
                 sdss_chkciv_setciv, state, /maybe
                 sdss_chkciv_next, state, KILL=kill
                 if kill eq 0 then begin
                    sdss_chkciv_setup, state
                    sdss_chkciv_update, state
                 endif ; else pause
              endif 
           end
           '`': begin           ; BAD (sneaky! to be left-hand)
              if (ev.press eq 1) then begin
                 sdss_chkciv_setciv, state, /bad
                 sdss_chkciv_next, state, KILL=kill
                 if kill eq 0 then begin
                    sdss_chkciv_setup, state
                    sdss_chkciv_update, state
                 endif ; else pause
              endif 
           end
           '0': begin           ; BAD (to match rating)
              if (ev.press eq 1) then begin
                 sdss_chkciv_setciv, state, /bad
                 sdss_chkciv_next, state, KILL=kill
                 if kill eq 0 then begin
                    sdss_chkciv_setup, state
                    sdss_chkciv_update, state
                 endif ; else pause
                 
              endif
           end
           'b': begin           ; BACK
              if (ev.press eq 1) then begin
                 sdss_chkciv_next, state, /back
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end	
           'n': begin           ; Next (one)
              if (ev.press eq 1) then begin
                 sdss_chkciv_next, state, kill=kill
                 if kill eq 0 then begin
                    sdss_chkciv_setup, state
                    sdss_chkciv_update, state
                 endif ; else pause
              endif 
           end	
           's': begin           ; Skip to unrated
              if (ev.press eq 1) then begin
                 sdss_chkciv_next, state, /skip, kill=kill
                 if kill eq 0 then begin
                    sdss_chkciv_setup, state
                    sdss_chkciv_update, state
                 endif ; else pause
                 if kill eq 1 then sdss_chkci_svciv, state
              endif 
           end	
           'q' : begin
              if (ev.press eq 1) then begin
                 sdss_chkciv_svciv, state
                 widget_control, ev.top, /destroy
                 return
              endif 
           end
           
           'w' : begin          ; Save (as in write)
              if (ev.press eq 1) then sdss_chkciv_svciv, state
           end
           'a': begin           ; MgII overplot
              if (ev.press eq 1) then begin
                 zalt = qalstr[state.curqal].wrest[idx_dblt]/state.dbltalt.wvI*(1.+state.zabs) - 1.
                 x_specplot, fx, sig, ytwo=conti[*,cindx], $
                             psym2=-3, wave=wv, zin=zalt, $
                             /lls, /block, inflg=4, $
                             xrange=[state.dbltalt.wvI,state.dbltalt.wvII]*$
                             (1.+zalt)+[-100,100],$
                             ythree=fx_lin, three_wave=wv_lin, $
                             title=qalstr[state.curqal].qso_name+$
                             string(qalstr[state.curqal].z_qso,format='(1x,"zqso = ",f7.5)')
              endif 
           end
           'f': begin           ; BALFLG
              if (ev.press eq 1) then begin
                 sdss_chkciv_chgbalflg, state
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end
           'd': begin           ; BLEND
              if (ev.press eq 1) then begin
                 sdss_chkciv_blend, state
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end
           'x': begin  ; adjust left bound of 1548 EW region
              if (ev.press eq 1) then begin
                 state.xcurs = ev.x
                 xnew = xgetx_plt(state, /strct)
                 if (total(civstr[state.curqal].wvlim_final) EQ 0) then $
                    civstr[state.curqal].wvlim_final = $
                    civstr[state.curqal].wvlim_orig
                 if (xnew * (1+state.zabs) LT $
                     civstr[state.curqal].wvlim_final[idx_dblt,1]) then begin
                    civstr[state.curqal].wvlim_final[idx_dblt,0] $
                       = xnew * (1+state.zabs)
                 endif else begin
                    print, "ERROR: left bound > right bound"
                    break
                 endelse

                 newcivstr = sdss_chkciv_ew(state, civstr[state.curqal])
                 civstr[state.curqal] = newcivstr[0]
                 sdss_chkciv_setewstr,state,civstr[state.curqal].ew_final,$
                                      civstr[state.curqal].sigew_final
                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end
           'l': begin  ; adjust left bound of 1548 EW region
              if (ev.press eq 1) then begin
                 state.xcurs = ev.x
                 xnew = xgetx_plt(state, /strct)
                 if (total(civstr[state.curqal].wvlim_final) EQ 0) then $
                    civstr[state.curqal].wvlim_final = $
                    civstr[state.curqal].wvlim_orig
                 if (xnew * (1+state.zabs) LT $
                     civstr[state.curqal].wvlim_final[idx_dblt,1]) then begin
                    civstr[state.curqal].wvlim_final[idx_dblt,0] $
                       = xnew * (1+state.zabs)
                 endif else begin
                    print, "ERROR: left bound > right bound"
                    break
                 endelse

                 newcivstr = sdss_chkciv_ew(state, civstr[state.curqal])
                 civstr[state.curqal] = newcivstr[0]
                 sdss_chkciv_setewstr,state,civstr[state.curqal].ew_final,$
                                      civstr[state.curqal].sigew_final
                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end
           'X': begin           ; adjust left bound of 1550 EW region
              if (ev.press eq 1) then begin
                 state.xcurs = ev.x
                 xnew = xgetx_plt(state, /strct)
                 if (total(civstr[state.curqal].wvlim_final) EQ 0) then $
                    civstr[state.curqal].wvlim_final = $
                    civstr[state.curqal].wvlim_orig
                 if (xnew * (1+state.zabs) LT $
                     civstr[state.curqal].wvlim_final[idx_dblt+1,1]) then begin
                    civstr[state.curqal].wvlim_final[idx_dblt+1,0] $
                       = xnew * (1+state.zabs)
                 endif else begin
                    print, "ERROR: left bound > right bound"
                    break
                 endelse

                 newcivstr = sdss_chkciv_ew(state, civstr[state.curqal])
                 civstr[state.curqal] = newcivstr[0]
                 sdss_chkciv_setewstr,state,civstr[state.curqal].ew_final,$
                                      civstr[state.curqal].sigew_final
                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end
           'L': begin           ; adjust left bound of 1550 EW region
              if (ev.press eq 1) then begin
                 state.xcurs = ev.x
                 xnew = xgetx_plt(state, /strct)
                 if (total(civstr[state.curqal].wvlim_final) EQ 0) then $
                    civstr[state.curqal].wvlim_final = $
                    civstr[state.curqal].wvlim_orig
                 if (xnew * (1+state.zabs) LT $
                     civstr[state.curqal].wvlim_final[idx_dblt+1,1]) then begin
                    civstr[state.curqal].wvlim_final[idx_dblt+1,0] $
                       = xnew * (1+state.zabs)
                 endif else begin
                    print, "ERROR: left bound > right bound"
                    break
                 endelse

                 newcivstr = sdss_chkciv_ew(state, civstr[state.curqal])
                 civstr[state.curqal] = newcivstr[0]
                 sdss_chkciv_setewstr, state, civstr[state.curqal].ew_final, $
                                       civstr[state.curqal].sigew_final
                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end
           'r': begin           ; adjust right bound of 1548 region
              if (ev.press eq 1) then begin
                 state.xcurs = ev.x
                 xnew = xgetx_plt(state, /strct)
                 if (total(civstr[state.curqal].wvlim_final) EQ 0) then $
                    civstr[state.curqal].wvlim_final = $
                    civstr[state.curqal].wvlim_orig
                 if (xnew * (1+state.zabs) GT $
                     civstr[state.curqal].wvlim_final[idx_dblt,0]) then begin
                    civstr[state.curqal].wvlim_final[idx_dblt,1] $
                       = xnew * (1+state.zabs)
                 endif else begin
                    print, "ERROR: Right bound < left bound"
                    break
                 endelse

                 newcivstr = sdss_chkciv_ew(state, civstr[state.curqal])
                 civstr[state.curqal] = newcivstr[0]
                 sdss_chkciv_setewstr, state, civstr[state.curqal].ew_final, $
                                       civstr[state.curqal].sigew_final
                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end
           'R': begin           ; adjust right bound of 1550 EW region
              if (ev.press eq 1) then begin
                 state.xcurs = ev.x
                 xnew = xgetx_plt(state, /strct)
                 if (total(civstr[state.curqal].wvlim_final) EQ 0) then $
                    civstr[state.curqal].wvlim_final = $
                    civstr[state.curqal].wvlim_orig
                 if (xnew * (1+state.zabs) GT $
                     civstr[state.curqal].wvlim_final[idx_dblt+1,0]) then begin
                    civstr[state.curqal].wvlim_final[idx_dblt+1,1] $
                       = xnew * (1+state.zabs)
                 endif else begin
                    print, "ERROR: Right bound < left bound"
                    break
                 endelse

                 newcivstr = sdss_chkciv_ew(state, civstr[state.curqal])
                 civstr[state.curqal] = newcivstr[0]
                 sdss_chkciv_setewstr, state, civstr[state.curqal].ew_final, $
                                       civstr[state.curqal].sigew_final

                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end                  ; 'R'
           
           'u': begin           ; restore bounds of 1548 line
              if (ev.press eq 1) then begin
                 civstr[state.curqal].zabs_final[idx_dblt] = $
                    qalstr[state.curqal].zabs_orig[idx_dblt]
                 civstr[state.curqal].ew_final[idx_dblt] = $
                    qalstr[state.curqal].ew_orig[idx_dblt]
                 civstr[state.curqal].sigew_final[idx_dblt] = $
                    qalstr[state.curqal].sigew_orig[idx_dblt]
                 civstr[state.curqal].wvlim_final[idx_dblt,*] = $
                    qalstr[state.curqal].wvlim_orig[idx_dblt,*]

                 sdss_chkciv_setewstr, state, civstr[state.curqal].ew_final, $
                                       civstr[state.curqal].sigew_final
                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif 
           end                  ; 'u'
           'U': begin           ; restore bounds of 1550 line
              if (ev.press eq 1) then begin
                 civstr[state.curqal].zabs_final[idx_dblt+1] = $
                    qalstr[state.curqal].zabs_orig[idx_dblt+1]
                 civstr[state.curqal].ew_final[idx_dblt+1] = $
                    qalstr[state.curqal].ew_orig[idx_dblt+1]
                 civstr[state.curqal].sigew_final[idx_dblt+1] = $
                    qalstr[state.curqal].sigew_orig[idx_dblt+1]
                 civstr[state.curqal].wvlim_final[idx_dblt+1,*] = $
                    qalstr[state.curqal].wvlim_orig[idx_dblt+1,*]

                 sdss_chkciv_setewstr, state, civstr[state.curqal].ew_final, $
                                       civstr[state.curqal].sigew_final
                 state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
                 
                 sdss_chkciv_setup, state
                 sdss_chkciv_update, state
              endif    
           end                  ; 'U'

           'i': begin           ; Info about other candidates in sightline
              if (ev.press eq 1) then begin
                 print,''
                 print,'Sightline Candidate Summary'
                 print,'index','zabs','EW','RTG','Notes',$
                       format='(a5,2x,a7,2x,a5,2x,a3,"; ",a)'
                 printcol,loscand,civstr[loscand].zabs_orig[idx_dblt],$
                          civstr[loscand].ew_orig[idx_dblt],$
                          civstr[loscand].rating[0],$
                          civstr[loscand].notes,$
                          format='(i5,2x,f7.5,2x,f5.2,2x,i3,": ",a)'
              endif 
           end                  ; 'i'


           'H': if (ev.press eq 1) then x_helpwidg, state.help
           ' ': if (ev.press eq 1) then $ ; plot
              x_specplot, fx, sig, ytwo=conti[*,cindx], $
                          psym2=-3, wave=wv, zin=state.zabs, $
                          /lls, /block, $
                          inflg=4, xrange=qalstr[state.curqal].wrest[idx_dblt:idx_dblt+1]*$
                          (1.+state.zabs)+[-100,100],$
                          ythree=fx_lin, three_wave=wv_lin, $
                          title=qalstr[state.curqal].qso_name+$
                          string(qalstr[state.curqal].z_qso,format='(1x,"zqso = ",f7.5)')
           't': if (ev.press eq 1) then $ ; plot 2 contis
              x_splot, wv, fx, psym1=10, /block, $
                       ytwo=sig, psym2=10, color2=clrsig, $ ; red
                       ythr=conti[*,sdss_getcflg(/eig,/index)], psym3=-3, $
                       color3=clrconti[sdss_getcflg(/eig,/index)], $
                       xmnx=qalstr[state.curqal].wrest[idx_dblt:idx_dblt+1]*(1.+state.zabs)+[-100,100], $
                       yfou=conti[*,sdss_getcflg(/spl,/index)], psym4=-3, $
                       color4=clrconti[sdss_getcflg(/spl,/index)], $
                       yfiv=conti[*,sdss_getcflg(/hyb,/index)], psym5=-3, $
                       color5=clrconti[sdss_getcflg(/hyb,/index)], $
                       xsix=wv_lin,ysix=fx_lin,psym6=10, color6=clr.blue, $
                       title=qalstr[state.curqal].qso_name+$
                       string(qalstr[state.curqal].z_qso,format='(1x,"zqso = ",f7.5)')+$
                       ': Using '+nameconti[cindx],$
                       lgnd=['Flux','Sigma',nameconti,'Lines']
           '-': if (ev.press eq 1) then sdss_pltsnr, state.filnam, snr=3.5

           ;; Auto-notes
           'E': begin          
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'ghost emission'
           end
           'B': begin         
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'bad conti'
           end
           'A': begin           ; x_specplot is AlIII 
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'AlII'
           end 
           'C': begin           
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'CIV'
           end 
           'D': begin           ; any DLA weirdness
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'DLA ion'
           end 
           'F': begin          
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'FeII' ; any of them
           end 
           'K': begin          
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'sky lines' ; any of them
           end 
           'M': begin       
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'MgII'
           end 
           'S': begin          
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'SiII' ; 1526
           end 
           'V': begin           ; matches x_specplot    
              if (ev.press eq 1) then $
                 sdss_chkciv_addnote, state, 'SiIV'
           end 

           else: ;print, 'sdss_chkciv: Not a valid key!'
           
        endcase                 ; eventch of 
        
     end                        ; 'TEXT' 
     else:                      ; nothing to do I guess
     
  endcase                       ; uval of 
  
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy 
end                             ; sdss_chkciv_event


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Recalcaulate equivalent widths 
;; Use same algorithm used previously (via sdss_ewciv)
function sdss_chkciv_ew, state, civstrct

  common sdss_chkciv_cmm

  ;; Do some checks
  limits_1548 = civstrct.wvlim_final[idx_dblt,*]
  limits_1550 = civstrct.wvlim_final[idx_dblt,*]

  in1548 = where(wv GE limits_1548[0] AND wv LE limits_1548[1], n1548)
  in1550 = where(wv GE limits_1550[0] AND wv LE limits_1550[1], n1550)

  if (n1548 EQ 0 OR n1550 EQ 0) then begin
     print, "ERROR: illegal region chosen?  Exiting"
     return, -1
  endif
  ;; Shouldn't matter putting in old redshift b/c will
  ;; calculate new redshift based on set (and fixed) ewlim_final
  ;; Also sets ncolm_final, etc and flags
  if civstrct.zabs_final[idx_dblt] gt 0. then $
     zabs = civstrct.zabs_final[idx_dblt] else zabs = civstrct.zabs_orig[idx_dblt]
  ;; This is the only place need state.dblt
  newcivstr = sdss_ewciv(wv,fx,sig,conti[*,cindx], state.dblt, zabs, $
                         istrct=civstrct,/final,/keepwvlim,debug=state.debug)
  state.zabs = newcivstr.zabs_final[idx_dblt]
  
  if keyword_set(state.debug) then begin
     ;; Recalculate the redshift for new limits the old way
     tmpciv = civstrct ; this is only one structure
     dwv = wv - shift(wv,1)
     dwv[0] = dwv[1]
     ;; New redshift is the average of both components, not separate.
     z_new = (total((1-fxnorm[in1548])*$
                    (wv[in1548]/tmpciv.wrest[idx_dblt]-1.0)) + $
              total((1-fxnorm[in1550])*$
                    (wv[in1550]/tmpciv.wrest[idx_dblt+1]-1.0))) / $
             (total(1-fxnorm[in1548]) + $
              total(1-fxnorm[in1550]))
     
     zmoment1 = sqrt((total((1-fxnorm[in1548])* $
                            (wv[in1548]/tmpciv.wrest[idx_dblt]-1.0-z_new)^2) + $
                      total((1-fxnorm[in1550])* $
                            (wv[in1550]/tmpciv.wrest[idx_dblt+1]-1.0-z_new)^2)) / $
                     (total(1-fxnorm[in1548]) + $
                      total(1-fxnorm[in1550])))
     
     tmpciv.zabs_final[idx_dblt:idx_dblt+1] = z_new
     state.zabs             = z_new

     ;; Calculate equivalent width
     ew_1548 = total((1-fxnorm[in1548]) * dwv[in1548])
     ew_1550 = total((1-fxnorm[in1550]) * dwv[in1550])
     tmpciv.ew_final[idx_dblt:idx_dblt+1] = [ew_1548, ew_1550] / (1+state.zabs)

     ;; Equivalent width error
     ;; WARNING: Assumes zero error in the continuum!
     sigew_1548 = sqrt(total((signorm[in1548] * dwv[in1548])^2))
     sigew_1550 = sqrt(total((signorm[in1550] * dwv[in1550])^2))
     tmpciv.sigew_final[idx_dblt:idx_dblt+1] = [sigew_1548, sigew_1550]/(1+state.zabs)

     print,'sdss_chkciv_ew() debug: comparing sdss_ewciv and alt measures:'
     print,'wrest','zabs','wvlo','wvhi','EW','sEW','dv',$
           format='(a7,1x,a7,2(1x,a8),2(1x,a6),1x,a4)'
     print,newcivstr.wrest[idx_dblt],$
           newcivstr.zabs_final[idx_dblt],newcivstr.wvlim_final[idx_dblt,0],$
           newcivstr.wvlim_final[idx_dblt,1],newcivstr.ew_final[idx_dblt],$
           newcivstr.sigew_final[idx_dblt],$
           format='(f7.2,1x,f7.5,2(1x,f8.3),2(1x,f6.4))'
     print,tmpciv.wrest[idx_dblt],$
           tmpciv.zabs_final[idx_dblt],tmpciv.wvlim_final[idx_dblt,0],$
           tmpciv.wvlim_final[idx_dblt,1],tmpciv.ew_final[idx_dblt],$
           tmpciv.sigew_final[idx_dblt],$
           format='(f7.2,1x,f7.5,2(1x,f8.3),2(1x,f6.4))'
     print,newcivstr.wrest[idx_dblt+1],$
           newcivstr.zabs_final[idx_dblt+1],newcivstr.wvlim_final[idx_dblt+1,0],$
           newcivstr.wvlim_final[idx_dblt+1,1],newcivstr.ew_final[idx_dblt+1],$
           newcivstr.sigew_final[idx_dblt+1],$
           round((newcivstr.zabs_final[idx_dblt+1]-newcivstr.zabs_final[idx_dblt])/$
                 (1+newcivstr.zabs_final[idx_dblt])*3.e5),$
           format='(f7.2,1x,f7.5,2(1x,f8.3),2(1x,f6.4),1x,i4)'
     print,tmpciv.wrest[idx_dblt+1],$
           tmpciv.zabs_final[idx_dblt+1],tmpciv.wvlim_final[idx_dblt+1,0],$
           tmpciv.wvlim_final[idx_dblt+1,1],tmpciv.ew_final[idx_dblt+1],$
           tmpciv.sigew_final[idx_dblt+1],$
           round((tmpciv.zabs_final[idx_dblt+1]-tmpciv.zabs_final[idx_dblt])/$
                 (1+tmpciv.zabs_final[idx_dblt])*3.e5),$
           format='(f7.2,1x,f7.5,2(1x,f8.3),2(1x,f6.4),1x,i4)'
  endif                         ; /debug

  return, newcivstr

end                             ; sdss_chkciv_ew()


pro sdss_chkciv_setewstr, state, ewarr, sigewarr
  ;; Make strings in state when EW changed
  ;; Ratio is of 1548/1550.

  common sdss_chkciv_cmm  

  state.ew = string(ewarr[idx_dblt],sigewarr[idx_dblt],$
                    format='(f4.2,"+/-",f4.2)')
  ewrto = ewarr[idx_dblt]/ewarr[idx_dblt+1]
  ;; This is the quickie error-propagation formula
  sigewrto = abs(ewrto)*sqrt((sigewarr[idx_dblt]/ewarr[idx_dblt])^2+$
                             (sigewarr[idx_dblt+1]/ewarr[idx_dblt+1])^2)
  state.ewrto = string(ewrto,sigewrto,$
                       format='(f4.2,"+/-",f4.2)')
end                             ;  sdss_chkciv_setewstr


pro sdss_chkciv_addnote, state, note
  ;; Add message and update
  common sdss_chkciv_cmm

  note = strtrim(note,2)
  if civstr[state.curqal].notes eq '' then $
        civstr[state.curqal].notes = note $
  else civstr[state.curqal].notes = $
     civstr[state.curqal].notes+'; ' + note 
  sdss_chkciv_update, state                 
 
end                             ; sdss_chkciv_addnote


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Metals Plot
pro sdss_chkciv_metals, state
  
  common sdss_chkciv_cmm

  ;; Set plot window
  if state.psfile NE 1 then begin
      widget_control, state.mdraw_id, get_value=wind
      wset, wind
  endif

  gdlin = where((state.velplt.flg MOD 2) EQ 1, ngd)

  ny = ngd / 2 + (ngd MOD 2 EQ 1)

  ;; Good lines
  !p.multi = [0,2,ny,0,1]
  for j=0L,ngd-1 do begin
      i = gdlin[j]
      pixmin = state.all_pmnx[0,i]
      pixmax = state.all_pmnx[1,i]
      ymnx = [(min(fxnorm[pixmin:pixmax],max=mx)-0.1) > state.velplt[i].ymnx[0],$
              (mx+0.1) < state.velplt[i].ymnx[1]]           ; customize

      ;; Plot
      if (j NE ny-1) AND (j NE ngd-1 ) then begin
          spaces = replicate('!17 ',30)
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            fxnorm[pixmin:pixmax], xrange=state.vmnx, $
            yrange=ymnx, xtickn=spaces, xmargin=[9,3], $
            ymargin=[0,0], yminor=1,$
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endif else begin
          plot, state.all_velo[0:state.all_pmnx[2,i],i], $
            fxnorm[pixmin:pixmax], xrange=state.vmnx, $
            yrange=ymnx, xmargin=[9,3], ymargin=[3,0], yminor=1, $
            charsize=1.8, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=1
      endelse

      ;; Labels
      xyouts, 0.07*(state.vmnx[1]-state.vmnx[0])+state.vmnx[0], $
        ymnx[0]+ $
        (ymnx[1]-ymnx[0])*0.05, $
        strtrim(state.velplt[i].name,2), $
        color=clr.black, charsize=1.5
      
      ;; Lines
      oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2
      oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3
      oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3
  endfor
  
  ;;;;;;;;;;;;;;
  ;; Top window
  !p.multi = [0,1,1]
  widget_control, state.tdraw_id, get_value=wind
  wset, wind
  gd = where(wv GT (qalstr[state.curqal].wrest[idx_dblt]-150.)*(1.+state.zabs) AND $
             wv LT (1.+state.zabs)*(qalstr[state.curqal].wrest[idx_dblt]+150.),ngd)

  ;; THE MAIN PLOT!!!
  ;; Plot flux +/- total flux error for two continua in rest frame
  ;; +/-6000 kms (all relative in case these options change)
  state.xymnx = [qalstr[state.curqal].wrest[idx_dblt]*(1.-2e-2),-0.2,$
                 qalstr[state.curqal].wrest[idx_dblt]*(1.+2e-2),1.4]
  dx = state.xymnx[2]-state.xymnx[0]
  dy = state.xymnx[3]-state.xymnx[1]
  plot, wv[gd]/(1.+state.zabs), fxnorm[gd], $
        xrange=state.xymnx[[0,2]], $ ;+/-50 used to be
        yrange=state.xymnx[[1,3]], $ ;[-0.2,1.4] used to be
        position=state.pos, $
        charsize=1.8, psym=10, background=clr.white, color=clr.black, $
        xstyle=1, ystyle=1


  oplot, wv[gd]/(1+state.zabs), signorm[gd], $
         color=clr.red, psym=10

  ;; Label with the rating, bottom left
  if keyword_set(state.label) then begin
     if state.ninitials eq 0 then $
        lbl = strtrim(civstr[state.curqal].rating[0],2) $
     else begin
        arr = strtrim(state.initials[0:state.ninitials-1],2) + ':' + $
              strtrim(civstr[state.curqal].rating[0:state.ninitials-1],2)
        lbl = strjoin(arr,'; ')
     endelse 
     xyouts, 0.0, 0.02, lbl, color=clr.red, charsize=3, /normal ; abs coord
  endif 

  ;; Plot doublet bounds; green for 1548, red for 1550
  ;; Fancy if-else 
  wvlim_plt = (total(civstr[state.curqal].wvlim_final) EQ 0) ? $
              qalstr[state.curqal].wvlim_orig : civstr[state.curqal].wvlim_final

  tmpy = state.xymnx[3]-0.0625*dy-[0.125*dy,0] ; [1.1,1.3]
  oplot, wvlim_plt[0,0]/([1.,1.]+state.zabs), tmpy,$
         color=clr.green, linestyle=0, thick=3
  oplot, [wvlim_plt[0,0]/(1.+state.zabs),wvlim_plt[0,1]/(1.+state.zabs)],$
         [tmpy[1],tmpy[1]], color=clr.green, linestyle=0, thick=3
  oplot, wvlim_plt[0,1]/([1.,1.]+state.zabs), tmpy,$
         color=clr.green, linestyle=0, thick=3

  tmpy = tmpy + 0.03125*dy ; [1.15,1.35]
  oplot, wvlim_plt[1,0]/([1.,1.]+state.zabs), tmpy,$
         color=clr.red, linestyle=0, thick=3
  oplot, [wvlim_plt[1,0]/(1.+state.zabs),wvlim_plt[1,1]/(1.+state.zabs)],$
         [tmpy[1],tmpy[1]], color=clr.red, linestyle=0, thick=3
  oplot, wvlim_plt[1,1]/([1.,1.]+state.zabs), tmpy,$
         color=clr.red, linestyle=0, thick=3

  ;; Indicate other candidates in sightline (even if mixed types)
  tmpy[0] += 0.03125*dy ; [1.25,1.35]
  for ii=0,nloscand-1 do begin
     if loscand[ii] eq state.curqal then continue ; don't mark current
     wvlim_plt = qalstr[state.curqal].wrest[idx_dblt:idx_dblt+1]*(1+qalstr[loscand[ii]].zabs_orig[idx_dblt])/$
                 (1+state.zabs)
     oplot, [wvlim_plt[0],wvlim_plt[0]], tmpy,color=clr.blue, thick=2, $
            linestyle=0
     oplot, wvlim_plt, [tmpy[1],tmpy[1]],color=clr.blue, linestyle=0, thick=2
     oplot, [wvlim_plt[1],wvlim_plt[1]], tmpy,color=clr.blue, thick=2, $
            linestyle=0
  endfor 
  
  ;; Indicate other metal lines in sightline
  tmpy = tmpy - 0.015625*dy     ; [1.225,1.325]
  for ii=0,nlosmls-1 do begin
     if losmls[ii] eq state.curqal then continue ; don't mark current
     ;; Do *not* want wrest[idx_dblt:idx_dblt+1]; assume mlsstr has
     ;; primary doublet in first elements
     wvlim_plt = mlsstr[losmls[ii]].wrest[0:1]*(1+mlsstr[losmls[ii]].zabs_orig[0])/$
                 (1+state.zabs)
     oplot, [wvlim_plt[0],wvlim_plt[0]], tmpy,color=clr.cyan, thick=2, $
            linestyle=0
     oplot, wvlim_plt, [tmpy[1],tmpy[1]],color=clr.cyan, linestyle=0, thick=2
     oplot, [wvlim_plt[1],wvlim_plt[1]], tmpy,color=clr.cyan, thick=2, $
            linestyle=0
  endfor 
  
  ;; Plot start/end of detection based on z of QSO
  oplot, ([1.,1.]+state.zqso)*(1.+state.dvqso*cinv)*$
         qalstr[state.curqal].wrest[idx_dblt+1]/(1+state.zabs),$ ; wvII  to maximize
         state.xymnx[[1,3]],color=clr.red,linestyle=1,thick=3
  oplot, state.minloswav/([1.,1.]+state.zabs),$ ; combination of wvI and dvgal
         state.xymnx[[1,3]],color=clr.green,linestyle=1,thick=3
  
  ;; Plot (at the bottom) the positions of lines if alternately
  ;; identified
  for aa=0,2 do begin
     tmpy = (1+aa)*0.03125*dy+state.xymnx[1] ; lower position
     case aa of
        0: begin
           if keyword_set(state.dbltchk1) then begin
              tmpdblt = state.dbltchk1
              tmpclr = clr.orange
              tmpx = 0.32*dx    ; tmpy = [-0.15,-0.05] 
           endif else tmpx = -1 ; exit flag
        end
        1: begin
           if keyword_set(state.dbltchk2) then begin
              tmpdblt = state.dbltchk2
              tmpclr = clr.darkgreen
              tmpx = 0.28*dx    ; tmpy = [-0.12,-0.02]
           endif else tmpx = -1 
        end
        2: begin
           if keyword_set(state.dbltchk3) then begin
              tmpdblt = state.dbltchk3
              tmpclr = clr.purple
              tmpx = 0.24*dx   ; tmpy = [-0.09,0.01]
           endif else tmpx = -1
        end
     endcase
     if tmpx eq -1 then continue
     
     zalt = [qalstr[state.curqal].wrest[idx_dblt]/tmpdblt.wvII,$ ; (1 + zalt)/(1+zabs)
             qalstr[state.curqal].wrest[idx_dblt]/tmpdblt.wvI]
     
     oplot, [1,1]*tmpdblt.wvI*zalt[0], [tmpy,tmpy+0.0625*dy], color=tmpclr, $
            thick=2, linestyle=0
     oplot, [1,1]*tmpdblt.wvI*zalt[1], [tmpy,tmpy+0.0625*dy], color=tmpclr, $
            thick=2, linestyle=0
     oplot, [1,1]*tmpdblt.wvII*zalt[1], [tmpy,tmpy+0.0625*dy], color=tmpclr, $
            thick=2, linestyle=0
     oplot, [tmpdblt.wvI*zalt[0],tmpdblt.wvII*zalt[1]], $
            replicate(tmpy,2), color=tmpclr, thick=2
     xyouts, qalstr[state.curqal].wrest[1]+tmpx, tmpy, tmpdblt.ion, $
             color=tmpclr, alignment=1, charsize=1.5

  endfor                        ; loop aa=0,2
  

  ;; Indicate centroids; dashed lines
  oplot, qalstr[state.curqal].wrest[idx_dblt]*[1.,1], [-9,9], $
         color=clr.blue, linestyle=2, thick=1.2
  oplot, qalstr[state.curqal].wrest[idx_dblt+1]*[1.,1], [-9,9], $
         color=clr.blue, linestyle=2, thick=1.2
  oplot, wv[gd]/(1.+state.zabs), wv[gd]*0. + 1., color=clr.green, linestyle=3

end ; sdss_chkciv_metals


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; Set BAL flag
pro sdss_chkciv_chgbalflg, state
  common sdss_chkciv_cmm
  
  mtch = where(strtrim(qalstr.qso_name,2) eq $
               strtrim(qalstr[state.curqal].qso_name,2),nmtch)
  if nmtch eq 0 then $
     stop,'sdss_chkciv_chgbalflg: No matches for QSO'

  if civstr[state.curqal].balflg eq 0 then begin
     ;; Revert to old version
     if qalstr[state.curqal].balflg eq 0 then $
        nwflg = sdss_getbalflg(/visual) $
     else nwflg = qalstr[state.curqal].balflg
  endif else begin
     nwflg = 0
  endelse 
  if keyword_set(state.debug) then $
     print,'sdss_chkciv_chgbalflg debug: setting QSO BALFLG = '+$
           strtrim(nwflg,2)+' for all instances, '+strtrim(nmtch,2)

  civstr[mtch].balflg = nwflg
end                             ; sdss_chkciv_chgbalflg


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; Set BLEND flag (rating[9])
pro sdss_chkciv_blend, state
  common sdss_chkciv_cmm
  ;; So typically rating[9] is equal to -1 when entered this field but
  ;; to indicate that the blending has been considered, will toggle
  ;; between 1 (blended) and 0 (not blended)
  ;; now rating[9] = sdss_getblendflg(/self) (or 6) is from
  ;; sdss_civsearch_broad() 
  if (civstr[state.curqal].rating[9] and sdss_getblendflg()) eq 0 then $
     nwflg = sdss_setblendflg(civstr[state.curqal].rating[9]) $ ; binary or
  else nwflg = sdss_getblendflg(/unblend)                       ; erase even self-belnd     

  if keyword_set(state.debug) then $
     print,'sdss_chkciv_blend debug: setting QSO.rating[9] = '+$
           strtrim(nwflg,2)+' for current # '+strtrim(state.curqal,2)

  civstr[state.curqal].rating[9] = nwflg
end                             ; sdss_chkciv_blend


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;
; UPDATE LYA+METALS
pro sdss_chkciv_update, state
  sdss_chkciv_updinfo, state
  sdss_chkciv_metals, state
  return
end                             ; sdss_chkciv_update


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; Setup data
pro sdss_chkciv_setup, state
  common sdss_chkciv_cmm
  
  if keyword_set(state.debug) then $
     print,'sdss_chkciv_setup debug: state.curqal = ',strtrim(state.curqal,2)
  if keyword_set(state.multi) then $
     sdss_chkciv_llist, state   ; changing line list

  ;; Shouldn't have to check if state.curqal ge state.nqal
  ;; because sdss_chkciv_next should do that

  ;; Read data
  if not keyword_set(state.flg_snglspec) then begin
     state.filnam = strtrim(qalstr[state.curqal].sdss_obs[0],2)
     if state.filnam ne state.sv_filnam then begin
        ;; New file
        if keyword_set(state.debug) then $
           print, 'sdss_chkciv_setup debug: Reading new spectrum ', $
                  qalstr[state.curqal].sdss_obs[0]
        state.sv_filnam = state.filnam
        
        ;; Spectrum passed around in common block
        parse_sdss, sdssdir+state.filnam, fx, wv, head=hd, sig=sig0, npix=npix 
        ;; Set common QSO info
        state.ra = sxpar(hd,'RA')
        state.dec = sxpar(hd,'DEC')
        mag = sxpar(hd,'MAG')   ; ugriz
        prs = strsplit(mag,/extract,count=nprs)
        state.mag = float(prs[2])
        
        ;; Figure out the other candidates in sightline
        loscand = where(strtrim(qalstr[state.curqal].qso_name,2) eq $
                        strtrim(qalstr.qso_name,2),nloscand)
        
        ;; Figure out other metal-line systems in sightline
        if keyword_set(mlsstr) then $
           losmls = where(strtrim(qalstr[state.curqal].qso_name,2) eq $
                       strtrim(mlsstr.qso_name,2),nlosmls) $
        else nlosmls = 0
        
        ;; Retrieve continuum which may not exist
        test = file_search(sdssdir+strtrim(qalstr[state.curqal].abslin_fil,2)+'*',count=ntest)
        if ntest gt 0 then begin
           cstrct = xmrdfits(sdssdir+strtrim(qalstr[state.curqal].abslin_fil,2),1,/silent)
           cindx = fix(alog(qalstr[state.curqal].cflg)/alog(2))
           
           conti = cstrct.conti[0:cstrct.npix-1,*] 
           fxnorm = fx* 0.
           bd = where(conti[*,cindx] eq 0. or finite(conti[*,cindx],/nan),complement=gd)
           fxnorm[gd] = fx[gd]/conti[gd,cindx]
           if bd[0] ne -1 then fxnorm[bd] = 0.
           ;; Set new error including continuum
           bdpix = where(sig0 eq 0.)
           ;; sig = sqrt(sig0^2*conti[*,cindx]^2 + $
           ;;            cstrct.sigconti[0:cstrct.npix-1,cindx]^2*fx^2)/$
           ;;       conti[*,cindx] ; unnorm-sigma + conti error
           ;; if bdpix[0] ne -1 then sig[bdpix] = 0. ; handled in
           ;; function below
           signorm = sdss_calcnormerr(fx,sig0,cstrct,cflg=qalstr[state.curqal].cflg)
           sig = sdss_calcnormerr(fx,sig0,cstrct,cflg=qalstr[state.curqal].cflg,/unnorm)
           if keyword_set(debug) then begin
              if qalstr[state.curqal].cflg eq sdss_getcflg(/eig) then $
                 print,'sdss_chkciv_setup debug: using eigen continuum'
              if qalstr[state.curqal].cflg eq sdss_getcflg(/hyb) then $
                 print,'sdss_chkciv_setup debug: using hybrid continuum'
           if qalstr[state.curqal].cflg eq sdss_getcflg(/spl) then $
              print,'sdss_chkciv_setup debug: using spline continuum'
        endif 
           
           ;; Set contrived "spectrum" to plot all detected lines, making a
           ;; square waveform
           sdss_pltabslin, fx, cstrct, /noplt, wave=wv, fx_lin=fx_lin, $
                           wv_lin=wv_lin
        endif else begin
           ;; abslin file does not exist so just make dummy files
           print,'sdss_chkciv_setup: abslin file DNE.'
           conti = replicate(0.,npix,3)
           wv_lin = [0.,0.]
           fx_lin = [0.,0.]
        endelse 
        
        ;; Other different info
        state.zqso = qalstr[state.curqal].z_qso
     endif                      ; different sdss_obs[0]
  endif

  ;; Things that do change
  ;; Set zabs
  if (civstr[state.curqal].zabs_final[idx_dblt] NE 0) then begin
     state.zabs = civstr[state.curqal].zabs_final[idx_dblt]
  endif else begin
     state.zabs = qalstr[state.curqal].zabs_orig[idx_dblt]
  endelse

  tmpdblt = dblt_retrieve(qalstr[state.curqal].wrest[idx_dblt])
  case tmpdblt.ion of
     'CIV': minloswav = 1310.
     'FkIV': minloswav = 1310.
     else: minloswav = 1230.
  endcase
  state.minloswav = minloswav*(1.+state.zqso) > $
                    qalstr[state.curqal].wrest[idx_dblt]*(1.+state.dvgal*cinv)
  
  ;; Set EW and dvabs
  if (civstr[state.curqal].ew_final[idx_dblt] NE 0) then begin
     sdss_chkciv_setewstr, state, civstr[state.curqal].ew_final, $
                           civstr[state.curqal].sigew_final
     state.dvabs = $
        (civstr[state.curqal].zabs_final[idx_dblt+1]-civstr[state.curqal].zabs_final[idx_dblt]) $
        /(1.+civstr[state.curqal].zabs_final[idx_dblt]) * 3.e5 ; km/s
  endif else begin
     sdss_chkciv_setewstr, state, qalstr[state.curqal].ew_orig, $
                           qalstr[state.curqal].sigew_orig
     state.dvabs = $
        (civstr[state.curqal].zabs_orig[idx_dblt+1]-civstr[state.curqal].zabs_orig[idx_dblt]) $
        /(1.+civstr[state.curqal].zabs_orig[idx_dblt]) * 3.e5 ; km/s
  endelse 

  ;; Set Rating
  state.rtg = civstr[state.curqal].rating[0]

  ;; xymnx
  state.xymnx[0] = 1215.6701*(1.+state.zabs) - 200.
  state.xymnx[2] = 1215.6701*(1.+state.zabs) + 200.
  gd = where(wv GT state.xymnx[0] AND wv LT state.xymnx[2], ngd)
  
  if ngd GT 1 then begin
     srt = sort(fx[gd])
     ymd = fx[gd[srt[round(0.9*ngd)<(ngd-1)]]]
  endif else ymd = 0.
  state.xymnx[1] = -1.
  state.xymnx[3] = ymd*1.5
  state.svxymnx = state.xymnx

  ;; Fit
  state.civ_conti = ymd
  fit = replicate(state.civ_conti, npix)

  ;; Set flg
  state.velplt[*].flg = 0
  gd = where(state.velplt.wrest*(state.zabs+1) GT min(wv) AND $
             state.velplt.wrest*(state.zabs+1) GT (state.zqso+1.)*1215.6701 AND $
             state.velplt.wrest*(state.zabs+1) LT max(wv), na)
  if na EQ 0 then $
     print,'sdss_chkciv_setup: no other lines to plot' $
  else begin
     state.nplt = na
     state.velplt[gd].flg = 1
     
     ;; Just do the< plots
     state.all_velo[*,gd] = x_allvelo(wv, state.zabs, $
                                      state.velplt[gd].wrest,$
                                      state.vmnx, $
                                      all_pmnx=all_pmnx, NPIX=5000L)
     state.all_pmnx[*,gd] = all_pmnx
  endelse 

  return
end                             ; sdss_chkciv_setup


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set Lines
pro sdss_chkciv_llist, state

  common sdss_chkciv_cmm
  tmpdblt = dblt_retrieve(qalstr[state.curqal].wrest[idx_dblt])
  case tmpdblt.ion of 
     'CIV': llist = getenv('XIDL_DIR')+'/SDSS/CIV/sdss_civ.lst'
     'FkIV': llist = getenv('XIDL_DIR')+'/SDSS/CIV/sdss_fkiv.lst'
     'SiIV': llist = getenv('XIDL_DIR')+'/SDSS/SiIV/sdss_siiv.lst'
     'FuIV': llist = getenv('XIDL_DIR')+'/SDSS/SiIV/sdss_fuiv.lst'
     'MgII': llist = getenv('XIDL_DIR')+'/SDSS/MgII/sdss_mgii.lst'
     'FkII': llist = getenv('XIDL_DIR')+'/SDSS/MgII/sdss_fkii.lst'
     'CaII': llist = getenv('XIDL_DIR')+'/SDSS/CaII/sdss_caii.lst'
     else: begin
        print,'sdss_chkciv_llist: state.dblt.ion DNE ',state.dblt.ion
        llist = getenv('XIDL_DIR')+'/SDSS/CIV/sdss_civ.lst' ; default
     endelse 
  endcase 
  lines = x_setllst(llist, 0)
  state.ntrans = n_elements(lines)
  state.velplt[0:state.ntrans-1].wrest = lines.wave
  state.velplt[0:state.ntrans-1].name = lines.name
  delvarx, lines
  state.velplt[0:state.ntrans-1].ymnx = [-0.11, 1.39]

  weak = where(abs(state.velplt.wrest - 1808.0130d) LT 0.01 OR $
               abs(state.velplt.wrest - 1611.2005d) LT 0.01 OR $
               abs(state.velplt.wrest - 2026.136d) LT 0.01 OR $
               abs(state.velplt.wrest - 2260.7805d) LT 0.01, nwk)
  
  if nwk NE 0 then state.velplt[weak].ymnx = [0.7, 1.1]
               
end                             ; sdss_chkciv_llist

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Move
pro sdss_chkciv_next, state, KILL=kill, skip=skip,  back=back
  common sdss_chkciv_cmm

  if keyword_set(back) then begin
     if state.curqal eq 0 then $
        print,'sdss_chkciv_next, /back: On first QSO, cannot go further back' $
     else state.curqal = state.curqal - 1
     
     if keyword_set(state.debug) then $
        print, "sdss_chkciv_next, /back debug: state.curqal=", strtrim(state.curqal,2), $
               "; state.nqal=", strtrim(state.nqal,2)
     ;; Shouldn't have to do a kill check or a save
     return
  endif 

  if keyword_set(skip) then begin
     ;; Skip to next unrated (in the future)
     tmp = lindgen(state.nqal)
     gd = where(civstr.rating[0] eq sdss_getrating(/unrated) and $
               tmp gt state.curqal)
     if gd[0] eq -1 then begin
        ;; See if any unrated at all
        gd = where(civstr.rating[0] eq sdss_getrating(/unrated))
        if gd[0] eq -1 then state.curqal = state.nqal-1 $ ; all done!
        else state.curqal = gd[0] ; in the past
     endif else state.curqal = gd[0] ; in the future
  endif else state.curqal = state.curqal + 1 ; all done!

  if keyword_set(state.debug) then $
     print, "sdss_chkciv_next debug: state.curqal=", strtrim(state.curqal,2), $
            "; state.nqal=", strtrim(state.nqal,2)

  if state.curqal ge state.nqal then begin 
     kill = 1
     print,'sdss_chkciv_next: End of candidates; you should quit (auto-saves).'

     if keyword_set(state.debug) then $
        print,'sdss_chkciv_next debug: kill set'
  endif else begin
     kill = 0
     ;; Now that state.curqal settled, can assess whether to change
     ;; state.dblt
     if keyword_set(state.multi) then begin
        state.dblt = dblt_retrieve(qalstr[state.curqal].wrest[idx_dblt])
        if keyword_set(state.debug) then $
           print,'sdss_chkciv_next debug: doublet now ',state.dblt.ion
     endif
  endelse
  
  ;; Always save
  if state.sav_ech EQ 1 or (state.curqal mod 100) eq 0 then $
        sdss_chkciv_svciv, state 

  return
end                             ; sdss_chkciv_next


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Next
pro sdss_chkciv_updinfo, state
  common sdss_chkciv_cmm

  ;; Name
  widget_control, state.name_id, $
    set_value=strtrim(qalstr[state.curqal].qso_name,2)

  ;; RA, DEC
  widget_control, state.mag_id, set_value=state.mag
  widget_control, state.ra_id, set_value=state.ra
  widget_control, state.dec_id, set_value=state.dec
  widget_control, state.rtg_id, set_value=state.rtg
  widget_control, state.num_id, set_value=state.curqal  
  widget_control, state.ew_id, set_value=state.ew
  widget_control, state.ewrto_id, set_value=state.ewrto
  widget_control, state.dvabs_id, set_value=state.dvabs

  ;; zabs
  widget_control, state.zabs_id, set_value=state.zabs
  widget_control, state.zqso_id, set_value=state.zqso
  widget_control, state.balflg_id, set_value=civstr[state.curqal].balflg ; not qalstr b/c may have been updated
  widget_control, state.blend_id, set_value=civstr[state.curqal].rating[9]
  widget_control, state.note_id, set_value=civstr[state.curqal].notes
  return

end                             ; sdss_chkciv_updateinfo


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Set CIV structure
pro sdss_chkciv_setciv, state, DEFINITE=definite, GOOD=good, MAYBE=maybe, BAD=bad

  common sdss_chkciv_cmm

  flg_civ = sdss_getrating(good=good,definite=definite, maybe=maybe, bad=bad)

   civstr[state.curqal].rating[0] = flg_civ
 
  if keyword_set(state.debug) then $
     print,'sdss_chkciv_setciv debug: state.curqal=',strtrim(state.curqal,2)
  return

end                             ; sdss_chkciv_setciv

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chkciv_svciv, state

  common sdss_chkciv_cmm

  ;; Save structure
  mwrfits, civstr, state.civfil, /create, /silent
  if not keyword_set(state.sav_ech) then $
  print,'sdss_chkciv_svciv: saved ',state.civfil,string(state.curqal,format='(1x,i5)')
  spawn,'gzip -f '+state.civfil

  if keyword_set(state.debug) then begin
     gd = where(civstr.rating[0] eq sdss_getrating(/unrated),ngd)
     gd2 = where(gd lt state.curqal,ngd2) ; current not rated
     print, 'sdss_chkciv_svciv debug: number left to rate', ngd
     print, 'sdss_chkciv_svciv debug: number skipped',ngd2
  endif 

  if keyword_set(state.debug) then print,''
  return
end                             ; sdss_chkciv_svciv

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; MAIN PROGRAM ;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sdss_chkciv, qalfil, civfil, dblt_name=dblt_name, multi=multi, $
                 hidealt=hidealt, altdblt_name=altdblt_name, xsize=xsize, $
                 ysize=ysize, dvqso=dvqso, debug=debug, sav_ech=sav_ech, $
                 label=label, initials=initials, svrtg=svrtg, mlsfil=mlsfil, $
                 inspec=inspec, _extra=extra

  if  N_params() LT 2  then begin 
     print,'Syntax - sdss_chkciv, qalfil, civfil, [dblt_name=, /multi, /hidealt, /debug, '
     print,'                  altdblt_name=, dvqso=, /debug, /sav_ech, xsize=, ysize=, /svrtg,'
     print,'                  /label, initials=, /svrtg, mlsfil=, inspec=, _extra=]'
     return
  endif 

  ;;  Optional Keywords
  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then ysize = ssz[1]-100
  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if not keyword_set(altdblt_name) then begin
     if dblt_name eq 'CIV' then altdblt_name = 'MgII' $
     else altdblt_name = 'CIV'
  endif 
  if not keyword_set(dvqso) then dvqso = -3000. ; km/s matches sdss_fndciv
  if not keyword_set(dvgal) then dvgal = 5000.  ; km/s matches sdss_fndciv
  if keyword_set(initials) then $
     ninitials = n_elements(initials) $
  else begin
     initials = 0
     ninitials = 0
  endelse 
  if not keyword_set(mlsfil) then mlsfil = ''

  ;; Allow qalfil to be structure for viewing
  if size(qalfil,/type) eq 8 then begin
     if not keyword_set(civfil) then civfil = qalfil
     ;; Going to have to make qalfil and civfil strings becasue
     ;; everything is set up for that but will return structures
     qalfil_strct = qalfil
     qalfil = 'sdss_chkciv_input.fit'
     mwrfits,qalfil_strct,qalfil,/create,/silent
     if size(civfil,/type) ne 7 then begin
        civfil_strct = civfil
        civfil = 'sdss_chkciv_output.fit'
        test = file_search(civfil+'*',count=ntest)
        if ntest ne 0 then $
           spawn,'\rm '+test
     endif else civfil_strct = 0
     svrtg = 1
  endif else begin
     qalfil_strct = 0
     civfil_strct = 0
  endelse 

  ;; Set structures for state structure
  tmp = { velpltstrct }
  dblt = dblt_retrieve(dblt_name)
  dbltalt = dblt_retrieve(altdblt_name)
  dbltchk1 = 0
  dbltchk2 = 0
  dbltchk3 = 0
  if not keyword_set(hidealt) then begin
     case dblt.ion of
        'CIV': begin
           ;; Want to indicate other redshifts of lines
           ;; Order from 1 = widest to 3 = narrowest
           dbltchk1 = dblt_retrieve('CaII')
           dbltchk2 = dblt_retrieve('SiIV')
           dbltchk3 = dblt_retrieve('FeII')
        end
        'FkIV': begin; Fake doublet
           dbltchk1 = dblt_retrieve('SiIV')
           dbltchk2 = dblt_retrieve('FeII')
           dbltchk3 = dblt_retrieve('CIV')
        end
        'SiIV': begin
           dbltchk1 = dblt_retrieve('CaII')
           dbltchk2 = dblt_retrieve('FeII')
        end
        'CaII': begin
           dbltchk1 = dblt_retrieve('SiIV')
           dbltchk2 = dblt_retrieve('FeII')
        end
        'MgII': begin
           dbltchk1 = dblt_retrieve('SiIV')
           dbltchk2 = dblt_retrieve('FeII')
           dbltchk3 = dblt_retrieve('CIV')
        end
        'FkII': begin ; Fake Doublet
           dbltchk1 = dblt_retrieve('SiIV')
           dbltchk2 = dblt_retrieve('FeII')
           dbltchk3 = dblt_retrieve('CIV')
        end
        else: begin
        end
     endcase
  endif 

  ;; STATE
  state = {             $
          nqal: 0L, $
          curqal: 0L, $         ; tracks position in input
          qalfil: strtrim(qalfil,2), $
          civfil: strtrim(civfil,2), $
          mlsfil: strtrim(mlsfil,2), $     ; could be structure
          save_ratings:keyword_set(svrtg), $
          zabs: 0., $           ; could use qalstr[state.curqal].zabs[idx_dblt] everywhere and eliminated this...
          dvabs: 0., $
          ew: '', $             ; value +/- error
          ewrto: '', $
          rtg: 0, $             ; Rating
          zqso: 0., $           ; QSO info
          dvqso:dvqso, $           
          dvgal:dvgal, $           
          minloswav:0., $ ;  based on doublet
          mag:0., $             ; ... from header
          ra:0., $              ; 
          dec:0., $             ; 
          filnam:'', $          ; Current spectrum
          flg_snglspec:keyword_set(inspec), $ ; lcok spec after sdss_chkciv_icmmn
          sv_filnam:'', $       ; Don't read in if same
          ntrans: 0L, $         ; PLOTTING LINES
          vmnx: [-800., 800.], $
          nplt: 0, $
          civ_conti: 0., $
          all_velo: fltarr(5000, 300, /nozero), $  
          all_pmnx: lonarr(3, 300, /nozero), $  
          velplt: replicate(tmp, 300), $
          dblt: dblt, $         ; from dblt_retrieve()
          dbltalt: dbltalt, $
          dbltchk1: dbltchk1, $
          dbltchk2: dbltchk2, $
          dbltchk3: dbltchk3, $
          multi:keyword_set(multi), $ ; 0: fix to state.dblt; 1: float to qalstr[state.curqal].wrest[idx_dblt]
          psfile: 0, $
          help: strarr(50), $
          debug: keyword_set(debug), $
          label: keyword_set(label), $
          initials: initials, $
          ninitials: ninitials, $
          sav_ech: keyword_set(sav_ech), $
          svxymnx: fltarr(4), $
          xymnx: fltarr(4), $
          base_id: 0L, $        ; Widgets
          mdraw_id: 0L, $       ; Spec Window
          mdrawbase_id: 0L, $
          zabs_id: 0L, $
          zqso_id: 0L, $
          blend_id: 0L, $
          balflg_id: 0L, $
          name_id: 0L, $
          top_id: 0L, $
          tdraw_id: 0L, $
          bottom_id: 0L, $
          info_id: 0L, $
          ra_id: 0L, $
          dec_id: 0L, $
          mag_id: 0L, $
          rtg_id: 0L, $          ; display the rating
          num_id:0L, $
          ew_id: 0L, $
          ewrto_id: 0L, $
          ntot_id:0L, $
          dvabs_id:0L, $
          note_id:0L, $
          xcurs: 0L, $
          ycurs: 0L, $
          pos: [0.1,0.1,0.95,1.0], $
          pltmnx: fltarr(4), $
          size: fltarr(2) $
          }

  ;;;;;;;;;;;;;;
  ;; Initialize common block
  sdss_chkciv_icmmn, state, fail=fail, inspec=inspec, _extra=extra ; idblt=
  if keyword_set(fail) then return


  ;;;;;;;;;;;;;;
  ;; SETUP LINES
  sdss_chkciv_llist, state

  ;; Other setup

  ;; WIDGET
  base = WIDGET_BASE( title = 'sdss_chkciv: Rate candidates', /column, $
                      UNAME='BASE', /tlb_size_events,xoffset=200L)
  state.base_id = base
  

  state.top_id = WIDGET_BASE( state.base_id, /column, $
                              /base_align_center,/align_center, $
                              xsize=xsize, ysize=round(1*ysize/3.), $
                              uvalue='TOP_BASE', frame=2)
  ;; Here's where it's set that TEXT only works in top window
  state.size  = [xsize,ysize]
  state.tdraw_id = WIDGET_DRAW(state.top_id, xsize=xsize, $
                               ysize=round(ysize/3), /frame, retain=2, $
                               uvalue='TEXT', /keyboard_events)
  state.bottom_id = WIDGET_BASE( state.base_id, /row, $
                                 /base_align_center,/align_center, $
                                 xsize=xsize, ysize=round(2*ysize/3.), $
                                 uvalue='TOP_BASE', frame=2)

  ;;;;;; Info window ;;;;;;;;;;;
  state.info_id = $
     WIDGET_BASE( state.bottom_id, /column, /base_align_center,/align_center, $
                  uvalue='INFO_BASE', frame=2, xsize=xsize/3.)
  ;; QSO Name, RA, Dec
  qsoinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.name_id = cw_field(qsoinf, title='QSO:', value=' ', xsize=15)
  state.ra_id = cw_field(qsoinf, title='RA:', value=0., xsize=7)
  state.dec_id = cw_field(qsoinf, title='DEC:', value=0., xsize=7)

  ;; RMAG, RA, DEC
  qsoinf2 = widget_base(state.info_id, /row, /align_center, frame=2)
  state.mag_id = cw_field(qsoinf2, title='MAG:', value=0., xsize=5)
  state.zqso_id = cw_field(qsoinf2, title='zqso:',value = state.zqso, xsize=7)
  state.ntot_id = cw_field(qsoinf2, title='Total #:',value = state.nqal, xsize=7)
  
  ;; Absorber info and Current #
  civinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.zabs_id = cw_field(civinf, title='zabs:', value=state.zabs, xsize=7)
  state.dvabs_id = cw_field(civinf, title='dvabs:', value=state.dvabs, xsize=6)
  state.num_id = cw_field(civinf, title='Cur #', value=0,xsize=7,$
                          /return_events,uvalue='CURNUM')

  ;; Equivalent width and ratio information
  ewinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.ew_id = cw_field(ewinf, title='EW:', value=state.ew, xsize=12)
  state.ewrto_id = cw_field(ewinf, title='EW ratio:', value=state.ewrto, xsize=12)

  ;; Rating, BAL flg and button 
  rtginf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.rtg_id = cw_field(rtginf, title='Rating:', value=sdss_getrating(/unrated),$
                          xsize=3)
  state.balflg_id = cw_field(rtginf, title='BALFLG', value=0,xsize=1)
  balflg = WIDGET_BUTTON(rtginf, value='Chg BAL (f)', uvalue='CHGBAL')
  

  ;; BUTTONS
  ;; Ratings
  butbase = widget_base(state.info_id, /row, /align_center, frame=2)
  bad  = WIDGET_BUTTON(butbase, value='BAD ('+strtrim(sdss_getrating(/bad),2)+')',$
                       uvalue='BAD')
  maybe = WIDGET_BUTTON(butbase, value='MAYBE ('+$
                        strtrim(sdss_getrating(/maybe),2)+')',uvalue='MAYBE')
  good = WIDGET_BUTTON(butbase, value='GOOD ('+$
                       strtrim(sdss_getrating(/good),2)+')',uvalue='GOOD')
  defi = WIDGET_BUTTON(butbase, value='DEFINITE ('+$
                       strtrim(sdss_getrating(/def),2)+')',uvalue='DEFINITE')
  ;; Niceties: Flag Blend, splt, alternate ID, Save
  butbase2 = widget_base(state.info_id, /row, /align_center, frame=2)
  state.blend_id = cw_field(butbase2, title='BLEND', value=0,xsize=1)
  blend = WIDGET_BUTTON(butbase2, value='BLEND (d)',uvalue='BLEND') 
  splt = WIDGET_BUTTON(butbase2, value='PLOT ( )',uvalue='SPLT') 
  mgii = WIDGET_BUTTON(butbase2, value=state.dbltalt.ion+' (a)', $
                       uvalue=state.dbltalt.ion)
  sav = WIDGET_BUTTON(butbase2, value='SAVE (w)', uvalue='SAVE')

  ;; Movement ones
  butbase3 = widget_base(state.info_id, /row, /align_center, frame=2)
  back = WIDGET_BUTTON(butbase3, value='BACK (b)', uvalue='BACK')
  next = WIDGET_BUTTON(butbase3, value='NEXT (n)', uvalue='NEXT')
  skip = WIDGET_BUTTON(butbase3, value='SKIP (s)', uvalue='SKIP')
  done = WIDGET_BUTTON(butbase3, value='DONE (q)',uvalue='DONE')

  ;; Notes
  noteinf = widget_base(state.info_id, /row, /align_center, frame=2)
  state.note_id = cw_field(noteinf, title='Notes (remember to hit Return):', $
                          /column, value='', xsize=xsize/20.,$
                          /return_events,uvalue='NOTES')

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;      Metals DRAW
  state.mdrawbase_id = $
     WIDGET_BASE( state.bottom_id, /row, /base_align_center,/align_center, $
                  uvalue='SDRAW_BASE', frame=2)

  state.mdraw_id = widget_draw(state.mdrawbase_id, xsize=round(xsize*2./3), $
                               ysize=round(2.*ysize/3), /frame, retain=2, $
                               uvalue='SDRAW')

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;      Define Help (make sure matches 'TEXT' in sdss_chkciv_ev and
  ;;      buttons above)
  state.help[0] = '       Help Menu      '
  state.help[1] = 'Keys only work when top window active.'
  state.help[2] = '3 -- Definite'
  state.help[3] = '2 -- Good'
  state.help[4] = '1 -- Maybe'
  state.help[5] = '0 or ` -- Bad'
  state.help[6] = 'f -- Toggle BAL flag for QSO'
  state.help[7] = 'b -- Back'
  state.help[8] = 'n -- Next'
  state.help[9] = 's -- Skip'
  state.help[10] = 'a -- '+state.dbltalt.ion+' overplot'
  state.help[11] = '  -- x_specplot with continuum and detected lines'  ; spacebar
  state.help[12] = 't -- x_specplot with both continua'  
  state.help[13] = '- -- sdss_pltsnr of spectrum' ; dash
  state.help[14] = 'i -- Print LOS candidate summary.'
  state.help[15] = 'E -- write "ghost emission" to notes.'
  state.help[16] = 'B -- write "bad conti" to notes.'
  state.help[17] = 'A -- write "AlII" to notes.'
  state.help[17] = 'C -- write "CIV" to notes.'
  state.help[18] = 'D -- write "DLA" to notes.'
  state.help[19] = 'F -- write "FeII" to notes.'
  state.help[20] = 'K -- write "sky lines" to notes.'
  state.help[21] = 'M -- write "MgII" to notes.'
  state.help[22] = 'S -- write "SiII" to notes.'
  state.help[23] = 'V -- write "SiIV" to notes.'
  state.help[24] = 'l or x -- Set left  bound of 1548 EW region.'
  state.help[25] = 'r -- Set right bound of 1548 EW region.'
  state.help[26] = 'L or X -- Set left  bound of 1550 EW region.'
  state.help[27] = 'R -- Set right bound of 1550 EW region.'
  state.help[28] = 'u -- Undo changed bounds of 1548 EW region.'
  state.help[29] = 'U -- Undo changed bounds of 1550 EW region.'
  state.help[30] = 'w -- Save (write)'
  state.help[31] = 'q -- Quit (saves)'
  state.help[33] = "Enter number in 'Cur #' field to jump."
  state.help[34] = 'Green/red brackets -- Wavelength bounds for measurements.'
  state.help[35] = 'Blue bracket -- Other candidates in LOS.'
  state.help[36] = 'Vertical, blue lines in x_spcplot -- All detected lines.'
  state.help[37] = 'Bottom brackets indicate positions of alternate (wide) doublets.'
  state.help[38] = 'Use Comment field to write to text file; remember to Return; no colons.'
  state.help[39] = 'Vertical, green/red dotted line indicates start/end of search region.'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Realize
  WIDGET_CONTROL, base, /realize

  ;; Load data
  sdss_chkciv_setup, state

  ;; PLOT
  sdss_chkciv_update, state

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  ;; Send to the xmanager
  xmanager, 'sdss_chkciv', base

  ;; Print final results
  sdss_printratedsumm,civfil

  ;; Handle if input were structures
  if keyword_set(qalfil_strct) then begin
     tmp = xmrdfits(qalfil,1,/silent) 
     spawn,'\rm '+qalfil+'*'
     qalfil = tmp
  endif 
  if keyword_set(civfil_strct) then begin
     spawn,'\rm '+civfil+'*'
     civfil = civstr
  endif 

  return
end
	
