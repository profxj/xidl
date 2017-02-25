;+ 
; NAME:
; mike_redux
;    Version 1.0
;
; PURPOSE:
;   Calls various routines to analyse the data
;
; CALLING SEQUENCE:
;   
;   mike_redux, wfccd, maskid, expsr, XSIZE=, YSIZE=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_redux, mike
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-May-2005 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
; Events
;;;;

pro mike_redux_event, ev

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'SHOWQA' : begin
          if state.qa_sel NE -1 then begin
              spawn, 'gv '+state.qa_fil[state.qa_sel]+'.gz'
          endif
      end
      'QALIST' : state.qa_sel = ev.index
      'DOSLITF' : mike_redux_doslitf, state
      'DOGAIN' : mike_redux_dogain, state
      'DOFLATS' : mike_redux_doflat, state
      'DOARCS' : mike_redux_doarcs, state
      'DOOBJ' : mike_redux_doobj, state
      'OBJLIST' : begin
          state.obj = ev.index + 1
          mike_redux_updset, state
          mike_redux_updobj, state, ev.index
      end
      'SETUPLIST' : mike_redux_changeset, state, ev.index
      'SETUPBUT' : begin
          if ev.select EQ 1 then state.flg_setup = 1 else $
            state.flg_setup = 0
          mike_redux_init, state
      end
      'OBJOPTION' : begin
          widget_control, state.objopt_id, get_value=tmp
          state.flg_objopt[0:n_elements(tmp)-1] = tmp
          state.clobber = tmp[3]
      end
      'SIDES' : begin
          widget_control, state.sides_id, get_value=tmp
          state.side = tmp
      end
      'EDITSTRCT' : begin
          tmp = state.struct
          mike_editstrct, tmp
          state.struct = tmp
      end
      'TOSSOUT' : begin
          mike_redux_tossout, state
          if flg NE 0 then begin
              mike_redux_updset, state
              mike_redux_updobj, state, ev.index
          endif
      end
      'WRITE' : begin
          mike_redux_write, state
      end
      'SPECPLOT' : mike_specplt
      'DONE' : begin
          mike_redux_write, state
          widget_control, ev.top, /destroy
          return
      end
      'DNSV' : begin
          widget_control, ev.top, /destroy
          return
      end
      else :
  endcase

  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_write, state

  print, 'mike_redux:  Writing to ', state.outfil
  mike_wrstrct, state.struct, FITS=state.outfil

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_tossout, state, flg

  flg = 0

  widget_control, state.arcs_id, get_value=arcs
  widget_control, state.flats_id, get_value=flats
  stop

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_dogain, state

  ;; Side
  side = (where(state.side,nside) + 1)

  ;; Run
  tmp = state.struct
  mike_setgain, tmp, state.setval, side
  state.struct = temporary(tmp)

  print, 'mike_redux:  It is best to WRITE out the structure now to save '
  print, ' the updated gain.'
  
  ;; Update setup
  mike_redux_updset, state

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_doflat, state

  ;; Side
  side = (where(state.side,nside) + 1)
  ;; Gain
  rdxi = where(state.rdx_steps EQ 'Gain')
  sum =0L
  for qq=0L,nside-1 do begin
      sum = sum + state.rdx_sel[side[qq]-1,rdxi]
  endfor
  flggain = (sum EQ nside)
  
  ;; Run
  tmp = state.struct
  if keyword_set(state.nogain) OR (flggain) then nogain = 1 else nogain = 0
  mike_allflat, tmp, state.setval, side, CLOBBER=state.clobber, NOGAIN=nogain
  state.struct = temporary(tmp)

  print, 'mike_redux:  It is best to WRITE out the structure now to save '
  print, ' the updated gain.'
  

  ;; Update setup
  mike_redux_updset, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_doslitf, state

  ;; Side
  side = (where(state.side,nside) + 1)

  ;; Run
  mike_slitflat, state.struct, state.setval, side, CLOBBER=state.clobber, $
    CHK=state.flg_objopt[1]

  ;; Update setup
  mike_redux_updset, state

  return
end
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_doarcs, state

  ;; Side
  side = (where(state.side,nside) + 1)
  
  ;; Run
  tmp = state.struct
  mike_allarc, tmp, state.setval, side, CLOBBER=state.clobber
  state.struct = temporary(tmp)

  ;; Update setup
  mike_redux_updset, state
  mike_redux_updobj, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_doobj, state

  ;; Side
  side = (where(state.side,nside) + 1)
  
  if state.obj EQ 0 then begin
      print, 'mike_redux:  Select an OBJ first!'
      return
  endif

  ;; Run
  mike_allobj, state.struct, state.setval, state.obj, side, /PROCALL, $
    DOCR=state.flg_objopt[2], CLOBBER=state.flg_objopt[3], $
    XOPT=state.flg_objopt[0], CHK=state.flg_objopt[1]

  ;; Update setup
  mike_redux_updset, state
  mike_redux_updobj, state

  return
end
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_changeset, state, indx

  ;; Setup
  if state.flg_setup EQ 1 then begin
      gd = where(state.flg_set)
      setup = state.setup[gd[indx]]
  endif else begin 
      setup = state.setup[indx]
  endelse
  state.setval = setup
  
  mtch = where(state.struct.setup EQ setup)
  idx = mtch[0]
  
  ;; Set Angles
;  widget_control, state.xdangl_id, set_value=state.struct[idx].xdangl
;  widget_control, state.echangl_id, set_value=state.struct[idx].echangl
 
  ;; Binning
  widget_control, state.binning_id, $
    set_value=strtrim(state.struct[idx].colbin,2)+ $
    'x'+strtrim(state.struct[idx].rowbin,2)

  ;; Decker
;  widget_control, state.decker_id, set_value=state.struct[idx].decker

  ;; Blocking
;  widget_control, state.block_id, set_value=state.struct[idx].block

  ;; Cross
;  widget_control, state.cross_id, set_value=state.struct[idx].cross

  ;; Obj
  state.obj = 0L
  gd = where(state.struct.setup EQ setup AND $
             state.struct.obj_id GT 0, ngd)
  if ngd NE 0 then begin
      obj_id = state.struct[gd[uniq(state.struct[gd].obj_id, $
                                    sort(state.struct[gd].obj_id))]].obj_id
      nobj = n_elements(obj_id)
      ;; Create string
      allobj = ['']
      for qq=0L,nobj-1 do begin
          mtch = where(state.struct.setup EQ setup AND $
                       state.struct.obj_id EQ obj_id[qq])
          allobj = [allobj, state.struct[mtch[0]].Obj]
      endfor
      allobj = allobj[1:*]
      ;; GUI
      widget_control, state.obj_id, set_value=strtrim(allobj,2)
  endif else widget_control, state.obj_id, set_value=''

  mike_redux_updset, state

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_updset, state

  setup = state.setval 
  
  ;; Arcs
  gd = where(state.struct.setup EQ setup AND state.struct.flg_anly NE 0 AND $
             state.struct.type EQ 'ARC',ngd)
  if ngd NE 0 then begin
      frame = state.struct[gd[uniq(state.struct[gd].frame, $
                                   sort(state.struct[gd].frame))]].frame
      nframe = n_elements(frame)
      ;; Create string
      allarc = ['']
      for qq=0L,nframe-1 do begin
          mtch = where(state.struct.setup EQ setup AND $
                       state.struct.frame EQ frame[qq])
          allarc = [allarc, state.struct[mtch[0]].img_root]
      endfor
      allarc = allarc[1:*]
      ;; GUI
      widget_control, state.arcs_id, set_value=strtrim(allarc,2)
      ;; Save
      narc = n_elements(allarc)
      state.narc = narc
      state.arc_fil[0:narc-1] = allarc
  endif else widget_control, state.arcs_id, set_value=''

  ;; Flats
  gd = where(state.struct.setup EQ setup AND state.struct.flg_anly NE 0 AND $
             state.struct.type EQ 'TFLT',ngd)
  if ngd NE 0 then begin
      frame = state.struct[gd[uniq(state.struct[gd].frame, $
                                   sort(state.struct[gd].frame))]].frame
      nframe = n_elements(frame)
      ;; Create string
      allflat = ['']
      for qq=0L,nframe-1 do begin
          mtch = where(state.struct.setup EQ setup AND $
                       state.struct.frame EQ frame[qq])
          allflat = [allflat, state.struct[mtch[0]].img_root]
      endfor
      allflat = allflat[1:*]
      ;; GUI
      widget_control, state.flats_id, set_value=strtrim(allflat,2)
      ;; Save
      nflat = n_elements(allflat)
      state.nflat = nflat
      state.flat_fil[0:nflat-1] = allflat
  endif else widget_control, state.flats_id, set_value=''

  ;;;;;;;;;;;
  ;; Steps

  state.nqa = 0L

  ;; Sides
  for qq=1L,2L do begin
      
      gd = where(state.struct.flg_anly NE 0 and state.struct.setup EQ setup $
                 and state.struct.side EQ qq, ngd)
      if ngd EQ 0 then continue
      ;; Gain
      rdxi = where(state.rdx_steps EQ 'Gain')
      if abs(100*state.struct[gd[0]].gain - $
        round(100*state.struct[gd[0]].gain)) GT 1e-3 then begin
          state.rdx_sel[qq-1,rdxi] = 1
      endif else state.rdx_sel[qq-1,rdxi] = 0

      ;; Milky flat
      fil = mike_getfil('mflat_fil', setup, SIDE=qq, /name, CHKFIL=chkf)
      rdxi = where(state.rdx_steps EQ 'MFlat Flat')
      if chkf then state.rdx_sel[qq-1,rdxi] = 1 $
      else state.rdx_sel[qq-1,rdxi] = 0


      ;; Proc flat
      fil = mike_getfil('tflat_fil', setup, SIDE=qq, /name, CHKFIL=chkf)
      rdxi = where(state.rdx_steps EQ 'TFlat Flat')
      if chkf then state.rdx_sel[qq-1,rdxi] = 1 $
      else state.rdx_sel[qq-1,rdxi] = 0

      ;; Edge flat
      fil = mike_getfil('ordr_str', setup, SIDE=qq, /name, CHKFIL=chkf)
      rdxi = where(state.rdx_steps EQ 'Edge Flat')
      if chkf then begin
          state.rdx_sel[qq-1,rdxi] = 1 
          ;; QA
          qafil = mike_getfil('qa_trcflat', setup, SIDE=qq)
          state.qa_fil[state.nqa] = qafil
          state.nqa = state.nqa + 1
      endif else state.rdx_sel[qq-1,rdxi] = 0

      ;; Slit flat
      fil = mike_getfil('qa_slitflat', setup, SIDE=qq, /name, CHKFIL=chkf)
      rdxi = where(state.rdx_steps EQ 'Slit Flat')
      if chkf then begin
          state.rdx_sel[qq-1,rdxi] = 1 
          qafil = mike_getfil('qa_slitflat', setup, SIDE=qq)
          state.qa_fil[state.nqa] = qafil
          state.nqa = state.nqa + 1
      endif else state.rdx_sel[qq-1,rdxi] = 0

      ;; GUI
      gd= where(state.rdx_sel[qq-1,*])
      widget_control, state.rdx_id[qq-1], set_list_select=gd
      if state.nqa NE 0 then $
        widget_control, state.qa_id, set_value=state.qa_fil[0:state.nqa-1]
  endfor

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_redux_updobj, state, indx

  setup = state.setval
  
  ;; Obj files
  gd = where(state.struct.setup EQ setup AND state.struct.flg_anly NE 0 AND $
             state.struct.obj_id EQ state.obj AND $
             (state.struct.type EQ 'OBJ' OR state.struct.type EQ 'STD'),ngd)
  if ngd NE 0 then begin
      frame = state.struct[gd[uniq(state.struct[gd].frame, $
                                   sort(state.struct[gd].frame))]].frame
      nframe = n_elements(frame)
      ;; Create string
      allobj = ['']
      for qq=0L,nframe-1 do begin
          mtch = where(state.struct.setup EQ setup AND $
                       state.struct.frame EQ frame[qq])
          allobj = [allobj, state.struct[mtch[0]].img_root]
      endfor
      allobj = allobj[1:*]
      ;; GUI
      widget_control, state.objf_id, set_value=strtrim(allobj,2)
  endif else widget_control, state.objf_id, set_value=''

  ;; Sides
  state.nqa = state.nqa < 2
  for qq=1L,3L do begin
      
      gd = where(state.struct.setup EQ setup AND $
                 state.struct.flg_anly NE 0 AND $
                 state.struct.side EQ qq AND $
                 state.struct.obj_id EQ state.obj AND $
                 (state.struct.type EQ 'OBJ' OR $
                  state.struct.type EQ 'STD'),ngd)
      if ngd EQ 0 then continue
      afil = state.struct[gd[0]].arc_fil
      rfil = state.struct[gd[0]].img_root
      frmobj = state.struct[gd[0]].frame
      pos = strpos(afil, '.fits')
      frame = long(strmid(afil,pos-4,4))
      
      ;; Proc Arc
      rdxi = where(state.rdx_steps EQ 'Proc Arc')
      if x_chkfil(afil+'*',/silent) then state.rdx_sel[qq-1,rdxi] = 1 $
      else state.rdx_sel[qq-1,rdxi] = 0

      ;; Fit Arc
      rdxi = where(state.rdx_steps EQ 'Fit Arc')
      out_fil = mike_getfil('arc_fit', setup, SIDE=qq, $
                             /name, CHKFIL=chkfil, SUBFIL=afil)
      if chkfil then begin
          state.rdx_sel[qq-1,rdxi] = 1 
          qafil = mike_getfil('qa_arcfit', setup, SIDE=qq, $
                                 /name, CHKFIL=chkfil, SUBFIL=afil)
          state.qa_fil[state.nqa] = qafil
          state.nqa = state.nqa + 1
      endif else state.rdx_sel[qq-1,rdxi] = 0
      
      ;; Fit2D Arc
      rdxi = where(state.rdx_steps EQ '2DFit Arc')
      out_fil = mike_getfil('arc_2Dfit', setup, SIDE=qq, $
                             /name, CHKFIL=chkfil, SUBFIL=afil)
      if chkfil then begin
          state.rdx_sel[qq-1,rdxi] = 1 
          qafil = mike_getfil('qa_arc2dfit', setup, SIDE=qq, $
                                 /name, CHKFIL=chkfil, SUBFIL=afil)
          state.qa_fil[state.nqa] = qafil
          state.nqa = state.nqa + 1
      endif $
      else state.rdx_sel[qq-1,rdxi] = 0
      
      ;; Trace Arc
      rdxi = where(state.rdx_steps EQ 'Trace Arc')
      out_fil = mike_getfil('arc_trc', setup, SIDE=qq, $
                             /name, CHKFIL=chkfil, SUBFIL=afil)
      if chkfil then begin
          state.rdx_sel[qq-1,rdxi] = 1 
          qafil = mike_getfil('qa_tracearc', setup, SIDE=qq, $
                                 /name, CHKFIL=chkfil, SUBFIL=afil)
          state.qa_fil[state.nqa] = qafil
          state.nqa = state.nqa + 1
      endif $
      else state.rdx_sel[qq-1,rdxi] = 0

      ;; Fit Trace Arc
      rdxi = where(state.rdx_steps EQ 'FitTrc Arc')

      out_fil = mike_getfil('arc_fittrc', setup, SIDE=qq, $
                             /name, CHKFIL=chkfil, SUBFIL=afil)
      if chkfil then begin
          state.rdx_sel[qq-1,rdxi] = 1 
          qafil = mike_getfil('qa_fittrcarc', setup, SIDE=qq, $
                                 /name, CHKFIL=chkfil, SUBFIL=afil)
          state.qa_fil[state.nqa] = qafil
          state.nqa = state.nqa + 1
      endif else state.rdx_sel[qq-1,rdxi] = 0

      ;; Arc Image
      rdxi = where(state.rdx_steps EQ 'Arc Image')
      out_fil = mike_getfil('arc_img', setup, SIDE=qq, $
                             /name, CHKFIL=chkfil, SUBFIL=afil)
      if chkfil then state.rdx_sel[qq-1,rdxi] = 1 $
      else state.rdx_sel[qq-1,rdxi] = 0

      ;; Proc Obj
      rdxi = where(state.rdx_steps EQ 'Proc Obj')
      outfil =mike_getfil('fin_fil', /name, SIDE=qq, $
                           CHKFIL=chkfil, SUBFIL=rfil)
      if chkfil then state.rdx_sel[qq-1,rdxi] = 1 $
      else state.rdx_sel[qq-1,rdxi] = 0

      ;; Trace Obj
      rdxi = where(state.rdx_steps EQ 'Trace Obj')
      objnm = mike_getfil('obj_fil', SIDE=qq, $
                            chkfil=chkfil,/name, SUBFIL=rfil)
      if chkfil then begin
          objstr = xmrdfits(objnm,1,/silent)
          state.rdx_sel[qq-1,rdxi] = 1 
          qafil = mike_getfil('qa_fntobj', setup, SIDE=qq, $
                               /name, CHKFIL=chkfil, SUBFIL=rfil)
          state.qa_fil[state.nqa] = qafil
          state.nqa = state.nqa + 1
      endif $
      else state.rdx_sel[qq-1,rdxi] = 0

      ;; Skysub
      rdxi = where(state.rdx_steps EQ 'Skysub')
      skyfil =mike_getfil('sky_fil', /name, SIDE=qq,$
                           CHKFIL=chkfil, SUBFIL=rfil)
      if chkfil then state.rdx_sel[qq-1,rdxi] = 1 $
      else state.rdx_sel[qq-1,rdxi] = 0

      ;; Extract
      rdxi = where(state.rdx_steps EQ 'Extract')
      if keyword_set(OBJSTR) then begin
          if total(objstr.fx) gt 0. then begin
              state.rdx_sel[qq-1,rdxi] = 1 
              qafil = mike_getfil('qa_extract', setup, SIDE=qq, $
                                   /name, CHKFIL=chkfil, SUBFIL=rfil)
              state.qa_fil[state.nqa] = qafil
              state.nqa = state.nqa + 1
          endif else state.rdx_sel[qq-1,rdxi] = 0
      endif

      ;; GUI
      gd= where(state.rdx_sel[qq-1,*])
      widget_control, state.rdx_id[qq-1], set_list_select=gd
      if keyword_set(OBJSTR) then delvarx, objstr
  endfor

  if state.nqa NE 0 then $
    widget_control, state.qa_id, set_value=state.qa_fil[0:state.nqa-1]

  return
end
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_redux_init, state

  ;; Setups
  if state.flg_setup EQ 1 then begin
      ;; Restict to obj
      gd = where(state.struct.obj_id GT 0,ngd)
      if ngd EQ 0 then stop
      gdset = state.struct[gd[uniq(state.struct[gd].setup, $
                                sort(state.struct[gd].setup))]].setup
      nset = n_elements(gdset)
      state.flg_set = 0B
      for jj=0L,nset-1 do begin
          a = where(state.setup EQ gdset[jj])
          state.flg_set[a] = 1B
      endfor
  endif else begin
      gdset = state.setup
  endelse
  ;; Gui
  widget_control, state.setup_id, set_value=strtrim(gdset,2)

  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro mike_redux, mike, OBJ=obj, HIRES_FIL=mike_fil


;
  if  N_params() LT 1 then begin 
    print,'Syntax - ' + $
      'mike_redux, mike, /OBJ  [v1.0]'
    return
  endif 

  if not keyword_set(mike_fil) then begin
      a = findfile('./mike_*fits*', count=count)
      if count EQ 0 then a = 'mike_save.fits'
      print, 'mike_redux: Writing to '+a[0]
      mike_fil = a[0]
  endif

  ;; Setups
  gdset = mike[uniq(mike.setup, sort(mike.setup))].setup
      
; STATE
  state = {             $
            rdx_steps: strarr(50), $
            qa_fil: strarr(50), $
            nqa: 0L, $
            qa_sel: -1L, $
            nsteps: 0L, $
            arc_fil: strarr(100), $
            narc: 0L, $
            flat_fil: strarr(500), $
            nflat: 0L, $
            blue_rdx: lonarr(50),$ 
            rdx_id: lonarr(3), $
            rdx_sel: lonarr(3,50), $
            struct: mike, $
            setup: gdset, $
            setval: 0L, $
            side: replicate(1B,2L), $
            outfil: mike_fil, $
            flg_objopt: intarr(10), $
            flg_set: bytarr(n_elements(gdset)), $
            path: '', $
            flg_setup: 0, $  ; 1=obj, 0=Fspec
            obj: 0L, $   ; OBJ ID
            nogain: 0, $
            clobber: 0, $
            base_id: 0L, $      ; Widgets
;            xdangl_id: 0L, $
;            echangl_id: 0L, $
            binning_id: 0L, $
;            block_id: 0L, $
;            decker_id: 0L, $
;            cross_id: 0L, $
            arcs_id: 0L, $
            flats_id: 0L, $
            qa_id: 0L, $
            setup_id: 0L, $
            obj_id: 0L, $
            objf_id: 0L, $
            sides_id: 0L, $
            objopt_id: 0L, $
            options_id: 0L, $
            bluerdx_id: 0L, $
            greenrdx_id: 0L, $
            redrdx_id: 0L, $
            exp_id: 0L $
          }

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    WIDGET
  state.base_id = WIDGET_BASE( title = 'mike_redux', /column, $
                    UNAME='BASE', xoffset=100L, yoffset=100L)
  base= state.base_id
  
;;;;;;;;;;;;;;
; 
  mainbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
  ;; Config
;  configbox = WIDGET_BASE( mainbar, /column, /frame, /base_align_center,$
;                         /align_center)
;  state.xdangl_id = cw_field(configbox, value=0.d, $
;                             xsize=7, title='XDANGL', UVALUE='XDANGL')
;  state.echangl_id = cw_field(configbox, value=0.d, $
;                             xsize=7, title='EANGL', UVALUE='EANGL')
;  state.decker_id = cw_field(configbox, value='', $
;                             xsize=2, title='DECKER', UVALUE='DECKER')
;  state.block_id = cw_field(configbox, value='', $
;                             xsize=5, title='BLOCK', UVALUE='BLOCK')
;  state.cross_id = cw_field(configbox, value='', $
;                             xsize=5, title='CROSS', UVALUE='CROSS')
  ;; Setup
  setupbox = WIDGET_BASE( mainbar, /column, /frame, /base_align_center,$
                         /align_center)
  setuplabel = WIDGET_LABEL(setupbox, value='Setups', /align_center)
  state.setup_id = WIDGET_LIST(setupbox, $
                               VALUE=strarr(10),  uvalue='SETUPLIST', $
                               xsize=2, ysize = 4)
  setupbut = CW_BGROUP(setupbox, ['OBJONLY'], /nonexclusive, $
                       /return_index, uvalue='SETUPBUT')
  state.binning_id = cw_field(setupbox, value='', $
                             xsize=3, title='BINNING', UVALUE='BINNING')
  
  ;; Objects
  objbox = WIDGET_BASE( mainbar, /column, /frame, /base_align_center,$
                         /align_center)
  objlabel = WIDGET_LABEL(objbox, value='Objects', /align_center)
  state.obj_id = WIDGET_LIST(objbox, $
                              VALUE=strarr(10),  uvalue='OBJLIST', $
                             xsize=15,ysize = 3)
  widget_control, state.obj_id, set_list_select=-1
  state.objf_id = WIDGET_LIST(objbox, xsize=15, $
                              VALUE=strarr(5),  uvalue='OBJLIST', ysize = 3)
  widget_control, state.objf_id, set_list_select=-1

  ;; CALIBS
  calibbox = WIDGET_BASE( mainbar, /column, /frame, /base_align_center,$
                         /align_center)
  ;; Arcs
  arcbox = WIDGET_BASE( calibbox, /column, /base_align_center,/align_center)
  objlabel = WIDGET_LABEL(arcbox, value='Arcs', /align_center)
  state.arcs_id = WIDGET_LIST(arcbox, xsize=15, /multiple, $
                              VALUE=strarr(10),  uvalue='ARCLIST', ysize = 3)
  widget_control, state.arcs_id, set_list_select=-1

  ;; Flats
  flatbox = WIDGET_BASE( calibbox, /column, /base_align_center,/align_center)
  flatlabel = WIDGET_LABEL(flatbox, value='Flats', /align_center)
  state.flats_id = WIDGET_LIST(flatbox, xsize=15,/multiple,  $
                              VALUE=strarr(10),  uvalue='FLATLIST', ysize = 3)
  widget_control, state.flats_id, set_list_select=-1
;  FLGANLY = WIDGET_BUTTON(calibbox, value='TOSSOUT',uvalue='TOSSOUT')

  ;; QA
  qabox = WIDGET_BASE( mainbar, /column, /base_align_center,/align_center)
  qalabel = WIDGET_LABEL(qabox, value='QA', /align_center)
  showqa = WIDGET_BUTTON(qabox, value='SHOW',uvalue='SHOWQA')
  state.qa_id = WIDGET_LIST(qabox, xsize=25, $
                              VALUE=state.qa_fil,  uvalue='QALIST', ysize = 9)
  widget_control, state.qa_id, set_list_select=-1

  ;; Steps
  procbox = WIDGET_BASE( state.base_id, /column, /frame, /base_align_center,$
                         /align_center)
  calbox = WIDGET_BASE(procbox, /row, /frame, /base_align_center,$
                         /align_center)
  state.sides_id = CW_BGROUP(calbox, ['B','R'], /nonexclusive, $
                       /return_index, uvalue='SIDES',/row)
  widget_control, state.sides_id, set_value=[1,1]
  
  dogain = WIDGET_BUTTON(calbox, value='DOGAIN',uvalue='DOGAIN')
  flats = WIDGET_BUTTON(calbox, value='DOFLATS',uvalue='DOFLATS')
  doarcs = WIDGET_BUTTON(calbox, value='DOARCS',uvalue='DOARCS')
  dosflat = WIDGET_BUTTON(calbox, value='DOSLITF',uvalue='DOSLITF')

;  doobjbox = WIDGET_BASE(procbox, /row, /frame, /base_align_center,$
;                         /align_center)
  doobj = WIDGET_BUTTON(calbox, value='DOOBJ',uvalue='DOOBJ')
  state.objopt_id = CW_BGROUP(calbox, ['XOPT','CHK','DOCR','CLOBBER'], $
                              /nonexclusive, $
                              /return_index, /row, uvalue='OBJOPTION')
  ;; OTHER
  otherbox = WIDGET_BASE(procbox, /row, /frame, /base_align_center,$
                         /align_center)
  state.options_id = cw_field(otherbox, value='', $
                             xsize=30, title='OPTIONS', UVALUE='OPTIONS')
  ESTRCT = WIDGET_BUTTON(otherbox, value='EDITSTRCT',uvalue='EDITSTRCT')
  WRITE = WIDGET_BUTTON(otherbox, value='WRITE',uvalue='WRITE')
  SPLOT = WIDGET_BUTTON(otherbox, value='SPECPLOT',uvalue='SPECPLOT')
  done = WIDGET_BUTTON(otherbox, value='DONE',uvalue='DONE')
  dnnsv = WIDGET_BUTTON(otherbox, value='DNSV',uvalue='DNSV')

  ;; Redux list
  rdx_stps = ['Gain', $
              'MFlat Flat', $
              'TFlat Flat', $
              'Edge Flat', $
              'Proc Arc', $
              'Fit Arc', $
              '2DFit Arc', $
              'Trace Arc',$ 
              'FitTrc Arc',$ 
              'Arc Image',$ 
              'Slit Flat',$ 
              'Proc Obj',$ 
              'Trace Obj',$ 
              'Skysub',$ 
              'Extract',$ 
              'Dum' ]
  state.nsteps = n_elements(rdx_stps)
  state.rdx_steps[0:state.nsteps-1] = rdx_stps
  
  ;; BGR
  bgrbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                         /align_center)
  ;; Blue
  bluebar = WIDGET_BASE( bgrbar, /column, /frame, /base_align_center,$
                         /align_center, title='Blue')
  bluelabel = WIDGET_LABEL(bluebar, value='Blue', /align_center)
  state.bluerdx_id = WIDGET_LIST(bluebar, xsize=15, $
                                 VALUE=state.rdx_steps[0:state.nsteps-1],  $
                                 uvalue='BLUESTPS',  ysize=15, $
                                /multiple)

  state.rdx_id[0] = state.bluerdx_id
  
  ;; Green
;  greenbar = WIDGET_BASE( bgrbar, /column, /frame, /base_align_center,$
;                         /align_center, title='Green')
;  greenlabel = WIDGET_LABEL(greenbar, value='Green', /align_center)
;  state.greenrdx_id = WIDGET_LIST(greenbar, xsize=15, $
;                                 VALUE=state.rdx_steps[0:state.nsteps-1],  $
;                                 uvalue='GREENSTPS',  ysize=15, $
;                                /multiple)
;  state.rdx_id[1] = state.greenrdx_id

  ;; Red
  redbar = WIDGET_BASE( bgrbar, /column, /frame, /base_align_center,$
                         /align_center, title='Red')
  redlabel = WIDGET_LABEL(redbar, value='Red', /align_center)
  state.redrdx_id = WIDGET_LIST(redbar, xsize=15, $
                                 VALUE=state.rdx_steps[0:state.nsteps-1],  $
                                 uvalue='REDSTPS',  ysize=15, $
                                /multiple)
  state.rdx_id[1] = state.redrdx_id

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      BUTTONS
;  butbase = widget_base(listbase, /column, /align_center)

; Realize
  WIDGET_CONTROL, base, /realize

  ;; Initialize
  mike_redux_init, state

  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

; Send to the xmanager
  xmanager, 'mike_redux', base, /no_block

  !P.MULTI= [0,1,1]
  print, 'mike_redux: WARNING!!!! --  It is likely that the mike structure '
  print, 'mike_redux: in memory is out of date.  '
  print, 'mike_redux: Read in the saved fits file after finishing !!'
  return
end

