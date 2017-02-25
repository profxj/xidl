;+ 
; NAME:
; hires_procarc   
;     Version 1.1
;
; PURPOSE:
;    Process (bias subtract, flatten) the arc files.  In addition, the
;    code first calls hires_arcxyoff which uses a 2D FFT to determine
;    the offset between the Template Arc and the current arc due to
;    thermal expansion in the instrument.
;   
;    This file also includes the routine hires_procarc_sngl which
;    allows the processing of a single Arc image given the filename of
;    the Raw image.
;
; CALLING SEQUENCE:
;   
;  hires_procarc, hires, setup, obj_id, [chip], ATEMPL=, /CLOBBER,
;  FLATFIL=
;  print, hires_procarc_sngl(raw_fil, hires, setup, [obj_id, chip,
;  exten], ATEMPL=, /CLOBBER)
;
; INPUTS:
;   hires     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [chip]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  One processed, combined Arc Image (e.g. Arcs/Arc_mb0539.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;   ATEMPL -- Index of the Template Arc image (default: 0L)
;   FLATFIL -- Filename of the milky flat (pixel to pixel correction)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_procarc, hires, 1, 1, 1, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_getfil
;  hires_arcxyoff
;  hires_proc
;
; REVISION HISTORY:
;   13-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hires_procarc_sngl, raw_fil, hires, setup, chip, exten, $
  TEMPL_FIL=templ_fil, FRAME=frame, $
  CLOBBER=clobber, FLATFIL=flatfil, CHK=chk, $
  XYOFF=xyoff, GAIN=gain, READNO=readno, ALIGN=align,$
  FITPRM=fitprm, _EXTRA=extra

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = hires_procarc_sngl( raw_fil, setup, chip, [exten], ' + $
        'FLATFIL=, /CLOBBER, /ALIGN ) [v1.1]'
      return, -1
  endif 
  
;  Optional Keywords
  compile_opt strictarr

  if not keyword_set( ATEMPL ) then atempl = 0L
  if not keyword_set(EXTEN) then exten = (chip > 0)  ;; Works for Single chip

  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = $
    '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa
  
  ;; Check for prior image
  if not keyword_set( FRAME ) then begin
      pos = strpos(raw_fil, '.fits')
      frame = long(strmid(raw_fil,pos-4,4))
  endif
  outfil = hires_getfil('arc_fil', FRAME=frame, CHIP=chip, $
                        /name, CHKFIL=chkfil)
  idx = where(hires.frame EQ frame AND hires.chip EQ chip AND hires.setup EQ setup)
  if CHKFIL NE 0 AND not keyword_set( CLOBBER ) then begin
      print, 'hires_procarc: Arc ', outfil, ' exists.  Returning..'
      return, outfil
  endif

  ;; Determine xy_off with the FFT
;  if keyword_set( TEMPL_FIL ) AND arg_present( XYOFF ) AND $
;    not keyword_set( ALIGN ) then begin
;      stop ;; Probably shouldnt be here
;      xyoff = hires_arcxyoff_work( templ_fil, raw_fil, chip )
;  endif

  ;; Process
  rslt = hires_proc_sngl(raw_fil, chip, exten, FLATFIL=flatfil, $
                         CLOBBER=clobber, $
                         REDOOV=clobber, /ARC, IOUTFIL=outfil, $
                         GAIN=hires[idx].gain, $
                         READNO=hires[idx].readno, SETUP=setup)

  ;; ALL DONE
  print, 'hires_procarc_sngl: All Done! '
  return, outfil

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_procarc, hires, setup, obj_id, chip, ATEMPL=atempl, $
                  CLOBBER=clobber, FLATFIL=flatfil, ALIGN=align

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_procarc, hires, setup, obj_id, [chip], FLATFIL=, /CLOBBER ' + $
        'ATEMPL= [v1.1]'
      return
  endif 
  
  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  if not keyword_set( CHIP ) then chip = [1L,2L,3L]
  if not keyword_set( ATEMPL ) then atempl = 0L
  
; Loop on chip
  nchip = n_elements(chip)
  for ii=0L,nchip-1 do begin

      qq = chip[ii]
      ;; CHIP
      case qq of 
          1: print, 'hires_procarc: Processing BLUE arc'
          2: print, 'hires_procarc: Processing GREEN arc'
          3: print, 'hires_procarc: Processing RED arc'
          else: stop
      endcase

      ;; Grab all Arc files
      arcs = where(hires.setup EQ setup AND hires.flg_anly NE 0 AND $
                   hires.chip EQ qq AND $
                   strtrim(hires.type,2) EQ 'ARC', narc)
      if narc EQ 0 then begin
          print, 'hires_procarc: No Arcs found! Returning' 
          return
      endif

      ;; Grab all obj files
      indx = where(hires.flg_anly NE 0 AND hires.chip EQ qq AND $
                   hires.obj_id EQ obj_id AND hires.setup EQ setup AND $
                   (strtrim(hires.type,2) EQ 'OBJ' OR $
                    strtrim(hires.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'hires_procarc: No Obj found!  Returning' 
          stop
          continue
      endif

      ;; LOOP
      for kk=0L,nindx-1 do begin

          print, 'hires_procarc: Index ', strtrim(kk,2)

          ;; Check for arc_fil name
          if strlen(strtrim(hires[indx[kk]].arc_fil,2)) EQ 0 then begin
              ;; Identify closest image in time
              diff = abs( (hires[indx[kk]].date + $
                           hires[indx[kk]].exp/2./3600./24.) - $
                          (hires[arcs].date + hires[arcs].exp/2./3600./24.) )
              mn = min(diff, gdarc)
              
              ;; Check for prior image
              outfil = hires_getfil('arc_fil', FRAME=hires[arcs[gdarc]].frame, $
                                    CHIP=qq, /name)
          endif else begin
              outfil=hires[indx[kk]].arc_fil 
              pos = strpos(outfil, '.fits')
              frame = long(strmid(outfil,pos-4,4))
              gdarc = where(hires[arcs].frame EQ frame,ngd)
              if ngd NE 1 then stop
          endelse
          
          if x_chkfil(outfil+'*',/silent) NE 0 $
            AND not keyword_set( CLOBBER ) then begin
              print, 'hires_procarc: Arc ', outfil, ' exists.  ' + $
                'Updating using previous offsets and then returning'
              hires[indx[kk]].arc_fil = outfil
              continue
          endif

          ;; Determine xy_off
;          hires_arcxyoff, hires, qq, arcs[atempl], arcs[gdarc], xyoff

          ;; Process
          hires_proc, hires, arcs[gdarc], FLATFIL=flatfil, CLOBBER=clobber, $
            REDOOV=clobber, /ARC, IOUTFIL=outfil

;          if keyword_set( ALIGN ) then begin
             ;; Grab the order structure
;             ordr_str = hires_getfil('ordr_str', setup, CHIP=chip)
             ;; Open the images
;             arc = xmrdfits(outfil, 0, /silent)
;             arci = xmrdfits(outfil, 1, /silent)
             ;; Here we go
;             xyoff = hires_arcalign_work(arc, arci, ordr_str, FITPRM=fitprm, /CHK)
;             hires[arcs[gdarc]].arc_xyoff = fitprm
;             hires[indx[kk]].arc_xyoff = fitprm
;           endif
      
          ;; Structure
          hires[indx[kk]].arc_fil = outfil

      endfor
  endfor

  ;; ALL DONE
  print, 'hires_procarc: All Done! '
  return

end

