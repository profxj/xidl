;+ 
; NAME:
; uves_procarc   
;     Version 1.1
;
; PURPOSE:
;    Process (bias subtract, flatten) the arc files.  In addition, the
;    code first calls uves_arcxyoff which uses a 2D FFT to determine
;    the offset between the Template Arc and the current arc due to
;    thermal expansion in the instrument.
;   
;    This file also includes the routine uves_procarc_sngl which
;    allows the processing of a single Arc image given the filename of
;    the Raw image.
;
; CALLING SEQUENCE:
;   
;  uves_procarc, uves, setup, obj_id, [chip], ATEMPL=, /CLOBBER,
;  FLATFIL=
;
; INPUTS:
;   uves     -  MIKE structure
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
;   uves_procarc, uves, 1, 1, 1, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_getfil
;  uves_arcxyoff
;  uves_proc
;
; REVISION HISTORY:
;   13-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function uves_procarc_sngl, raw_fil, setup, side, uves, $
  TEMPL_FIL=templ_fil, FRAME=frame, EXTEN=exten, $
  CLOBBER=clobber, FLATFIL=flatfil, CHK=chk, $
  XYOFF=xyoff, GAIN=gain, READNO=readno, ALIGN=align,$
  FITPRM=fitprm, _EXTRA=extra

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'rslt = uves_procarc_sngl( raw_fil, setup, side, uves' + $
        'FLATFIL=, /CLOBBER, /ALIGN ) [v1.1]'
      return, -1
  endif 
  
;  Optional Keywords
  compile_opt strictarr

  if not keyword_set( ATEMPL ) then atempl = 0L

  ;; Assuming Sara naming convention
  wcen = uves_getwcen(raw_fil,uves,FRAME=frame, SARA=sara)
  cwcen = strtrim(wcen,2)

  ;; QA
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = $
    '0'+strtrim(setup,2) 
  fqa = 'QA/Arcs'+cwcen+'_'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa
  
  ;; Check for prior image
  outfil = uves_getfil('arc_fil', setup, WCEN=wcen, FRAME=frame, $
                        /name, CHKFIL=chkfil)
  if CHKFIL NE 0 AND not keyword_set( CLOBBER ) then begin
      print, 'uves_procarc: Arc ', outfil, ' exists.  Returning..'
      return, outfil
  endif

  ;; Determine xy_off with the FFT
;  if keyword_set( TEMPL_FIL ) AND arg_present( XYOFF ) AND $
;    not keyword_set( ALIGN ) then begin
;      stop ;; Probably shouldnt be here
;      xyoff = uves_arcxyoff_work( templ_fil, raw_fil, chip )
;  endif

  ;; Process
  rslt = uves_proc_sngl(raw_fil, setup, side, EXTEN=exten, /NOFLAT, $
                        CLOBBER=clobber, $
                        REDOOV=clobber, /ARC, IOUTFIL=outfil)
;                        GAIN=1., READNO=3., $  ; Hard-wired because not important

  ;; ALL DONE
  print, 'uves_procarc_sngl: All Done! '
  return, outfil

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_procarc, uves, setup, obj_id, chip, ATEMPL=atempl, $
                  CLOBBER=clobber, FLATFIL=flatfil, ALIGN=align

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'uves_procarc, uves, setup, obj_id, [chip], FLATFIL=, /CLOBBER ' + $
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
          1: print, 'uves_procarc: Processing BLUE arc'
          2: print, 'uves_procarc: Processing GREEN arc'
          3: print, 'uves_procarc: Processing RED arc'
          else: stop
      endcase

      ;; Grab all Arc files
      arcs = where(uves.setup EQ setup AND uves.flg_anly NE 0 AND $
                   uves.chip EQ qq AND $
                   strtrim(uves.type,2) EQ 'ARC', narc)
      if narc EQ 0 then begin
          print, 'uves_procarc: No Arcs found! Returning' 
          return
      endif

      ;; Grab all obj files
      indx = where(uves.flg_anly NE 0 AND uves.chip EQ qq AND $
                   uves.obj_id EQ obj_id AND uves.setup EQ setup AND $
                   (strtrim(uves.type,2) EQ 'OBJ' OR $
                    strtrim(uves.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'uves_procarc: No Obj found!  Returning' 
          stop
          continue
      endif

      ;; LOOP
      for kk=0L,nindx-1 do begin

          print, 'uves_procarc: Index ', strtrim(kk,2)

          ;; Check for arc_fil name
          if strlen(strtrim(uves[indx[kk]].arc_fil,2)) EQ 0 then begin
              ;; Identify closest image in time
              diff = abs( (uves[indx[kk]].date + $
                           uves[indx[kk]].exp/2./3600./24.) - $
                          (uves[arcs].date + uves[arcs].exp/2./3600./24.) )
              mn = min(diff, gdarc)
              
              ;; Check for prior image
              outfil = uves_getfil('arc_fil', FRAME=uves[arcs[gdarc]].frame, $
                                    CHIP=qq, /name)
          endif else begin
              outfil=uves[indx[kk]].arc_fil 
              pos = strpos(outfil, '.fits')
              frame = long(strmid(outfil,pos-4,4))
              gdarc = where(uves[arcs].frame EQ frame,ngd)
              if ngd NE 1 then stop
          endelse
          
          if x_chkfil(outfil+'*',/silent) NE 0 $
            AND not keyword_set( CLOBBER ) then begin
              print, 'uves_procarc: Arc ', outfil, ' exists.  ' + $
                'Updating using previous offsets and then returning'
              uves[indx[kk]].arc_fil = outfil
              continue
          endif

          ;; Determine xy_off
;          uves_arcxyoff, uves, qq, arcs[atempl], arcs[gdarc], xyoff

          ;; Process
          uves_proc, uves, arcs[gdarc], FLATFIL=flatfil, CLOBBER=clobber, $
            REDOOV=clobber, /ARC, IOUTFIL=outfil

;          if keyword_set( ALIGN ) then begin
             ;; Grab the order structure
;             ordr_str = uves_getfil('ordr_str', setup, CHIP=chip)
             ;; Open the images
;             arc = xmrdfits(outfil, 0, /silent)
;             arci = xmrdfits(outfil, 1, /silent)
             ;; Here we go
;             xyoff = uves_arcalign_work(arc, arci, ordr_str, FITPRM=fitprm, /CHK)
;             uves[arcs[gdarc]].arc_xyoff = fitprm
;             uves[indx[kk]].arc_xyoff = fitprm
;           endif
      
          ;; Structure
          uves[indx[kk]].arc_fil = outfil

      endfor
  endfor

  ;; ALL DONE
  print, 'uves_procarc: All Done! '
  return

end

