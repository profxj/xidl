;+ 
; NAME:
; apf_getfil   
;     Version 1.1
;
; PURPOSE:
;   Pass back a structure and/or filename for a specified file type in
;   the HIRES pipeline.
;
; CALLING SEQUENCE:
;   
;  rslt = apf_getfil('type', [setup], CHIP=, /NAME, SUBFIL=, CHKFIL=,
;  SZ=, INDX=, FIL_NM=, HEAD=, FRAME=)
;
; INPUTS:
;   [setup]   -  Setup identifier 
;
; RETURNS:
;  Structure, image, name, etc.
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /name   -- Only return resolved name (string)
;   CHKFIL  -- Value equal to the number of files matching name
;   SZ      -- Image size.  Used to determine binning
;   SUBFIL  -- Image name generally used to parse the root name of the
;              image  (e.g.  'Arcs/arc_mb0539.fits').  Required in
;              many cases.
;   CHIP    -- Specify chip [1,2,3]
;   INDX    -- Image extension in the fits file (generally 0, 1, or 2)
;   FRAME=  -- Frame number for the image
;
; OPTIONAL OUTPUTS:
;   FIL_NM  -- Filename of the file 
;   HEAD    -- Image header
;
; COMMENTS:
;
; EXAMPLES:
;   ordr_str = apf_getfil('ordr_str', 1, CHIP=1)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------

function apf_getfil, type, setup, CHIP=chip, NAME=name, SUBFIL=subfil, $
                      CHKFIL=chkfil, SZ=sz, INDX=indx, FIL_NM=fil_cchip, $
                      HEAD=head, FRAME=frame


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = apf_getfil(type, [setup], CHIP=, /NAME, SUBFIL=, CHKFIL=,'
      print, '     SZ=, INDX=, FIL_NM=, HEAD=, FRAME=) [v1.1]'
      return, -1
  endif 

  ;; Setup
  if keyword_set( SETUP ) then begin
      if setup GE 10 then c_s = strtrim(setup,2) $
      else c_s = '0'+strtrim(setup,2) 
  endif


  ;; Size
  if keyword_set( SZ ) then begin
      cbin = round(2048. / sz[0])
      rbin = round(4096. / sz[1])
  endif

  ;;  Assume apf files have this type of file name:  apf0000.fits 

  if keyword_set( SUBFIL ) then begin
      dot = rstrpos( SUBFIL[0], '.')
      
      expname = strmid(subfil,usepos,dot-usepos)
      mname  = 'm' + expname
      mfits  = 'm' + expname + '.fits'
  endif

  ;; Side
  cchip = 'S'

  if keyword_set( FRAME ) then begin
      dum = frame
      cframe = ''
      while( DUM LT 10000 ) do begin
          dum = dum*10
          cframe = cframe+'0'
      endwhile
      cframe = cframe+strtrim(frame,2)
  endif

  ;; Big case statement
  case type of
      'arc_fil': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Arcs/Arc_'+cchip+cframe+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'arc_fit': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Arcs/Fits/Arc_'+cchip+cframe+'_fit.idl'
          NAME = 1
      end
      'arc_2Dfit': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Arcs/Fits/Arc_'+cchip+cframe+'_fit2D.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_trc': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Arcs/TRC/Arc_'+cchip+cframe+'_T.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_fittrc': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Arcs/TRC/Arc_'+cchip+cframe+'_F.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_mkaimg': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Arcs/Arc_'+cchip+cframe+'_I.fits'
      end
      ;; Bias
      'bias_fil': begin
          if not keyword_set( SZ ) then stop
          ;if not keyword_set( CHIP ) then stop
          fil_cchip = 'Bias/Bias'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+cchip+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'ov_fil': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'OV/ov_'+cchip+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Final file
      'fin_fil': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Final/f_apf'+cframe+cchip+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'fspec_fil': begin
          if not keyword_set( SUBFIL ) then stop
          ;if not keyword_set( CHIP ) then stop
          fil_cchip = 'FSpec/'+strtrim(subfil)+cchip+'.fits' 
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Obj file
      'obj_fil': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Extract/Obj_'+cframe+cchip+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'objN_fil': begin
          if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Extract/Obj_'+cframe+cchip+'N.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'sens_fil': begin
          if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Extract/sens_'+cframe+cchip+'.idl'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Sky file
      'sky_fil': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Sky/sky_apf'+cframe+cchip+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'scatt_fil': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'Sky/scatt_apf'+cframe+cchip+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; FLATS
      'qtz_fil': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cchip = 'Flats/TFlat_'+cchip+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'nqtz_fil': begin  ;; Normalized flat
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cchip = 'Flats/TFlat_N'+cchip+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'pixflt_fil': begin  ;; Normalized flat
          if not keyword_set( CHIP ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cchip = 'Flats/PFlat_N'+cchip+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'tflat_str': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cchip = 'Flats/TStr_'+cchip+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Order structure
      'ordr_str': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cchip = 'Flats/OStr_'+cchip+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Logs
      'flat_log': begin
          NAME = 1
          fil_cchip = 'Logs/flat_redux.log'
      end
      ;; QA
      'qa_trcflat': begin
          NAME = 1
          fil_cchip = 'QA/Flats'+c_s+'/qa_trcflt_'+c_s+cchip+'.ps'
      end
      'qa_slitflat': begin
          NAME = 1
          fil_cchip = 'QA/Flats'+c_s+'/qa_slitflat_'+c_s+cchip+'.ps'
      end
      'qa_fittflat': begin
          NAME = 1
          fil_cchip = 'QA/Flats'+c_s+'/qa_fittflt_'+c_s+cchip+'.ps'
      end
      'qa_arcfit': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'QA/Arcs'+c_s+'/qa_arcfit_'+cchip+cframe+'.ps'
          NAME = 1
      end
      'qa_arc2dfit': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'QA/Arcs'+c_s+'/qa_arc2dfit_'+cchip+cframe+'.ps'
          NAME = 1
      end
      'qa_tracearc': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'QA/Arcs'+c_s+'/qa_tracearc_'+cchip+cframe+'.ps'
          NAME = 1
      end
      'qa_fittrcarc': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cchip = 'QA/Arcs'+c_s+'/qa_fitrcarc_'+cchip+cframe+'.ps'
          NAME = 1
      end
      'qa_fntobj': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_cchip = 'QA/Obj'+c_s+'/qa_fntobj_'+cframe+cchip+'.ps'
      end
      'qa_skysub': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_cchip = 'QA/Obj'+c_s+'/qa_skysub_'+cframe+cchip+'.ps'
      end
      'qa_extract': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_cchip = 'QA/Obj'+c_s+'/qa_extract_'+cframe+cchip+'.ps'
      end
      'qa_chkflux': begin
          ;if not keyword_set( CHIP ) then stop
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_cchip = 'QA/Obj'+c_s+'/qa_flux_'+cframe+cchip+'.ps'
      end
      else: stop
  endcase

  ;; CHKFIL
  if arg_present( CHKFIL ) then a = findfile(fil_cchip+'*', count=CHKFIL)

  ;; Return
  if keyword_set( NAME ) then return, fil_cchip $
  else begin
      ;; Check for file
      if not keyword_set( CHKFIL ) then a = findfile(fil_cchip+'*', count=CHKFIL)
      if CHKFIL NE 0 then return, xmrdfits(fil_cchip, indx, head, /silent) $
      else begin
          case type of
              'arc_trc': $
                print, 'File ', fil_cchip, ' does not exist. Run apf_tracearc!'
              'arc_2Dfit': $
                print, 'File', fil_cchip, ' does not exist. Run apf_fit2darc!'
              'arc_fittrc': $
                print, 'File', fil_cchip, ' does not exist. Run apf_fittrcarc!'
              'flat_fil': $
                print, 'File', fil_cchip, ' does not exist. Run apf_mktflat!'
              else: print, 'apf_getfil: File does not exist! ', fil_cchip
          endcase
          stop
      endelse
      return, -1
  endelse

  return, -1
end

