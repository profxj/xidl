;+ 
; NAME:
; hamspec_getfil   
;     Version 1.1
;
; PURPOSE:
;   Pass back a structure and/or filename for a specified file type in
;   the HIRES pipeline.
;
; CALLING SEQUENCE:
;   
;  rslt = hamspec_getfil('type', [setup], /NAME, SUBFIL=, CHKFIL=,
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
;   ordr_str = hamspec_getfil('ordr_str', 1)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2011 Written by JXP
;-
;------------------------------------------------------------------------------

function hamspec_getfil, type, setup, NAME=name, SUBFIL=subfil, $
                      CHKFIL=chkfil, SZ=sz, INDX=indx, FIL_NM=fil_name, $
                      HEAD=head, FRAME=frame


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = hamspec_getfil(type, [setup], /NAME, SUBFIL=, CHKFIL=,'
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
     stop
     ; Should set binning from keywords
  endif

  ;;  Assume hamspec files have this type of file name:  hamspec0000.fits 

  if keyword_set( SUBFIL ) then begin
      dot = rstrpos( SUBFIL[0], '.')
      
      expname = strmid(subfil,usepos,dot-usepos)
      mname  = 'm' + expname
      mfits  = 'm' + expname + '.fits'
   endif


  ;;
  if keyword_set( FRAME ) then begin
      dum = frame
      cframe = ''
      while( DUM LT 1000 ) do begin
          dum = dum*10
          cframe = cframe+'0'
      endwhile
      cframe = cframe+strtrim(frame,2)
  endif

  ;; Big case statement
  case type of
      'arc_fil': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Arcs/Arc_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'arc_fit': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Arcs/Fits/Arc_'+c_s+'_fit.idl'
          NAME = 1
      end
      'arc_2Dfit': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Arcs/Fits/Arc_'+c_s+'_fit2D.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_trc': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Arcs/TRC/Arc_'+c_s+'_T.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_fittrc': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Arcs/TRC/Arc_'+c_s+'_F.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_mkaimg': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Arcs/Arc_'+c_s+'_I.fits'
      end
      ;; Bias
      'bias_fil': begin
          if not keyword_set( SZ ) then stop
          fil_name = 'Bias/Bias'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'ov_fil': begin
          if not keyword_set( FRAME ) then stop
          fil_name = 'OV/ov_'+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Final file
      'fin_fil': begin
          if not keyword_set( FRAME ) then stop
          fil_name = 'Final/f_hamspec'+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'fspec_fil': begin
          if not keyword_set( SUBFIL ) then stop
          fil_name = 'FSpec/'+strtrim(subfil)+'.fits' 
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Obj file
      'obj_fil': begin
          if not keyword_set( FRAME ) then stop
          fil_name = 'Extract/Obj_'+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'objN_fil': begin
          if not keyword_set( FRAME ) then stop
          fil_name = 'Extract/Obj_'+cframe+'N.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'sens_fil': begin
          if not keyword_set( FRAME ) then stop
          fil_name = 'Extract/sens_'+cframe+'.idl'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Sky file
      'sky_fil': begin
          if not keyword_set( FRAME ) then stop
          fil_name = 'Sky/sky_hamspec'+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'scatt_fil': begin
          if not keyword_set( FRAME ) then stop
          fil_name = 'Sky/scatt_hamspec'+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; FLATS
      'qtz_fil': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Flats/TFlat_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'nqtz_fil': begin  ;; Normalized flat
          if not keyword_set( SETUP ) then stop
          fil_name = 'Flats/TFlat_N'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'pixflt_fil': begin  ;; Normalized flat
          if not keyword_set( SETUP ) then stop
          fil_name = 'Flats/PFlat_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'npixflt_fil': begin  ;; Normalized flat
          if not keyword_set( SETUP ) then stop
          fil_name = 'Flats/PFlat_N'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'mflt_fil': begin  ;; Milky flat, already normalized (yippee)
          if not keyword_set( SETUP ) then stop
          fil_name = 'Flats/Flat_'+c_s+'_M.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'tflat_str': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Flats/TStr_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Order structure
      'ordr_str': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'Flats/OStr_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Logs
      'flat_log': begin
          NAME = 1
          fil_name = 'Logs/flat_redux.log'
      end
      ;; QA
      'qa_trcflat': begin
          NAME = 1
          fil_name = 'QA/Flats'+c_s+'/qa_trcflt_'+c_s+'.ps'
      end
      'qa_slitflat': begin
          NAME = 1
          fil_name = 'QA/Flats'+c_s+'/qa_slitflat_'+c_s+'.ps'
      end
      'qa_fittflat': begin
          NAME = 1
          fil_name = 'QA/Flats'+c_s+'/qa_fittflt_'+c_s+'.ps'
      end
      'qa_arcfit': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'QA/Arcs'+c_s+'/qa_arcfit_'+c_s+'.ps'
          NAME = 1
      end
      'qa_arc2dfit': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'QA/Arcs'+c_s+'/qa_arc2dfit_'+c_s+'.ps'
          NAME = 1
      end
      'qa_tracearc': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'QA/Arcs'+c_s+'/qa_tracearc_'+c_s+'.ps'
          NAME = 1
      end
      'qa_fittrcarc': begin
          if not keyword_set( SETUP ) then stop
          fil_name = 'QA/Arcs'+c_s+'/qa_fitrcarc_'+c_s+'.ps'
          NAME = 1
      end
      'qa_fntobj': begin
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_name = 'QA/Obj'+c_s+'/qa_fntobj_'+cframe+'.ps'
      end
      'qa_skysub': begin
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_name = 'QA/Obj'+c_s+'/qa_skysub_'+cframe+'.ps'
      end
      'qa_extract': begin
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_name = 'QA/Obj'+c_s+'/qa_extract_'+cframe+'.ps'
      end
      'qa_chkflux': begin
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_name = 'QA/Obj'+c_s+'/qa_flux_'+cframe+'.ps'
      end
      else: stop
  endcase

  ;; CHKFIL
  if arg_present( CHKFIL ) then a = findfile(fil_name+'*', count=CHKFIL)

  ;; Return
  if keyword_set( NAME ) then return, fil_name $
  else begin
      ;; Check for file
      if not keyword_set( CHKFIL ) then a = findfile(fil_name+'*', count=CHKFIL)
      if CHKFIL NE 0 then return, xmrdfits(fil_name, indx, head, /silent) $
      else begin
          case type of
              'arc_trc': $
                print, 'File ', fil_name, ' does not exist. Run hamspec_tracearc!'
              'arc_2Dfit': $
                print, 'File', fil_name, ' does not exist. Run hamspec_fit2darc!'
              'arc_fittrc': $
                print, 'File', fil_name, ' does not exist. Run hamspec_fittrcarc!'
              'flat_fil': $
                print, 'File', fil_name, ' does not exist. Run hamspec_mktflat!'
              else: print, 'hamspec_getfil: File does not exist! ', fil_name
          endcase
          stop
      endelse
      return, -1
  endelse

  return, -1
end

