;+ 
; NAME:
; uves_getfil   
;     Version 1.1
;
; PURPOSE:
;   Pass back a structure and/or filename for a specified file type in
;   the MIKE code.
;
; CALLING SEQUENCE:
;   
;  rslt = uves_getfil('type', [setup], WCEN=, /NAME, SUBFIL=, CHKFIL=,
;  SZ=, INDX=, FIL_NM=, HEAD=)
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
;   SZ      -- Image size
;   SUBFIL  -- Image name generally used to parse the root name of the
;              image  (e.g.  'Arcs/arc_mb0539.fits').  Required in
;              many cases.
;   WCEN    -- Specify camera: Blue (1) or Red (2)
;   INDX    -- Image extension in the fits file (generally 0, 1, or 2)
;
; OPTIONAL OUTPUTS:
;   FIL_NM  -- Filename of the file 
;   HEAD    -- Image header
;
; COMMENTS:
;
; EXAMPLES:
;   ordr_str = uves_getfil('ordr_str', 1, WCEN=340)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------

function uves_getfil, type, setup, WCEN=wcen, NAME=name, SUBFIL=subfil, $
                      CHKFIL=chkfil, SZ=sz, INDX=indx, FIL_NM=fil_cwcen, $
                      HEAD=head, FRAME=frame, OBJN=objn


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = uves_getfil(type, [setup], WCEN=, /NAME, SUBFIL=, CHKFIL=,'
      print, '     SZ=, INDX=, FIL_NM=, HEAD=, FRAME=) [v1.1]'
      return, -1
  endif 

  ;; Setup
  if keyword_set( SETUP ) then begin
      if setup GE 10 then c_s = strtrim(setup,2) $
      else c_s = '0'+strtrim(setup,2) 
  endif


  ;; Size
;  if keyword_set( SZ ) then begin
;      cbin = round(2048. / sz[0])
;      rbin = round(4096. / sz[1])
;  endif
  

  ;;  Assume uves files have this type of file name:  uves0000.fits 

  if keyword_set( SUBFIL ) then begin
      dot = rstrpos( SUBFIL[0], '.')
      
      expname = strmid(subfil,usepos,dot-usepos)
      mname  = 'm' + expname
      mfits  = 'm' + expname + '.fits'
  endif

  ;; WCEN
  if keyword_set( WCEN ) then begin
      cwcen = strtrim(round(wcen),2)
  endif

  if keyword_set( FRAME ) then begin
      dum = frame
      cframe = ''
      while( DUM LT 10 ) do begin
          dum = dum*10
          cframe = cframe+'0'
      endwhile
      cframe = cframe+strtrim(frame,2)
  endif

  ;; Big case statement
  case type of
      'arc_fil': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Arcs/Arc_'+cwcen+'_'+c_s+'_'+cframe+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'arc_img': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Arcs/Arc_'+cwcen+'_'+c_s+'_'+cframe+'I.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'arc_fit': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Arcs/Fits/Arc_'+cwcen+'_'+c_s+'_'+cframe+'_fit.idl'
          NAME = 1
      end
      'arc_2Dfit': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Arcs/Fits/Arc_'+cwcen+'_'+c_s+'_'+cframe+'_fit2D.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_trc': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Arcs/TRC/Arc_'+cwcen+'_'+c_s+'_'+cframe+'_T.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_fittrc': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Arcs/TRC/Arc_'+cwcen+'_'+c_s+'_'+cframe+'_F.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_mkaimg': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Arcs/Arc_'+cwcen+'_'+c_s+'_'+cframe+'_I.fits'
      end
      ;; Bias
;      'bias_fil': begin
;          if not keyword_set( SZ ) then stop
;          if not keyword_set( WCEN ) then stop
;          fil_cwcen = 'Bias/Bias'+strtrim(cbin,2)+ $
;            'x'+strtrim(rbin,2)+cwcen+'.fits'
;          if not keyword_set( INDX ) then indx = 0L
;      end
      'ov_fil': begin
;          if not keyword_set( WCEN ) then stop
          if not keyword_set( OBJN ) then stop
          fil_cwcen = 'OV/ov_'+objn
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Final file
      'fin_fil': begin
          if not keyword_set( OBJN ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Final/f_'+objn+'_'+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'fspec_fil': begin
          if not keyword_set( SUBFIL ) then stop
          if not keyword_set( WCEN ) then stop
          fil_cwcen = 'FSpec/'+strtrim(subfil)+cwcen+'.fits' 
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Obj file
      'obj_fil': begin
          if not keyword_set( OBJN ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Extract/Obj_'+objn+'_'+cframe+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'objN_fil': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Extract/Obj_'+cframe+cwcen+'N.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'sens_fil': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Extract/sens_'+cframe+cwcen+'.idl'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Sky file
      'sky_fil': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Sky/sky_uves_'+cframe+'_'+cwcen+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'scatt_fil': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'Sky/scatt_'+cframe+'_'+cwcen+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; FLATS
      'qtz_fil': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cwcen = 'Flats/TFlat_'+cwcen+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'nqtz_fil': begin  ;; Normalized flat
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cwcen = 'Flats/TFlat_N'+cwcen+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'pixflt_fil': begin  ;; Normalized flat
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cwcen = 'Flats/PFlat_N'+cwcen+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'tflat_str': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cwcen = 'Flats/TStr_'+cwcen+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Order structure
      'ordr_str': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          fil_cwcen = 'Flats/OStr_'+cwcen+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Logs
      'flat_log': begin
          NAME = 1
          fil_cwcen = 'Logs/flat_redux.log'
      end
      ;; QA
      'qa_trcflat': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          NAME = 1
          fil_cwcen = 'QA/Flats'+cwcen+'_'+c_s+'/qa_trcflt_'+cwcen+'_'+c_s+'.ps'
      end
      'qa_slitflat': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          NAME = 1
          fil_cwcen = 'QA/Flats'+cwcen+'_'+c_s+'/qa_slitflat'+cwcen+'_'+c_s+'.ps'
      end
      'qa_fittflat': begin
          NAME = 1
          fil_cwcen = 'QA/Flats'+cwcen+'_'+c_s+'/qa_fittflt'+cwcen+'_'+c_s+'.ps'
      end
      'qa_arcfit': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'QA/Arcs'+cwcen+'_'+c_s+'/qa_arcfit_'+cwcen+'_'+c_s+'_'+cframe+'.ps'
          NAME = 1
      end
      'qa_arc2dfit': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'QA/Arcs'+cwcen+'_'+c_s+'/qa_arc2dfit_'+cwcen+'_'+c_s+'_'+cframe+'.ps'
          NAME = 1
      end
      'qa_tracearc': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'QA/Arcs'+cwcen+'_'+c_s+'/qa_tracearc_'+cwcen+'_'+c_s+'_'+cframe+'.ps'
          NAME = 1
      end
      'qa_fittrcarc': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( FRAME ) then stop
          fil_cwcen = 'QA/Arcs'+cwcen+'_'+c_s+'/qa_fitrcarc_'+cwcen+'_'+c_s+'_'+cframe+'.ps'
          NAME = 1
      end
      'qa_fntobj': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( SETUP ) then stop
          if not keyword_set( OBJN ) then stop
          fil_cwcen = 'QA/Obj'+cwcen+'_'+c_s+'/qa_fntobj'+objn
          NAME = 1
      end
      'qa_skysub': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_cwcen = 'QA/Obj'+cwcen+'_'+c_s+'/' + $
                      'qa_skysub_'+cframe+'_'+cwcen+'.ps'
      end
      'qa_extract': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_cwcen = 'QA/Obj'+c_s+'/qa_extract_'+cframe+cwcen+'.ps'
      end
      'qa_chkflux': begin
          if not keyword_set( WCEN ) then stop
          if not keyword_set( FRAME ) then stop
          NAME = 1
          fil_cwcen = 'QA/Obj'+c_s+'/qa_flux_'+cframe+cwcen+'.ps'
      end
      else: stop
  endcase

  ;; CHKFIL
  if arg_present( CHKFIL ) then a = findfile(fil_cwcen+'*', count=CHKFIL)

  ;; Return
  if keyword_set( NAME ) then return, fil_cwcen $
  else begin
      ;; Check for file
      if not keyword_set( CHKFIL ) then a = findfile(fil_cwcen+'*', count=CHKFIL)
      if CHKFIL NE 0 then return, xmrdfits(fil_cwcen, indx, head, /silent) $
      else begin
          case type of
              'arc_trc': $
                print, 'File ', fil_cwcen, ' does not exist. Run uves_tracearc!'
              'arc_2Dfit': $
                print, 'File', fil_cwcen, ' does not exist. Run uves_fit2darc!'
              'arc_fittrc': $
                print, 'File', fil_cwcen, ' does not exist. Run uves_fittrcarc!'
              'flat_fil': $
                print, 'File', fil_cwcen, ' does not exist. Run uves_mktflat!'
              else: print, 'uves_getfil: File does not exist! ', fil_cwcen
          endcase
          stop
      endelse
      return, -1
  endelse

  return, -1
end

