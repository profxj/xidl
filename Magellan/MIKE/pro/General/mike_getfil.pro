;+ 
; NAME:
; mike_getfil   
;     Version 1.1
;
; PURPOSE:
;   Pass back a structure and/or filename for a specified file type in
;   the MIKE code.
;
; CALLING SEQUENCE:
;   
;  rslt = mike_getfil('type', [setup], SIDE=, /NAME, SUBFIL=, CHKFIL=,
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
;   SIDE    -- Specify camera: Blue (1) or Red (2)
;   INDX    -- Image extension in the fits file (generally 0, 1, or 2)
;
; OPTIONAL OUTPUTS:
;   FIL_NM  -- Filename of the file 
;   HEAD    -- Image header
;
; COMMENTS:
;
; EXAMPLES:
;   ordr_str = mike_getfil('ordr_str', 1, SIDE=1)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

function mike_getfil, type, setup, SIDE=side, NAME=name, SUBFIL=subfil, $
                      CHKFIL=chkfil, SZ=sz, INDX=indx, FIL_NM=fil_nm, $
                      HEAD=head


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = mike_getfil(type, [setup], SIDE=, /NAME, SUBFIL=, CHKFIL=,'
      print, '     SZ=, INDX=, FIL_NM=, HEAD=) [v1.1]'
      return, -1
  endif 

  ;; Setup
  if keyword_set( SETUP ) then begin
      if setup GE 10 then c_s = strtrim(setup,2) $
      else c_s = '0'+strtrim(setup,2) 
  endif


  ;; Size
  if keyword_set( SZ ) then begin
      cbin = 2048L / sz[0]
      rbin = 4096L / sz[1]
  endif

  ;;  Assume mike files have this type of file name:  aa0000.fits or a0000.fits
  ;;  Assume the last letter is b or r and designates blue or red

  if keyword_set( SUBFIL ) then begin
      dot = rstrpos( SUBFIL[0], '.')
      rpos   = rstrpos( SUBFIL[0], 'r', dot)
      bpos   = rstrpos( SUBFIL[0], 'b', dot)
      
      if keyword_set(side) then begin
          if side EQ 1 AND bpos EQ -1 then $
            print, "This doesn't look like a blue MIKE file"
          
          if side EQ 2 AND rpos EQ -1 then $
            print, "This doesn't look like a red MIKE file"
      endif
      
      usepos = rpos > bpos      ; max of these two
      if bpos EQ rpos OR usepos EQ -1 then begin
          print,  "Unlikely to have found correct name, help!"
          return, ''
      endif
      
      ;; Deal with 'm'
      if strmid(subfil[0], (usepos-1)>0,1) EQ 'm' then flg_m = 1 else flg_m = 0

      if NOT keyword_set(side) then begin
          if rpos GE 0 AND rpos GT bpos then side = 2L
          if bpos GE 0 AND bpos GT rpos then side = 1L
      endif
      
      if dot LE usepos then begin
          print,  "Can't find the file exposure name, help!"
          return, ''
      endif
      
      expname = strmid(subfil,usepos,dot-usepos)
      mname  = 'm' + expname
      mfits  = 'm' + expname + '.fits'
      midl  = 'm' + expname + '.idl'
  endif

  ;; Side
  if keyword_set( SIDE ) then begin
      if side EQ 1 then nm = 'B' else nm = 'R' 
  endif

  ;; Big case statement
  case type of
      'arc_fil': begin
          fil_nm = 'Arcs/Arc_'+mfits
          if not keyword_set( INDX ) then indx = 0L
      end
      'arc_fit': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Arcs/Fits/'+ mname +'_fit.idl'
          NAME = 1
      end
      'arc_img': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Arcs/Arc_' + mname + 'I.fits'
      end
      'arc_fittrc': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Arcs/TRC/Arc_'+ mname +'_F.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_trc': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Arcs/TRC/Arc_'+mname+'_T.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_2Dfit': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Arcs/Fits/Arc_'+mname+'_fit2D.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'guess_arc': begin
          if not keyword_set( SZ ) then stop
          if not keyword_set( SIDE ) then stop
          ;; CCD change
          if (side EQ 1) AND $
            (sxpar(head,'UT-DATE') LT '2009-03-15') then $
            pref = 'old_' else pref = ''
          fil_nm = getenv('MIKE_DIR')+'/pro/Arcs/'+pref+'templ_arc_'+$
                   strtrim(cbin,2)+'x'+$
                   strtrim(rbin,2)+nm+'.idl' 
      end
      ;; Bias
      'bias_fil': begin
          if not keyword_set( SZ ) then stop
          if not keyword_set( SIDE ) then stop
          fil_nm = 'Bias/Bias'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+nm+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'ov_fil': begin
          if not keyword_set( SUBFIL ) then stop
          if flg_m EQ 1 then fil_nm = 'OV/ov_'+mname+'.fits' $
          else fil_nm = 'OV/ov_'+expname+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Final file
      'fin_fil': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Final/f_'+mfits
          if not keyword_set( INDX ) then indx = 0L
      end
      'fspec_fil': begin
          if not keyword_set( SUBFIL ) then stop
          if not keyword_set( SIDE ) then stop
          fil_nm = 'FSpec/'+strtrim(subfil)+nm+'.fits' 
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Obj file
      'obj_fil': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Extract/Obj_'+mfits
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Sky file
      'sky_fil': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Sky/sky_'+mfits
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Sky file
      'scatt_fil': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Sky/scatt_'+mfits
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Sky file
      'sens_fil': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Extract/Sens_'+midl
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; FLATS
      'mflat_fil': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/Flat_'+nm+'_'+c_s+'_M.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'tflat_fil': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/Flat_'+nm+'_'+c_s+'_T.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'tflat_str': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/TStr_'+nm+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Order structure
      'ordr_str': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/OStr_'+nm+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; QA
      'qa_trcflat': begin
          NAME = 1
          fil_nm = 'QA/Flats'+c_s+'/qa_trcflt_'+c_s+nm+'.ps'
      end
      'qa_slitflat': begin
          NAME = 1
          fil_nm = 'QA/Flats'+c_s+'/qa_slitflat_'+c_s+nm+'.ps'
      end
      'qa_fittflat': begin
          NAME = 1
          fil_nm = 'QA/Flats'+c_s+'/qa_fittflt_'+c_s+nm+'.ps'
      end
      'qa_arcfit': begin
          if not keyword_set( SUBFIL ) then stop
          NAME = 1
          fil_nm = 'QA/Arcs'+c_s+'/qa_arcfit_'+mname+'.ps'
      end
      'qa_arc2dfit': begin
          if not keyword_set( SUBFIL ) then stop
          NAME = 1
          fil_nm = 'QA/Arcs'+c_s+'/qa_arc2dfit_'+mname+'.ps'
      end
      'qa_tracearc': begin
          if not keyword_set( SUBFIL ) then stop
          NAME = 1
          fil_nm = 'QA/Arcs'+c_s+'/qa_tracearc_'+mname+'.ps'
      end
      'qa_fittrcarc': begin
          if not keyword_set( SUBFIL ) then stop
          NAME = 1
          fil_nm = 'QA/Arcs'+c_s+'/qa_fittrcarc_'+mname+'.ps'
      end
      'qa_fntobj': begin
          if not keyword_set( SUBFIL ) then stop
          NAME = 1
          fil_nm = 'QA/Obj'+c_s+'/qa_fntobj_'+mname+'.ps'
      end
      'qa_skysub': begin
          if not keyword_set( SUBFIL ) then stop
          NAME = 1
          fil_nm = 'QA/Obj'+c_s+'/qa_skysub_'+mname+'.ps'
      end
      else: stop
  endcase

  ;; CHKFIL
  if arg_present( CHKFIL ) then a = findfile(fil_nm+'*', count=CHKFIL)

  ;; Return
  if keyword_set( NAME ) then return, fil_nm $
  else begin
      ;; Check for file
      if not keyword_set( CHKFIL ) then a = findfile(fil_nm+'*', count=CHKFIL)
      if CHKFIL NE 0 then return, xmrdfits(fil_nm, indx, head, /silent) $
      else begin
          case type of
              'arc_trc': $
                print, 'File ', fil_nm, ' does not exist. Run mike_tracearc!'
              'arc_2Dfit': $
                print, 'File', fil_nm, ' does not exist. Run mike_fit2darc!'
              'arc_fittrc': $
                print, 'File', fil_nm, ' does not exist. Run mike_fittrcarc!'
              'mflat_fil': $
                print, 'File', fil_nm, ' does not exist. Run mike_mkmflat!'
              'tflat_fil': $
                print, 'File', fil_nm, ' does not exist. Run mike_mktflat!'
              else: print, 'mike_getfil: File does not exist! ', fil_nm
          endcase
          stop
      endelse
      return, -1
  endelse

  return, -1
end

