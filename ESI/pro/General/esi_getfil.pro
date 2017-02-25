;+ 
; NAME:
; esi_getfil   
;     Version 1.0
;
; PURPOSE:
;   Pass back a structure (or filename) given simple inputs
;
; CALLING SEQUENCE:
;   
;  rslt = esi_getfil('file', setup, side)
;
; INPUTS:
;   setup     -  
;   [side]    -  
;
; RETURNS:
;
; OUTPUTS:
;  Structure, image, name, etc.
;
; OPTIONAL KEYWORDS:
;   /name   - Only return resolved name (string)
;   CHKFIL  - Value equal to the number of files matching name
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ordr_str = esi_getfil('ordr_str', 1, 1)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

function esi_getfil, type, mode, SIDE=side, NAME=name, SUBFIL=subfil, $
                      CHKFIL=chkfil, SZ=sz, INDX=indx, FIL_NM=fil_nm, $
                      HEAD=head, CBIN=cbin, RBIN=rbin, SLIT=slit, ORDR=ordr, $
                      OBJID = OBJID


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = esi_getfil(type, [setup], SIDE=, /NAME, SUBFIL=, CHKFIL=,'
      print, '     SZ=) [v1.0]'
      return, -1
  endif 

  if keyword_set( SLIT ) then c_s = esi_slitnm( slit )

  ;; Size
  if keyword_set( SZ ) then begin
      cbin = 2048L / sz[0]
      rbin = 4096L / sz[1]
  endif

  ;; Big case statement
  case type of
      'arc_fil': begin
          fil_nm = 'Arcs/ArcECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'arc_trc': begin
          if ordr LT 10 then cordr = '0'+string(ordr, FORMAT='(i1)') $
          else cordr = string(ordr, FORMAT='(i2)')
          fil_nm = 'Arcs/TRC/ArcECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'trc'+cordr+'.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_trcps': begin
          if ordr LT 10 then cordr = '0'+string(ordr, FORMAT='(i1)') $
          else cordr = string(ordr, FORMAT='(i2)')
          fil_nm = 'Arcs/TRC/ArcECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'trc'+cordr+'.ps'
          NAME=1
      end
      'arc_fit': begin
          fil_nm = 'Arcs/ArcECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'fit.idl'
          ;; IDL file (name only)
          NAME = 1
      end
      'arc_2Dfit': begin
          fil_nm = 'Arcs/ArcECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'2Df.fits'
          ;; IDL file (name only)
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_psfil': begin
          fil_nm = 'Arcs/ArcECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'fit.ps'
          NAME=1
      end
      'arc_img': begin
          fil_nm = 'Arcs/ArcECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'IMG.fits'
          indx = 0L
      end
;      'arc_fittrc': begin
;          if not keyword_set( SUBFIL ) then stop
;          spos = strpos(subfil, '.fits')
;          fil_nm = 'Arcs/TRC/'+strmid(subfil,5,spos-5)+'_F.fits'
;          if not keyword_set( INDX ) then indx = 1L
;      end
;      'guess_arc': begin
;          if not keyword_set( SZ ) then stop
;          if not keyword_set( SIDE ) then stop
;          fil_nm = getenv('MIKE_DIR')+'/pro/Arcs/templ_arc_'+$
;            strtrim(cbin,2)+'x'+$
;            strtrim(rbin,2)+nm+'.idl' 
;      end
      ;; Bias
      'bias_fil': begin
          if n_elements(MODE) EQ 0 THEN STOP
;          if not keyword_set( MODE ) then stop
          if mode EQ 0 then b_s = 'I' else b_s = 'S'
          fil_nm = 'Bias/Bias'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+b_s+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; OV
      'ov_fil': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, 'm')
          fil_nm = 'OV/ov_'+strmid(subfil, spos)
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Final file
      'fin_fil': begin
          if not keyword_set( SUBFIL ) then stop
          fil_nm = 'Final/f_'+subfil
          if not keyword_set( INDX ) then indx = 0L
       end
      ;; Fring file
      'fringe_fil': begin
         if not keyword_set( SUBFIL ) then stop
         fil_nm = 'Final/fringe_'+subfil
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
          spos = strpos(subfil, 'm')
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; Standards
      'std_trc': begin
          fil_nm = 'Extract/STD_ECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'TRC.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'sens_fil': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, '.fits')
          fil_nm = 'Extract/sens_'+strmid(subfil,0,spos)+'.idl'
          if not keyword_set( INDX ) then indx = 1L
      end
      'sens_fil_save': BEGIN
           if not keyword_set( SUBFIL ) then stop
           spos = strpos(subfil, '.fits')
           fil_nm = 'Extract/sens_saveall'+strmid(subfil, 0, spos)+'.idl'
           if not keyword_set( INDX ) then indx = 1L
       end
      ;; FLATS
      'echflat_fil': begin
         IF mode EQ 0 OR mode EQ 1 THEN BEGIN 
            IF mode EQ 0 THEN f_s = '_I' $
            ELSE IF mode EQ 1 THEN f_s = '_D'
            fil_nm = 'Flats/FlatECH'+c_s+'_'+strtrim(cbin, 2)+ $
                     'x'+strtrim(rbin, 2)+f_s+'.fits'
         ENDIF ELSE IF mode EQ 2 THEN $
            fil_nm = 'Flats/TwiflatECH'+c_s+'_'+strtrim(cbin, 2)+ $
                     'x'+strtrim(rbin, 2)+'.fits'
         if not keyword_set( INDX ) then indx = 0L
      end
      'finflat_fil': begin
          fil_nm = 'Flats/FlatECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'N.fits'
          if not keyword_set( INDX ) then indx = 0L
       end
      'twiflat_fil': begin
         fil_nm = 'Flats/TwiflatECH'+c_s+'_'+strtrim(cbin, 2)+ $
                  'x'+strtrim(rbin, 2)+'N.fits'
         if not keyword_set( INDX ) then indx = 0L
      end
      'sedg_fil': begin
          fil_nm = 'Flats/SEdgECH'+c_s+'_'+strtrim(cbin,2)+ $
            'x'+strtrim(rbin,2)+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      'trc_flat': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/TStr_'+nm+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Order structure
      'ordr_str': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/OStr_'+nm+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
       end
      ;; Profile file
;      'prof_fil': begin
;         if not keyword_set( SUBFIL ) then stop
;         spos = strpos(subfil, '.fits')
;         fil_nm = 'Profile/prof_'+strmid(subfil, 0, spos) + $
;                  '-' + objid + '.fits'
;         fil_nm = 'Profile/prof_'+subfil
;         if not keyword_set( INDX ) then indx = 1L
;      end
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
                print, 'File ', fil_nm, ' does not exist. Run esi_tracearc!'
              'arc_2Dfit': $
                print, 'File', fil_nm, ' does not exist. Run esi_fit2darc!'
              'arc_fittrc': $
                print, 'File', fil_nm, ' does not exist. Run esi_fittrcarc!'
              'mflat_fil': $
                print, 'File', fil_nm, ' does not exist. Run esi_mkmflat!'
              'tflat_fil': $
                print, 'File', fil_nm, ' does not exist. Run esi_mktflat!'
              else: print, 'esi_getfil: File does not exist! ', fil_nm
          endcase
          stop
      endelse
      return, -1
  endelse

  return, -1
end

