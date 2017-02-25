;+ 
; NAME:
; imacsls_getfil   
;     Version 1.0
;
; PURPOSE:
;   Pass back a structure (or filename) given simple input
;
; CALLING SEQUENCE:
;   
;  rslt = imacsls_getfil(type, setup, side)
;
; INPUTS:
;   type   --  Type of file
;   setup  --  Setup ID 
;
; RETURNS:
;
; OUTPUTS:
;  Structure, image, name, etc.
;
; OPTIONAL KEYWORDS:
;   /name   - Only return resolved name (string)
;   SIDE=   - Side (CCD to the blue=1 or red=2)
;   SUBFIL= - Sub file requried to grab the right file
;
; OPTIONAL OUTPUTS:
;   FIL_NM  - Name of file
;   CHKFIL  - Value equal to the number of files matching name
;
; COMMENTS:
;
; EXAMPLES:
;   ordr_str = imacsls_getfil('ordr_str', 1, 1)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------
function imacsls_getfil, type, setup, SIDE=side, NAME=name, SUBFIL=subfil, $
                      CHKFIL=chkfil, SZ=sz, INDX=indx, FIL_NM=fil_nm, $
                      HEAD=head


  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'rslt = imacsls_getfil(type, [setup], SIDE=, /NAME, SUBFIL=, CHKFIL=,'
      print, '     SZ=) [v1.0]'
      return, -1
  endif 

  ;; Setup
  if keyword_set( SETUP ) then begin
      if setup GE 10 then c_s = strtrim(setup,2) $
      else c_s = '0'+strtrim(setup,2) 
  endif

  ;; Side
  if keyword_set( SIDE ) then begin
      if side EQ 1 then nm = 'B' else nm = 'R' 
  endif

  ;; Size
  if keyword_set( SZ ) then begin
      cbin = 2048L / sz[0]
      rbin = 4096L / sz[1]
  endif

  ;; Big case statement
  case type of
      'arc_fil': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Arcs/Arc_'+nm+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      'arc_fit': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, '.fits')
          fil_nm = 'Arcs/Fits/'+strmid(subfil,5,spos-5)+'_fit.idl'
          ;; IDL file (name only)
          NAME = 1
      end
      'arc_img': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, '.fits')
          fil_nm = strmid(subfil,0,spos)+'I.fits'
      end
      'arc_fittrc': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, '.fits')
          fil_nm = 'Arcs/TRC/'+strmid(subfil,5,spos-5)+'_F.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_trc': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, '.fits')
          fil_nm = 'Arcs/TRC/'+strmid(subfil,5,spos-5)+'_T.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'arc_2Dfit': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, '.fits')
          fil_nm = 'Arcs/Fits/'+strmid(subfil,5,spos-5)+'_fit2D.fits'
          if not keyword_set( INDX ) then indx = 1L
      end
      'guess_arc': begin
          if not keyword_set( SZ ) then stop
          if not keyword_set( SIDE ) then stop
          fil_nm = getenv('MIKE_DIR')+'/pro/Arcs/templ_arc_'+$
            strtrim(cbin,2)+'x'+$
            strtrim(rbin,2)+nm+'.idl' 
      end
      'arc_psfil': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, '.fits')
          fil_nm = 'Arcs/Fits/'+strmid(subfil,5,spos-5)+'fit.ps'
          if not keyword_set( INDX ) then indx = 1L
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
          spos = strpos(subfil, 'm')
          fil_nm = 'OV/ov_'+strmid(subfil, spos)
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Final file
      'fin_fil': begin
          if not keyword_set( SUBFIL ) then stop
          spos = strpos(subfil, 'm')
          fil_nm = 'Final/f_'+strmid(subfil, spos)
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
          fil_nm = 'Extract/Obj_'+strmid(subfil,spos)
          if not keyword_set( INDX ) then indx = 1L
      end
      ;; FLATS
      'flat_fil': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/Flat_'+nm+'_'+c_s+'.fits' 
          if not keyword_set( INDX ) then indx = 0L
      end
      ;; Order structure
      'ordr_str': begin
          if not keyword_set( SIDE ) then stop
          if not keyword_set( SETUP ) then stop
          fil_nm = 'Flats/OStr_'+nm+'_'+c_s+'.fits'
          if not keyword_set( INDX ) then indx = 1L
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
                print, 'File ', fil_nm, ' does not exist. Run imacsls_tracearc!'
              'arc_2Dfit': $
                print, 'File', fil_nm, ' does not exist. Run imacsls_fit2darc!'
              'arc_fittrc': $
                print, 'File', fil_nm, ' does not exist. Run imacsls_fittrcarc!'
              'mflat_fil': $
                print, 'File', fil_nm, ' does not exist. Run imacsls_mkmflat!'
              'tflat_fil': $
                print, 'File', fil_nm, ' does not exist. Run imacsls_mktflat!'
              else: print, 'imacsls_getfil: File does not exist! ', fil_nm
          endcase
          stop
      endelse
      return, -1
  endelse

  return, -1
end

