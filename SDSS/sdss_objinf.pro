;+ 
; NAME:
; sdss_objinf   
;    Version 1.0
;
; PURPOSE:
;    Given an array of strings (QSO names), return info and plot 
;
; CALLING SEQUENCE:
;   
; 
;
; INPUTS:
;   strings - Array of strings
;
; RETURNS:
;   uniq  - Array of unique members
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   sort - Sort the output
;
; OPTIONAL OUTPUTS:
;  COUNT - number of unique strings
;
; COMMENTS:
;
; EXAMPLES:
;   flg = x_chkfil( lbls, count=count)
;
;sdss_objinf, '0571-52286-276', /plot
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_objinf, names, SUMMF=SUMMF, PLOT=plot, DATDIR=datdir, OUTFIL=outfil

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'sdss_objinf, names, /PLOT (v1.0)'
      return
  endif 
  
  if not keyword_set( SUMMF ) then $
    summf = '/net/profx/data4/SDSS/DR1_QSO/summ_qso_dr1_R195.fit'
  if not keyword_set( DATDIR ) then $
    datdir = '/net/profx/data4/SDSS/DR1_QSO/spectro/1d_20/'

;;  Summary list
  sdss = xmrdfits(summf, 1, /silent)

;; Loop
  nnm = n_elements(names)
  print, 'Name', 'zobj', 'Mag(R)', 'RA(2000)', 'DEC(2000)', $
    FORMAT='(4x,a4,8x,a4,2x,a6,3x,a8,4x,a9)'
  for q=0L,nnm-1 do begin

      ;; Parse the name
      pnm = strmid(names[q],0,4)
      plate = long( pnm )
      mjd = long( strmid(names[q],5,5) )
      fnm = strmid(names[q],11,3)
      fid = long( fnm )

      ;; Find the ID
      iobj = where(sdss.plate EQ plate AND sdss.mjd EQ mjd and $
                   sdss.fiberid EQ fid, nobj)
      case nobj of
          0: begin
              flg = 0
              print, 'sdss_objinf: Target '+names[q]+' was not found!'
              print, 'sdss_objinf: continuing..'
          end
          1: flg=1 
          else: stop
      endcase
      
      if flg EQ 0 then continue
          

      ;; Print to screen
      x_radec, ra, dec, sdss[iobj].raobj, sdss[iobj].decobj, /flip
      print, names[q], ra, dec, sdss[iobj].mag_r,  sdss[iobj].z, $
        FORMAT='(a14,6x,a11,1x,a9,2x, f5.2,f5.2)'

      ;; Plot
      if keyword_set( PLOT ) then begin
          ;; File name
          sfil = datdir+pnm+'/1d/'+'spSpec-'+$
            strtrim(mjd,2)+'-'+pnm+'-'+fnm+'.fit'
          ;; Plot
          x_specplot, sfil, inflg=5, /block
      endif
          
  endfor
      
  

return

end
