;+ 
; NAME:
; mike_chkarctrc   
;     Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;  mike_fitarc, mike, slit, /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, ORDRS=, PIXSHFT=, /PINTER
;
; INPUTS:
;   mike     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   /PINTER   - Perform fit only for the pre-identified lines
;              (recommended)
;   /INTER    - Identify lines interactively and then fit
;   LINLIST=  -  Arc line list (default:
;                                $XIDL_DIR/Spec/Arcs/Lists/ESI_.lst)
;   GUESSARC= -  Initial guess for wavelength solution (default:
;                                $XIDL_DIR/ESI/CALIBS/ESI_arcfit.idl)
;   /CHK      - Manually check steps along the way
;   /DEBUG    - More intensive checking
;   ORDRS=    - Only fit select orders (default: [0L,9L] )
;   PIXSHFT=  - Manually set pixel shift from calib file (deafult:
;                Let the program determine this using FFT formalism)
;   SIGREJ=   - Rejection sigma for outliers in arc line fitting
;              (default: 2.)
;   /CLOBBER  - Overwrite previous fits
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   There is some fine tuning in a couple of orders to help with the
;    fits.  Only setup for 1x1 binning right now
;
; EXAMPLES:
;   mike_fitarc, mike, 0.5, /CHK, /PINTER, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Apr-2003 Written by SB
;
;  Usage:
;   mike_trcflat, flatfil, trc_set
;   mike_fittflat, trc_set, oset
;   mike_chkarctrc, arcfile, oset, aset 
;
;  WARNING: Gives a segmentation fault, but I have yet to trap
;      the fault.  It must be in a call_external
;
;
pro mike_chkarctrc, mike, indx, NOSTOP=nostop
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_chkarctrc, mike, indx, /NOSTOP [v1.0]'
      return
  endif 

;  Optional Keywords

; Setup
  setup = mike[indx].setup
  if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 

; Loop on side
      ;; SIDE
  side = mike[indx].side
  case side of 
      1:  begin
          print, 'mike_chkarctrc: Tracing BLUE arc'
          nm = 'B'
          nord = 2
      end
      2: begin
          print, 'mike_chkarctrc: Tracing RED arc'
          nm = 'R' 
          nord = 1
      end
      else: stop
  endcase

   ;; 
  spos = strpos(mike[indx].img_root, '.fits')
  xset_fil = 'Arcs/TRC/Arc_m'+strmid(mike[indx].img_root,0,spos)+'_T.fits'
  
  ;; Warn about xatv!
  arc_fil = 'Arcs/Arc_'+mike[indx].img_root
  if not keyword_set( NOSTOP ) then begin
      print, 'mike_chkarctrc:  You need to run this command before '+ $
        'running this program:  '
      print, '  IDL>  xatv, ''', arc_fil, ''''
      print, 'If you havent done this yet, return out of this program'
      print, '  run the command and rerun this program with /NOSTOP'
      print, 'If you have done this command, just continue'
      stop
  endif
  
  xatverase
  ;; Restore
;  restore, filnm

  ;; Open order structure
  ordr_fil = 'Flats/OStr_'+nm+'_'+c_s+'.fits'
  ordr_str = xmrdfits(ordr_fil,1)
 
  ;; XFIL
;  xset_fil = 'Arcs/TRC/Arc_'+nm+'_'+c_s+'T.fits'
  xset = xmrdfits(xset_fil, 1, /silent)

  ;; Loop on order
  for q=0L, n_elements(ordr_str)-1 do begin

      ;; Fit for slope variation
      fitstr = x_setfitstrct(NORD=nord, FLGREJ=1, HSIG=2.5, LSIG=2.5)
      fit = x_fitrej(xset[q].xguess[0:xset[q].ngood-1], $
                     xset[q].coeff[1,0:xset[q].ngood-1], FITSTR=fitstr)
      if fit[0] NE -1 then begin
          print, strtrim(q,2), ': RMS = ', fitstr.rms/fit[0]
          
          ;; Loop
          nlin = xset[q].ngood
          for i=0L,nlin-1 do begin
              ;; Slope and b value
              m = x_calcfit(xset[q].xguess[i], FITSTR=fitstr) / xset[q].xmax
;              m = 2.*x_calcfit(xset[q].xguess[i], FITSTR=fitstr) / xset[q].xmax
              x0 = 0.5*(ordr_str[q].lhedg[xset[q].xguess[i]] + $
                        ordr_str[q].rhedg[xset[q].xguess[i]]) $
                + mike[indx].arc_xyoff[0]
              b = xset[q].xguess[i] - m*x0
              
              ;; Plot
              xval = ordr_str[q].lhedg[xset[q].xguess[i]] + $
                findgen(round(ordr_str[q].rhedg[xset[q].xguess[i]]- $
                              ordr_str[q].lhedg[xset[q].xguess[i]])+1) + $
                mike[indx].arc_xyoff[0]
              xatvplot, xval, xval*m + b
          endfor
      endif
  endfor

  return
end
