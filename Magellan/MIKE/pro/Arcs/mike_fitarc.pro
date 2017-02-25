;+ 
; NAME:
; mike_fitarc   
;     Version 1.1
;
; PURPOSE:
;  To identify and centroid arc lines in each order.  There is
;  actually no fitting done at this stage other than to reject bad
;  lines.   The main program calls mike_fitarc_work as its driver.
;  The algorithm mike_fitarc_work does the following (via x_fitarc):
;    1) Input the arc image from mike_mkarc 
;    2) Input an archived arc solution appropriate to the side
;    3) Extract 1D (boxcar) spectra down each order :: extract_arc
;    4) Cross correlate (FFT) against the archived 1D arc spectrum,
;    this gives the order number and the pixel offset
;    5) Automatically identify a set of lines (x_templarc)
;    6) Perform a low order fit to these lines
;    7) Reidentify all lines with this fit and refit 
;    8) Write arc solutions (one per order) to a fits file 
;    9) If the orders extend beyond the archived solution, attempt to
;    extrapolate to the remaining orders.  The idea is to use the
;    known wavelengths from the good orders to extrapolate a guess at
;    the solution.  Then use this guess to peak up on the arc lines in
;    these additional orders.
;
; CALLING SEQUENCE:
;   
;  mike_fitarc, mike, setup, obj_id, [side], /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, IORDR=, /PINTER 
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   SHFTPRM=  - Fit structure for shifting the orders of the arc
;   /PINTER   - Perform fit for pre-identified lines using x_identify
;   /INTER    - Identify lines interactively and then fit
;   LINLIST   -  Arc line list (default: $XIDL_DIR/Spec/Arcs/Lists/hires_thar.lst
;   /CHK      - Manually check steps along the way
;   /DEBUG    - Debugging
;   SIGREJ=   - Rejection sigma for outliers in arc line fitting
;              (default: 2.)
;   IORDR     - Initial order for analysis
;   FORDR     - Final order for analysis
;   /CLOBBER  - Overwrite previous fits
;   SHFTPRM=  - Fit structure for shifting the orders of the arc
;   /FWEIGHT  - Centroid arc lines using a flux-weighting scheme
;   IPSIG     - Array which sets the sigma significance for arc lines
;               as a function of order#.  This is for the initial
;               solution and it is recommended to use a higher value
;               than IFSIG.  An example:  [ [0,40,10],
;               [40,80,20]] -- This sets the significance to 10 for
;                              orders 0 to 40 and 20 for orders 40 to
;                              80.
;   IFSIG     - Similar to IPSIG except this is for the final line
;               list.  One may tend to use weaker lines (lower
;               significance) for this step. 
;   /NOLOG    - Indicates the template file does not have Log
;               wavelengths
;   SATUR     - Saturation level for the Arc [default: 30000.]
;   /NOECHO   - Throw out orders 77 and redder for the blue side to
;               deal with the echo in the old dichroic.  (JXP 9/05)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_fitarc, mike, 1, 1
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_fitarc
;   mike_crossarc
;   mike_fitarc_work
;   extract_arc
;   x_templarc
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   21-Aug-2002 Streamlined + Added ps output
;   01-Feb-2003 Polished (JXP)
;   01-Jan-2004 Added a guess at unmatched orders (SB)
;   Apr-2005 Ported to XIDL via x_fitarc (JXP)
;-
;------------------------------------------------------------------------------

function mike_fitarc_work, arc_fil, setup, side, $
                      INTER=inter, LINLIST=linlist, NOECHO=noecho, $
                      CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
                      IORDR=iordr, PINTER=pinter, IPSIG=ipsig, IFSIG=ifsig, $
                      GUESSARC=guessarc, SHFTPRM=shftprm, _EXTRA=extra

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'a = mike_fitarc_work(arc_fil, setup, side, /INTER, IORDR=, /DEBUG, /CHK '
      print, '     SHFTPRM=, /clobber, /NOECHO) [v2.0]'
      return, -1
  endif 

;  Optional Keywords
  if not keyword_set( SIGREJ ) then sigrej = 2.
  if not keyword_set( TRCNSIG ) then trcnsig = 3.
  if not keyword_set( RADIUS ) then radius = 2.0
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/mike_thar.lst' 
;    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' 
  if not keyword_set( PKWDTH ) then pkwdth = 3L
  if not keyword_set( NOLOG ) then FLOG = 1 else FLOG = 0
;  if not keyword_set( XOFF ) then xoff = 0.

  print, 'mike_fitarc: Using line list ', linlist

  ;; Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side, fil_nm=ordr_fil)
  nordr = n_elements(ordr_str)
  temp_ordr = ordr_str
  if keyword_set( SHFTPRM ) then $  ; Shift as requesetd
    xoff = mike_shifti(shftprm, OSTR=temp_ordr)

  ;; OUTFIL
  out_fil = mike_getfil('arc_fit', setup, SIDE=side, $
                        /name, SUBFIL=arc_fil, CHKFIL=chkfil)
  if keyword_set(IORDR) AND CHKFIL NE 0 $
    AND not keyword_set( CLOBBER ) then begin
      print, 'mike_fitarc: Arc fits file exists!  ' + $
        'Use /CLOBBER to overwrite'
      return, out_fil
  endif 

  ;; Set orders
  if side EQ 1 then begin
      ostep = 1 
      BCKWD=1
  endif else ostep = -1

  ;; PSIG
  if side EQ 1 then begin
      if not keyword_set(IPSIG) then ipsig = [ [0,105,30], $
                                               [105,109,20], $
                                               [109,200,5] ]
      if not keyword_set(IFSIG) then ifsig = [ [0,103,20], $
                                               [103,200,5] ]
  endif else begin
      if not keyword_set(IPSIG) then ipsig = [ [0,40,10], $
                                               [40,200,30] ]
      if not keyword_set(IFSIG) then ifsig = [ [0,45,10], $
                                               [45,200,30] ]
  endelse

  qafil = mike_getfil('qa_arcfit', setup, SUBFIL=arc_fil) 

  ;; Open archived solutions
  ;; Need size
  arc = xmrdfits(arc_fil, 0, head, /silent)
  sz_arc = size(arc, /dimensions)
  delvarx, arc

  if not keyword_set( GUESSARC ) then begin
      guessarc = mike_getfil('guess_arc', SIDE=side, SZ=sz_arc, /name, $
                             chkfil=chkfil, HEAD=head)
      if chkfil EQ 0 then $
        guessarc = mike_getfil('guess_arc', SIDE=side, SZ=[1024,2048], /name, $
                               chkfil=chkfil, HEAD=head)

      if CHKFIL EQ 0 then begin
          print, "mike_fitarc: Cannot find arc template: " , guessarc
          stop
      endif
  endif
  print, 'mike_fitarc: Using this template file: ', guessarc


  rslt = x_fitarc( arc_fil, ordr_str, out_fil, $
                   INTER=inter, LINLIST=linlist,$ 
                   CHK=chk, CLOBBER=clobber, IPSIG=ipsig, IFSIG=ifsig, $
                   PINTER=pinter, DEBUG=debug, SIGREJ=sigrej, $
                   IORDR=iordr, GUESSARC=guessarc, OSTEP=ostep, $
                   QAFIL=qafil, NOLOG=nolog, SATUR=30000., $
                   FIT_MSK=fit_msk, BCKWD=bckwd, _EXTRA=extra)

  ;;
  ordr_str.flg_anly = 1

  ;; Zero out 77 and redder for the original data
  ;; These have an 'echo' in the arc lines (each line has a pair)
  head = xheadfits(arc_fil)
  if ((side EQ 1) AND (sxpar(head,'UT-DATE') LE '2004-05-14')) OR $
    keyword_set(NOECHO) then begin
      restore, out_fil
      mt = where(guess_ordr LT 78L, nmt)
      if nmt NE 0 then begin
          print, 'mike_fitarc:  Zeroing out the echo in order #77 and redder !!'
          sv_lines[mt].nlin = 0
          save, guess_ordr, sv_aspec, all_arcfit, sv_lines, rejstr, $
                filename=out_fil
      endif
      ;; Zero out the orders too
      bad = where(ordr_str.order LE 77, nbad)
      if nbad NE 0 then ordr_str[bad].flg_anly = 0
  endif
      
  ;; Write order structure
  mwrfits, ordr_str, ordr_fil, /create, /silent

  return, arc_fil
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_fitarc, mike, setup, obj_id, side, $
                 INTER=inter, LINLIST=linlist, $
                 CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
                 IORDR=iordr, PINTER=pinter


;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_fitarc, mike, setup, obj_id, [side], /INTER, IORDR=, /DEBUG,'+ $
        '/CHK, SIGREJ=, /CLOBBER, LINLIST= [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]

; Loop on side
  for kk=0L,n_elements(side)-1 do begin
      qq = side[kk]
      ;; SIDE
      if qq EQ 1 then print, 'mike_fitarc: Fitting BLUE arc' $
      else print, 'mike_fitarc: Fitting RED arc'

      ;; Grab all obj indices
      indx = where(mike.flg_anly NE 0 AND mike.side EQ qq AND $
                   mike.obj_id EQ obj_id AND mike.setup EQ setup AND $
                   (strtrim(mike.type,2) EQ 'OBJ' OR $
                   strtrim(mike.type,2) EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'mike_fitarc: No Obj found!  Skipping..' 
          continue
      endif

      arcfil_all= mike[indx].arc_fil
      asort = sort(arcfil_all)
      auniq = uniq(arcfil_all[asort])

      ;; LOOP
      for mm=0L,n_elements(auniq)-1 do begin
 
          arcfil = arcfil_all[asort[auniq[mm]]]
          idx = indx[asort[auniq[mm]]]

          ;; Check for arcfil
          if x_chkfil(arcfil+'*',/silent) EQ 0 then begin
              print, 'mike_fitarc: Use mike_procarc to create arc file first!!'
              continue
          endif

          rslt = mike_fitarc_work( arcfil, setup, qq, INTER=inter, $
                                   LINLIST=linlist,$ 
                                   CHK=chk, CLOBBER=clobber, $
                                   PINTER=pinter, DEBUG=debug, SIGREJ=sigrej, $
                                   IORDR=iordr)
;                                   XOFF=mike[idx].arc_xyoff[0], $ 
          if size(rslt,/tname) NE 'STRING' then stop
      endfor
  endfor

; All done
  print, 'mike_fitarc: All done!'

  return
end
