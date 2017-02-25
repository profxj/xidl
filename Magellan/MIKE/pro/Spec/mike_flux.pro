;+ 
; NAME:
; mike_flux
;     Version 1.1
;
; PURPOSE:
;    Fluxes the spectrum using a standard star calibration.  Default
;    is to use the one in CALIBS which is probably good enough for
;    relative fluxing.  Definitely not good enough for absolute.
;
; CALLING SEQUENCE:
;   
;  mike_flux, mike, setup, obj_id, side, FLUXFIL=, /CLOBBER, /STD
;
; INPUTS:
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   side    -  Blue (1) or Red (2) side
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLUXFIL -  File name of standard star for fluxing (default:
;               $MIKE_DIR/pro/Std/Archive/sens_blue#.fits)
;   /STD    - Extraction should be set for a standard star
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_flux, mike, obj_id, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2004 JXP
;   09-Jan-2009 Moustakas - boxcar keyword added
;-
;------------------------------------------------------------------------------

pro mike_flux, mike, setup, obj_id, side, exp, FLUXFIL=fluxfil, $
               CLOBBER=clobber, STD=std, OLD=old, boxcar=boxcar

  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_flux, mike, setup, obj_id, side, [exp], FLUXFIL=, ' 
      print, '/CLOBBER, /STD, /OLD     [v1.1]'
      return
  endif

;  Optional Keywords


; Grab exposure
  if not keyword_set( STD ) then begin
      allexp = where((mike.type EQ 'OBJ' OR mike.type EQ 'STD') $
                     AND mike.flg_anly NE 0 AND $
                     mike.setup EQ setup AND mike.obj_id EQ obj_id $
                     AND mike.side EQ side, nexp)
      if nexp EQ 0 then begin
          print, 'mike_flux: No objects found!'
          return
      endif
      ;;  Exposures
      if n_elements(exp) NE 0 then allexp = allexp[exp]
      nexp = n_elements(allexp)
  endif else begin
      allexp = obj_id[0]
      nexp = 1
  endelse

  if not keyword_set( FLUXFIL ) then begin
      stop  ;; Not recommended to use the defaults right now (JXP: 6/26/05)
      ;; But it still works ok!  (JXP: 1/24/14)
      if side EQ 1 then begin
          fluxfil = getenv('MIKE_DIR')+'/pro/Std/Archive/sens_blue'+ $
            strtrim(mike[allexp[0]].rowbin,2)+'.fits'
      endif else begin
          fluxfil = getenv('MIKE_DIR')+'/pro/Std/Archive/sens_red'+ $
            strtrim(mike[allexp[0]].rowbin,2)+'.fits'
      endelse
      if x_chkfil(fluxfil) EQ 0 then begin
          print, 'mike_flux:  Archived file', fluxfil, 'does not exist!'
          print, 'mike_flux:  You will need to create your own..'
          return
      endif else print, 'mike_flux:  Using archived solution ', fluxfil
  endif

  if not keyword_set(old) then begin
      restore, fluxfil  
;      if n_elements(tot_fit) NE n_elements(objstr) then stop
  endif

; Flux func
  print, 'mike_flux: Fluxing...'
  for qq=0,nexp-1 do begin
      idx = allexp[qq]
      ;; Open objfil
      objfil = mike[idx].obj_fil
      objstr = xmrdfits(objfil, 1, /silent)

     if keyword_set(old) then begin
;        flux = objstr.fx * 0
;        sig  = objstr.var * 0

         ;; Expects bsplines
         isv = 0L
         for jj=0L,n_elements(objstr)-1 do begin
            if objstr[jj].npix EQ 0 then continue
             print, 'mike_flux: Order ', objstr[jj].order
             gpx = lindgen(objstr[jj].npix)
             ;; Read bset
             flg_flx = 1
             for ii=isv,99 do begin
                 bset = xmrdfits(fluxfil, ii+1, /silent)
                 if size(bset,/type) NE 8 then flg_flx = 0
                 if size(bset,/type) NE 8 then break
                 isv = ii
                 if bset.ordr EQ objstr[jj].order then break
             endfor
             if flg_flx EQ 1 then begin
                 ;; Get values
                 full_rtio = bspline_valu(objstr[jj].wave[gpx], bset)
                 
                 objstr[jj].flux[gpx] = objstr[jj].fx[gpx]*full_rtio/mike[idx].exp
                 b = where(objstr[jj].var[gpx] GT 0.)
                 objstr[jj].sig[gpx[b]] = sqrt(objstr[jj].var[gpx[b]]) $
                   * full_rtio[b] / mike[idx].exp
                 objstr[jj].flg_flux = 1
             endif else begin
                 print, 'mike_flux:  Unable to flux this order ', $
                   objstr[jj].order
                 isv = 0L
             endelse
         endfor
         
      endif else begin
          ;; Legendre polynomials
          for jj=0L,n_elements(objstr)-1 do begin
              ;; Boxcar
              if keyword_set(BOXCAR) then begin
                  objstr[jj].npix = objstr[jj].nrow
                  objstr[jj].wave = objstr[jj].box_wv
                  objstr[jj].fx = objstr[jj].box_fx
                  objstr[jj].var = objstr[jj].box_var
              endif
              ;; Match up
              mtch = where(objstr[jj].order EQ ordr_fit,nmt)
              if nmt EQ 0 then begin
                  print, 'mike_flux: No match to order = ', objstr[jj].order
                  continue
              endif
              mtch = mtch[0]
              if (objstr[jj].npix gt 0L) then begin ; jm09jan05nyu
                 gpx = lindgen(objstr[jj].npix)
                 objstr[jj].flux = 0.0 ; zero-out
                 objstr[jj].sig = 0.0  ; zero-out
                 ;; Sensitivity function
;                 gds = where(sens_wv[*,mtch] GT 0.)
;                 linterp, sens_wv[gds,mtch], sens_func[gds,mtch], $
;                   objstr[jj].wave[gpx], full_rtio
                 ;; Apply
                 full_rtio = x_calcfit(objstr[jj].wave[gpx], fitstr=tot_fit[mtch])
                 objstr[jj].flux[gpx] = objstr[jj].fx[gpx]/full_rtio/mike[idx].exp
                 b = where(objstr[jj].var[gpx] GT 0.)
                 objstr[jj].sig[gpx[b]] = sqrt(objstr[jj].var[gpx[b]]) $
                   / full_rtio[b] / mike[idx].exp
                 objstr[jj].flg_flux = 1
              endif
           endfor
      endelse
 
      print, 'mike_flux: Writing fluxed file ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil
  endfor
  
  print, 'mike_flux: All done'
  return

end

