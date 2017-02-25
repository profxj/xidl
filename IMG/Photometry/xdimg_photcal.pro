;+ 
; NAME:
; xdimg_photcal   
;   Version 1.1
;
; PURPOSE:
;    Drives the program x_intphotcal which allows the user to 
;      interactively create photometric solutions.
;
; CALLING SEQUENCE:
;   
;   xdimg_photcal, struct, LNDTFIL=, SETAM=, XSIZE=, YSIZE=, MIN_NOBS=
;     MIN_MOBS=, SETAM=, OUTROOT=, MAGROOT=, /NCLR
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest.  This
;             program focuses on the STD frames only.
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  LNDTFIL  = Filename of landolt fits file [default:
;        $XIDL_DIR/IMG/Photometry/Lists/nlandolt.fits ]
;  SETAM    = AM terms (one per filter: UBVRI)
;  MIN_MOBS = Minimum number of nights observed by Landolt (default=2)
;  MIN_NOBS = Minimum number of observations by Landolt (default=4)
;  XSIZE=   - Size of gui in screen x-pixels (default = 1000)
;  YSIZE=   - Size of gui in screen y-pixels (default = 600)
;  OUTROOT  - Root name of Twilight flats (default is 'Photo/soln_')
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_photcal, dimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_photcal, struct, LNDTFIL=lndtfil, XSIZE=xsize, YSIZE=ysize, $
                   MIN_NOBS=min_nobs, MIN_MOBS=min_mobs, SETAM=setam, $
                   OUTROOT=outroot, MAGROOT=magroot, NCLR=nclr

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'xdimg_photcal, struct, LNDTFIL=, XSIZE=, YSIZE=, MIN_NOBS= '
      print, '      MIN_MOBS=, SETAM=, OUTROOT=, /NCLR, MAGROOT= (v1.1)'
      return
  endif 

  close, /all

;  Optional Keywords

  if not keyword_set( LNDTFIL ) then $
    lndtfil = getenv('XIDL_DIR')+'/IMG/Photometry/Lists/nlandolt.fits'

  if n_elements( MIN_NOBS ) EQ 0 then min_nobs = 4
  if n_elements( MIN_MOBS ) EQ 0 then min_mobs = 2
  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200

  if not keyword_set( OUTROOT ) then outroot = 'Photo/soln_'
  if not keyword_set( MAGROOT ) then magroot = 'Photo/mag_'
  outroot = strtrim(outroot,2)

;  Temp arrarys

  tmp_am = fltarr(10000)
  tmp_mag = fltarr(10000)
  tmp_sig = fltarr(10000)
  tmp_nam = strarr(10000)
  tfilt = ' '
  tlist = ' '
  tnam = ' '
  tmp = {stdstruct}

; Read in Landolt

  landolt = mrdfits(lndtfil, 1, STRUCTYP='lndltstr')
  landolt.Name = strtrim(landolt.Name,2)

;  Find the Standard Stars

  stds = where(struct.type EQ 'STD' AND struct.flg_anly NE 0 AND $
               struct.flg_final NE 0, nstds)
  if nstds EQ 0 then begin
      print, 'xdimg_photcal: No standard stars to analyse!'
      return
  endif

;  Loop on Filter

  x_filters, struct[stds].filter, filt, nfilt

  for q=0,nfilt-1 do begin
      
      wfilt = where(strtrim(struct[stds].filter, 2) EQ filt[q], numi)

      
;  Read in the Mag data

      nstrs = 0
      for w=0,numi-1 do begin
          jj = stds[wfilt[w]]
          magfil = strjoin([strtrim(magroot,2), 'f_', $
                            strmid(strtrim(struct[jj].img_root,2),0, $
                                   strlen(strtrim(struct[jj].img_root,2))-5), '.dat'])
          nlin = numlines(magfil)
          

          openr, 1, magfil
          readf, 1, FORMAT='(a)', tfilt
          readf, 1, FORMAT='(a)', tlist
          readf, 1, tAM

          for i=0,nlin-4 do begin
              readf, 1, FORMAT='(a15,f,f)', tnam, tmg, tsg
              tmp_AM[nstrs] = tAM
              tmp_nam[nstrs] = strtrim(tnam, 2)
              tmp_mag[nstrs] = tmg 
              tmp_sig[nstrs] = tsg
              nstrs = nstrs + 1
          endfor
          close, 1
      endfor
      
      obs = replicate(tmp, nstrs)

      obs.flg_anly = 1
      obs.Name = tmp_nam[0:nstrs-1]
      obs.Mag = tmp_mag[0:nstrs-1]
      obs.filter = filt[q]
      obs.sig_Mag = tmp_sig[0:nstrs-1]
      obs.AM = tmp_AM[0:nstrs-1]

      ;  Set Color term

      if nfilt GT 1 then begin
          if q LT nfilt-1 then obs.sCLR = strjoin([filt[q],filt[q+1]]) $
          else obs.sCLR = strjoin([filt[q-1],filt[q]]) 
      endif else obs.sCLR = 'NC'

      ; AM term
      if keyword_set( SETAM ) then int_setam = setam[q] else int_setam=0

      outfil = strjoin([outroot, filt[q], '_', obs[0].sCLR, '.dat'])
      ; NCLR
      if not keyword_set(NCLR) then begin
          if nfilt GT 1 then nclr=0 else nclr=1
      endif

      x_intphotcal, obs, landolt, outfil, XSIZE=xsize, YSIZE=ysize,  $
        MIN_NOBS=min_nobs, MIN_MOBS=min_mobs, SETAM=int_setam, NCLR=nclr

      print, 'Output is in ', outfil

      delvarx, ans, sig
  endfor

  delvarx, obs, tmp_am, tmp_mag, tmp_sig, tmp_nam, tfilt, tlist, tnam

end
