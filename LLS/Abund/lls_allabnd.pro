;+ 
; NAME:
; lls_allabd
;  V1.1
;
; PURPOSE:
;    Given a LLS structure, calculate abundances for various
;    transitions using the AODM.  Also fill up elemental abundances
;    using dominant ions.
;
; CALLING SEQUENCE:
;   lls_allabd, lls, nn
;
; INPUTS:
;   lls -- DLA structure
;   nn  -- Index of the structure
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;  /NOCLM -- Dont bother to do column densities.  Useful for
;            plotting only
;  /SKIP -- Skip systems without a datafile provided
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lls_allabd, lls
;
;
; PROCEDURES CALLED:
;  upd_fd
;
; REVISION HISTORY:
;   08-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function lls_updfd, val, flg

  if val EQ 0 then return, flg
  
  if (val MOD 2*flg) GT (flg-1) then return, val $
  else return, val + flg
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lls_allabnd, lls, nn, sys, FLG_PLT=flg_plt, TYP_FLG=typ_flg, ROOT=root, $
                 NOCLM=noclm, LUN=lun, FLGNN=flgnn, SKIP=skip

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lls_allabd, lls, [nn],[sys], FLG_PLT=, ABNDFIL= [v1.1]'
    return
  endif 

  if not keyword_set(ROOT) then root = ''

  ;; Fine-structure lines are a pain..
  if not keyword_set( FINEDAT ) then $
    finedat = getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'
  readcol, finedat, fin_Z, fin_ion, fin_j, fin_wav, fin_E, fin_A, $
           FORMAT='I,I,F,D,F,F', /silent
  c2wv = fin_wav[where(fin_Z EQ 6 and fin_ion EQ 2)]
  si2wv = fin_wav[where(fin_Z EQ 14 and fin_ion EQ 2)]
  o1wv = fin_wav[where(fin_Z EQ 8 and fin_ion EQ 1)]

  if not keyword_set(FLGNN) then flgnn = 0
;  if size(nn,/type) EQ 0 then begin
;;      sys=0
;      lls = {llsstruct}
;      nn = 0L
;      if not keyword_set(ABNDFIL) then return
;      lls[nn].systems[sys].abndfil = abndfil
;      flgnn = 1
;  endif

  if not keyword_set( LLS_SIG ) then lls_sig = 3.
;  if not keyword_set( CONTI_ERR ) then conti_err = 0.
;  resolve, 'ew_to_colm'

  ;; Initialize
  lls[nn].systems[sys].ion.state.flgclm = -1
  lls[nn].systems[sys].ion.state[*,0].flgclm = 0

  dumc = ''
  dumf = 0.0
  dumf2 = 0.0
  dumf3 = 0.0
  dumd = 0.d
  dumd2 = 0.d
  dumi = 0

  ;; Open data file
  openr, lun, root+strtrim(lls[nn].systems[sys].abndfil,2), /get_lun
  readf, lun, dumc
  llsnm = dumc
  ;; Data file(s)
  readf, lun, dumi
  lls[nn].systems[sys].ndfil = dumi

  ;; HIRES
  if lls[nn].systems[sys].ndfil MOD 2 GT 0 then begin
      readf, lun, dumc
      lls[nn].systems[sys].Hfil = dumc
  endif
  ;; ESI 
  if lls[nn].systems[sys].ndfil MOD 4 GT 1 then begin
      readf, lun, dumc
      lls[nn].systems[sys].Efil = dumc
  endif
  ;; UVES
  if lls[nn].systems[sys].ndfil MOD 8 GT 3 then begin
      readf, lun, dumc
      lls[nn].systems[sys].Ufil = dumc
  endif
  ;; X?
  if lls[nn].systems[sys].ndfil MOD 16 GT 7 then begin
      readf, lun, dumc
      lls[nn].systems[sys].Xfil = dumc
  endif
  ;; MIKEB
  if lls[nn].systems[sys].ndfil MOD 32 GT 15 then begin
      readf, lun, dumc
      lls[nn].systems[sys].MBfil = dumc
  endif
  ;; MIKER
  if lls[nn].systems[sys].ndfil MOD 64 GT 31 then begin
      readf, lun, dumc
      lls[nn].systems[sys].MRfil = dumc
  endif
  ;; Redshift
  readf, lun, dumd
  if dumd NE lls[nn].systems[sys].zabs and flgnn EQ 0 then begin
      print, 'lls_allabnd: Update zabs fool!'
      print, dumd, lls[nn].systems[sys].zabs
      stop
  endif else lls[nn].systems[sys].zabs = dumd

  ;; Tables
  readf, lun, dumc
  lls[nn].systems[sys].tab_fil = strtrim(dumc,2)

  ;; HI
  readf, lun, dumf, dumd2
  if dumf NE lls[nn].systems[sys].NHI and flgnn EQ 0 then begin
      print, 'lls_allabnd: Update NHI fool!', lls[nn].systems[sys].NHI, dumf
      wait, 1
  endif ;else lls[nn].systems[sys].NHI = dumf

  lls[nn].systems[sys].ion[1].indx[1] = 1
  lls[nn].systems[sys].ion[1].state[1].flgclm = 1
  lls[nn].systems[sys].ion[1].state[1].clm = dumd
  lls[nn].systems[sys].ion[1].state[1].lambda = 1215.6701d
  lls[nn].systems[sys].ion[1].state[1].sigclm = dumd2
  lls[nn].systems[sys].ion[1].state[1,1].flgclm = 1
  lls[nn].systems[sys].ion[1].state[1,1].clm = dumd
  lls[nn].systems[sys].ion[1].state[1,1].lambda = 1215.6701d
  lls[nn].systems[sys].ion[1].state[1,1].sigclm = dumd2
  lls[nn].systems[sys].elm[1].clm = dumd
  lls[nn].systems[sys].elm[1].sigclm = dumd2

  ;; Extras?
  nhand = 0
  readf, lun, nhand
  for i=0L,nhand-1 do begin
      Znmb = 0
      readf, lun, Znmb
      cin = 0.
      ein = 0.
      fins = 0
      readf, lun, cin, ein, fins
      lls[nn].systems[sys].elm[Znmb].flginst = fins
      ein = 10^(cin+ein) - 10^cin
      ;; Calc the column (and fill)
      dla_calcclm, lls[nn].systems, sys, ZNmb, 10^cin, ein
      lls[nn].systems[sys].elm[Znmb].flgclm = 1
  endfor

  if keyword_set(NOCLM) then begin
     ;; Remember to close the file!
     return
  endif

  while not eof(lun) do begin
      ;; Ion flag
      readf, lun, dumi
      ionflg = dumi
      ;; Ion info
      readf, lun, dumd, dumf, dumf2, dumi;, dumf3
      lambda = dumd
      vmin = dumf
      vmax = dumf2
      flg_fil = dumi
;      contierr = dumf3
      contierr = 0.

      ;; Get atomic info
      getion, lambda, ion, elm, Z=atnmb

      ;; CII*
      if min(abs(c2wv-lambda)) LT 1e-3 then ion = 6
;      if lambda EQ 1335.7077d then ion = 6

      ;; SiII*
      if min(abs(si2wv-lambda)) LT 1e-3 then ion = 9

      ;; OI* and OI**
      if min(abs(o1wv-lambda)) LT 1e-3 then ion = 9

      ;; Increment + basic info
      idx = lls[nn].systems[sys].ion[atnmb].indx[ion]++ + 1
      lls[nn].systems[sys].ion[atnmb].state[ion,idx].lambda = lambda
      lls[nn].systems[sys].ion[atnmb].state[ion,idx].flginst = flg_fil
      lls[nn].systems[sys].ion[atnmb].state[ion,idx].vmnx[0] = vmin
      lls[nn].systems[sys].ion[atnmb].state[ion,idx].vmnx[1] = vmax

      ;; Deal with by hand columns
      if ionflg GE 8 then begin
         if ionflg LT 16 then begin
            lls[nn].systems[sys].ion[atnmb].state[ion,idx].clm = 10^vmin
            lls[nn].systems[sys].ion[atnmb].state[ion,idx].sigclm = $
               (10^(vmin+vmax) - 10^(vmin-vmax))/2	      
            lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = ionflg
         endif  ;; If GE 16, then this is an EW!
         continue
      endif

      ;; Read in data
      ;print, 'lambda: ', lambda
      ;stop
      if not keyword_set( SV_FIL ) then sv_fil = -1
      mdflg = flg_fil mod 64
      if flg_fil GT 128 then iflg=2 else iflg=0
      if mdflg NE SV_FIL then begin
          case mdflg of
              1: dat = x_readspec(lls[nn].systems[sys].Hfil, /struct, /autofsig, NOSTOP=SKIP)
              2: dat = x_readspec(lls[nn].systems[sys].Efil, /struct, /autofsig, NOSTOP=SKIP)
              4: dat = x_readspec(lls[nn].systems[sys].Ufil, /struct, /autofsig, INFLG =iflg, NOSTOP=SKIP)
              8: dat = x_readspec(lls[nn].systems[sys].Xfil, /struct, /autofsig, NOSTOP=SKIP, INFLG=iflg)
              16: dat = x_readspec(lls[nn].systems[sys].MBfil, /struct, /autofsig,INFLG=iflg, NOSTOP=SKIP)
              32: dat = x_readspec(lls[nn].systems[sys].MRfil, /struct, /autofsig, NOSTOP=SKIP)
              else: stop
           endcase
      endif else sv_fil = mdflg
      if (size(dat))[0] EQ 0 then break

      ;; Velocity array
      if vmin GT vmax then begin
          print, 'lls_allabnd: Fix vmin, vmax', vmin, vmax
          stop
      endif
      dumv = (x_allvelo(dat.wv, lls[nn].systems[sys].zabs, lambda, $
                            [vmin, vmax], ALL_PMNX=pmnx))[*]
      dat.velo[pmnx[0]:pmnx[1]] = dumv[0:pmnx[2,0]]

      ;; AODM
      if (pmnx[1]-pmnx[0]) LE 1 then stop

      if n_elements(dat.sig) NE n_elements(dat.wv) then begin
            dat.sig = dat.fx
            gd = where(dat.sig GT 0)
            dat.sig[gd] = sqrt(dat.sig[gd])
      endif

      x_aodm, dat.wv, dat.fx[pmnx[0]:pmnx[1]], $
        dat.sig[pmnx[0]:pmnx[1]], lambda, clm, sigc, $
        VELO=dumv[0:pmnx[2,0]], FLG_SAT=flg_sat
      
      ;; Continuum error
      if CONTIERR GT 0. then begin
          dw = dat.wv[pmnx[1]] - dat.wv[pmnx[0]]
          sig_conti = contierr * dw * 1000. ;; mA
          ;; Convert to column
          sig_ci = ew_to_colm([lambda], [sig_conti])
          ;; Add in quadrature
          sigc = sqrt(sigc^2 + sig_ci[0]^2)
      endif
      lls[nn].systems[sys].ion[atnmb].state[ion,idx].clm = clm
      lls[nn].systems[sys].ion[atnmb].state[ion,idx].sigclm = sigc

      ;;;;;;;
      ;; Flag

      ;; CII*
      if lambda EQ 1335.7077d then begin
          ;; Saturation
          if flg_sat GT 0 then begin
              lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = 2
              continue
          endif
          if clm LT lls_sig*sigc then begin
              lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = 4
              lls[nn].systems[sys].ion[atnmb].state[ion,idx].clm = lls_sig*sigc
              continue
          endif
          if ionflg GE 4  then begin
              lls[nn].systems[sys].ion[atnmb].state[ion,idx].clm = clm > lls_sig*sigc
              lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = ionflg
              continue
          endif
      endif

      ;; Standard
      if flg_sat GT 0 OR ionflg EQ 2 OR ionflg EQ 3  then begin
          ;; Lower limits
          if (ionflg MOD 8) gt 1 then begin ; Keep flg given by hand
              lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = ionflg
              if keyword_set(VERBOSE) then print, 'lls_allabnd: Taking value in the file'
          endif else lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = ionflg + 2
      endif else begin

          ;; Upper limit?  ion_flg = 4,5
          if clm GT lls_sig*sigc then begin
              lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = ionflg
              ;; Blends
              if ionflg EQ 4 OR ionflg EQ 5 then $
                lls[nn].systems[sys].ion[atnmb].state[ion,idx].clm = lls_sig*sigc > clm 
          endif else begin
              ;; Upper limit!!
              if ionflg LT 4 then begin
                  lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = ionflg + 4
                  lls[nn].systems[sys].ion[atnmb].state[ion,idx].clm = lls_sig*sigc 
              endif else begin
                  lls[nn].systems[sys].ion[atnmb].state[ion,idx].flgclm = ionflg 
                  lls[nn].systems[sys].ion[atnmb].state[ion,idx].clm = $
                    lls_sig*sigc > clm
              endelse
          endelse
      endelse
  endwhile

  ;;ccccccccccccccccccccccccccccccccccccc
  ;;  Update elemental column if appropriate
  ;;ccccccccccccccccccccccccccccccccccccc

  ;; Loop on Z
  imx = 100 < ((size(lls.systems[sys].elm))[1]-1)
  for i=2,imx do begin
      ;; Loop on ion
      for j=1,6 do begin
          ionflg = 0
          ;; Deal with by hand!
          if lls[nn].systems[sys].elm[i].flgclm EQ 1 then ionflg=1
          ;; CII*
          if i EQ 6 AND j EQ 6 then continue
          tidx = lls[nn].systems[sys].ion[i].indx[j]
          ;; Count for elemental
          ;; Loop on the index
          ;; Normal: lls[nn].systems[sys].ion[atnmb].indx
          for k=1,tidx do begin
              if lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 1 OR $
                lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 9 then begin
                  ;; Good measurement
                  ionflg = 1
                  lls_calcclm, lls, nn, sys, i, $
                               lls[nn].systems[sys].ion[i].state[j,k].clm, $
                               lls[nn].systems[sys].ion[i].state[j,k].sigclm
                  lls[nn].systems[sys].elm[i].flgclm = 1
                  lls[nn].systems[sys].elm[i].flginst = $
                    lls_updfd(lls[nn].systems[sys].elm[i].flginst, $
                              lls[nn].systems[sys].ion[i].state[j,k].flginst)
              endif
          endfor

          ;; Saturated
          if ionflg NE 1 then begin
              for k=1,tidx do begin
                  if lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 3 OR $
                    lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 11 then begin
                      ionflg = 2
                      lls[nn].systems[sys].elm[i].flgclm = 2
                      lls_calcclm, lls, nn, sys, i, $
                                   lls[nn].systems[sys].ion[i].state[j,k].clm, $
                                   lls[nn].systems[sys].ion[i].state[j,k].sigclm
                      lls[nn].systems[sys].elm[i].flginst = $
                        lls_updfd(lls[nn].systems[sys].elm[i].flginst, $
                                  lls[nn].systems[sys].ion[i].state[j,k].flginst)
                  endif
              endfor

              ;; Upperlimit
              if ionflg NE 2 then begin
                  for k=1,tidx do begin
                      if lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 5 OR $
                        lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 13 then begin
                          ionflg = 3
                          if lls[nn].systems[sys].elm[i].clm LE 0. then $
                            lls[nn].systems[sys].elm[i].clm = $
                              lls[nn].systems[sys].ion[i].state[j,k].clm $
                          else $
                            lls[nn].systems[sys].elm[i].clm = $
                            lls_sig*lls[nn].systems[sys].ion[i].state[j,k].sigclm < $
                            lls[nn].systems[sys].ion[i].state[j,k].clm
                          ;;
                          lls[nn].systems[sys].elm[i].flgclm = 3
                          lls[nn].systems[sys].elm[i].flginst = $
                            lls_updfd(lls[nn].systems[sys].elm[i].flginst, $
                                      lls[nn].systems[sys].ion[i].state[j,k].flginst)
                      endif
                  endfor
              endif
          endif
          
      endfor
  endfor
	   
  ;;ccccccccccccccccccccccccccccccccccccc
  ;; Update All Ions columns
  ;;cccccccccccccccccccccccccccccccccccc

  ;; Loop on Z
  for i=2,imx do begin
      ;; Loop on ion
      for j=1,6 do begin
          ionflg = 0
          ;; CII*
          if i EQ 6 AND j EQ 6 then continue
          ;; Loop on the index
          tidx = lls[nn].systems[sys].ion[i].indx[j]
          ;; Normal
          for k=1,tidx do begin
              if lls[nn].systems[sys].ion[i].state[j,k].flgclm LT 2 OR $
                lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 9 OR $
                lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 8 then begin
                  ionflg = 1
                  lls_allclm, lls, nn, sys, i, j, $
                              lls[nn].systems[sys].ion[i].state[j,k].clm, $
                    lls[nn].systems[sys].ion[i].state[j,k].sigclm
                  lls[nn].systems[sys].ion[i].state[j].flgclm  = 1
                  lls[nn].systems[sys].ion[i].state[j].flginst = $
                    lls_updfd(lls[nn].systems[sys].ion[i].state[j].flginst, $
                              lls[nn].systems[sys].ion[i].state[j,k].flginst)
              endif
          endfor
          ;; Saturated
          if ionflg NE 1 then begin
              for k=1,tidx do begin
                  if lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 3 OR $
                    lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 2 OR $
                    lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 10 OR $
                    lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 11 then begin
                      ionflg = 2
                      lls[nn].systems[sys].ion[i].state[j].flgclm  = 2
                      lls_allclm, lls, nn, sys, i, j, $
                                  lls[nn].systems[sys].ion[i].state[j,k].clm, $
                        lls[nn].systems[sys].ion[i].state[j,k].sigclm
                      lls[nn].systems[sys].ion[i].state[j].flginst = $
                        lls_updfd(lls[nn].systems[sys].ion[i].state[j].flginst, $
                                  lls[nn].systems[sys].ion[i].state[j,k].flginst)
                  endif
              endfor
              ;; Upperlimit
              if ionflg NE 2 then begin
                  for k=1,tidx do begin
                      if lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 4 OR $
                        lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 5 OR $
                        lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 12 OR $
                        lls[nn].systems[sys].ion[i].state[j,k].flgclm EQ 13 then begin
                          ionflg = 3
                          ;;  Allow for blends for upper limits
                          if lls[nn].systems[sys].ion[i].state[j].clm LE 0. then $
                            lls[nn].systems[sys].ion[i].state[j].clm = $
                            lls[nn].systems[sys].ion[i].state[j,k].clm $
                          else $
                            lls[nn].systems[sys].ion[i].state[j].clm = $
                            lls[nn].systems[sys].ion[i].state[j,k].clm < $
                            lls[nn].systems[sys].ion[i].state[j].clm 
                          ;;
                          lls[nn].systems[sys].ion[i].state[j].flgclm  = 3
                          lls[nn].systems[sys].ion[i].state[j].flginst = $
                            lls_updfd(lls[nn].systems[sys].ion[i].state[j].flginst, $
                                      lls[nn].systems[sys].ion[i].state[j,k].flginst)
                      endif
                  endfor
              endif
          endif
      endfor
  endfor
	   

  ;; Close shop
  free_lun, lun
      
  
  return
end
