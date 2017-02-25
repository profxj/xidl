;+ 
; NAME:
; ism_allabd
;  V1.1
;
; PURPOSE:
;    Calculates the abundances for a series of absorption lines given
;    an abundance structure.  Similar to dla_allabd
;
; CALLING SEQUENCE:
;   ism_allabd, ism, nn
;
; INPUTS:
;   ism -- DLA structure
;   nn  -- Index of the structure
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ism_allabd, ism
;
;
; PROCEDURES CALLED:
;  upd_fd
;
; REVISION HISTORY:
;   08-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function ism_updfd, val, flg

  if val EQ 0 then return, flg
  
  if (val MOD 2*flg) GT (flg-1) then return, val $
  else return, val + flg
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ism_allabnd, ism, nn, FLG_PLT=flg_plt, TYP_FLG=typ_flg

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'ism_allabd, ism, nn, [v1.1]'
    return
  endif 

  if not keyword_set( ISM_SIG ) then ism_sig = 3.
;  if not keyword_set( CONTI_ERR ) then conti_err = 0.
;  resolve, 'ew_to_colm'

  ;; Initialize
  ism[nn].ion.state.flgclm = -1
  ism[nn].ion.state[*,0].flgclm = 0

  dumc = ''
  dumf = 0.0
  dumf2 = 0.0
  dumf3 = 0.0
  dumd = 0.d
  dumd2 = 0.d
  dumi = 0

  ;; Open data file
  openr, lun, strtrim(ism[nn].abndfil,2), /get_lun
  readf, lun, dumc
  ismnm = dumc
  ;; Data file(s)
  readf, lun, dumi
  ism[nn].ndfil = dumi

  ;; HIRES
  if ism[nn].ndfil MOD 2 GT 0 then begin
      readf, lun, dumc
      ism[nn].Hfil = dumc
  endif
  ;; ESI 
  if ism[nn].ndfil MOD 4 GT 1 then begin
      readf, lun, dumc
      ism[nn].Efil = dumc
  endif
  ;; UVES
  if ism[nn].ndfil MOD 8 GT 3 then begin
      readf, lun, dumc
      ism[nn].Ufil = dumc
  endif
  ;; X?
  if ism[nn].ndfil MOD 16 GT 7 then begin
      readf, lun, dumc
      ism[nn].Xfil = dumc
  endif
  ;; INFLG=2
  if ism[nn].ndfil MOD 64 GT 31 then begin
      readf, lun, dumc
      ism[nn].Xfil = dumc
  endif
  ;; Redshift
  readf, lun, dumd
  if dumd NE ism[nn].zabs then begin
      print, 'ism_allabnd: Update zabs fool!'
      stop
  endif

  ;; Tables
  readf, lun, dumc
  ism[nn].tab_fil = strtrim(dumc,2)

  ;; HI
  readf, lun, dumd, dumd2
  if dumd NE ism[nn].NHI then begin
      print, 'ism_allabnd: Update NHI fool!'
      stop
  endif
  ism[nn].ion[1].indx[1] = 1
  ism[nn].ion[1].state[1].flgclm = 1
  ism[nn].ion[1].state[1].clm = dumd
  ism[nn].ion[1].state[1].lambda = 1215.6701d
  ism[nn].ion[1].state[1].sigclm = dumd2
  ism[nn].ion[1].state[1,1].flgclm = 1
  ism[nn].ion[1].state[1,1].clm = dumd
  ism[nn].ion[1].state[1,1].lambda = 1215.6701d
  ism[nn].ion[1].state[1,1].sigclm = dumd2
  ism[nn].elm[1].clm = dumd
  ism[nn].elm[1].sigclm = dumd2

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
      ism[nn].elm[Znmb].flginst = fins
      ein = 10^(cin+ein) - 10^cin
      ;; Calc the column (and fill)
      dla_calcclm, ism, nn, ZNmb, 10^cin, ein
      ism[nn].elm[Znmb].flgclm = 1
  endfor

  while not eof(lun) do begin
      ;; Ion flag
      readf, lun, dumi
      ionflg = dumi
      ;; Ion info
      readf, lun, dumd, dumf, dumf2, dumi, dumf3
      lambda = dumd
      vmin = dumf
      vmax = dumf2
      flg_fil = dumi
      contierr = dumf3

      ;; Get atomic info
      getion, lambda, ion, elm, Z=atnmb

      ;; CII*
      if lambda EQ 1335.7077d then ion = 6

      ;; Increment + basic info
      idx = ism[nn].ion[atnmb].indx[ion]++ + 1
      ism[nn].ion[atnmb].state[ion,idx].lambda = lambda
      ism[nn].ion[atnmb].state[ion,idx].flginst = flg_fil
      

      ;; Deal with by hand columns
      if ionflg GE 8 then begin
          ism[nn].ion[atnmb].state[ion,idx].clm = 10^vmin
          ism[nn].ion[atnmb].state[ion,idx].sigclm = $
            (10^(vmin+vmax) - 10^(vmin-vmax))/2	      
          ism[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
          continue
      endif

      ;; Read in data
      if not keyword_set( SV_FIL ) then sv_fil = -1
      mdflg = flg_fil mod 64
      if mdflg NE SV_FIL then begin
          case mdflg of
              1: dat = x_readspec(ism[nn].Hfil, /struct, /autofsig)
              2: dat = x_readspec(ism[nn].Efil, /struct, /autofsig)
              4: dat = x_readspec(ism[nn].Ufil, /struct, /autofsig)
              8: dat = x_readspec(ism[nn].Xfil, /struct, /autofsig)
              32: dat = x_readspec(ism[nn].Xfil, /struct, /autofsig,inflg=2)
          endcase
      endif else sv_fil = mdflg

      ;; Velocity array
      dumv = (x_allvelo(dat.wv, ism[nn].zabs, lambda, $
                            [vmin, vmax], ALL_PMNX=pmnx))[*]
      dat.velo[pmnx[0]:pmnx[1]] = dumv[0:pmnx[2,0]]

      ;; AODM
      x_aodm, dat.wv, dat.fx[pmnx[0]:pmnx[1]], $
        dat.sig[pmnx[0]:pmnx[1]], lambda, clm, sigc, $
        VELO=dumv[0:pmnx[2,0]], FLG_SAT=flg_sat
      
      ;; Continuum error
      if CONTIERR GT 0. then begin
          dw = dat.wv[pmnx[1]] - dat.wv[pmnx[0]]
          sig_conti = contierr * dw * 1000. ;; mA
          ;; Convert to column
          sig_ci = ew_to_colm([lambda], [sig_conti], /silent)
          ;; Add in quadrature
          sigc = sqrt(sigc^2 + sig_ci[0]^2)
      endif
      ism[nn].ion[atnmb].state[ion,idx].clm = clm
      ism[nn].ion[atnmb].state[ion,idx].sigclm = sigc

      ;;;;;;;
      ;; Flag

      ;; CII*
      if lambda EQ 1335.7077d then begin
          ;; Saturation
          if flg_sat GT 0 then begin
              ism[nn].ion[atnmb].state[ion,idx].flgclm = 2
              continue
          endif
          if clm LT ism_sig*sigc then begin
              ism[nn].ion[atnmb].state[ion,idx].flgclm = 4
              ism[nn].ion[atnmb].state[ion,idx].clm = ism_sig*sigc
              continue
          endif
          if ionflg GE 4  then begin
              ism[nn].ion[atnmb].state[ion,idx].clm = clm > ism_sig*sigc
              ism[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
              continue
          endif
      endif

      ;; Standard
      if flg_sat GT 0 OR ionflg EQ 2 OR ionflg EQ 3  then begin
          ;; Lower limits
          if (ionflg MOD 8) gt 1 then begin ; Keep flg given by hand
              ism[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
              print, 'ism_allabnd: Taking value in the file'
          endif else ism[nn].ion[atnmb].state[ion,idx].flgclm = ionflg + 2
      endif else begin

          ;; Upper limit?  ion_flg = 4,5
          if clm GT ism_sig*sigc then begin
              ism[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
              ;; Blends
              if ionflg EQ 4 OR ionflg EQ 5 then $
                ism[nn].ion[atnmb].state[ion,idx].clm = ism_sig*sigc > clm 
          endif else begin
              ;; Upper limit!!
              if ionflg LT 4 then begin
                  ism[nn].ion[atnmb].state[ion,idx].flgclm = ionflg + 4
                  ism[nn].ion[atnmb].state[ion,idx].clm = ism_sig*sigc 
              endif else begin
                  ism[nn].ion[atnmb].state[ion,idx].flgclm = ionflg 
                  ism[nn].ion[atnmb].state[ion,idx].clm = $
                    ism_sig*sigc > clm
              endelse
          endelse
      endelse
  endwhile

  ;;ccccccccccccccccccccccccccccccccccccc
  ;;  Update elemental column if appropriate
  ;;ccccccccccccccccccccccccccccccccccccc

  ;; Loop on Z
  imx = 100 < ((size(ism.elm))[1]-1)
  for i=2,imx do begin
      ;; Loop on ion
      for j=1,6 do begin
          ionflg = 0
          ;; Deal with by hand!
          if ism[nn].elm[i].flgclm EQ 1 then ionflg=1
          ;; CII*
          if i EQ 6 AND j EQ 6 then continue
          tidx = ism[nn].ion[i].indx[j]
          ;; Count for elemental
          ;; Loop on the index
          ;; Normalism[nn].ion[atnmb].indx
          for k=1,tidx do begin
              if ism[nn].ion[i].state[j,k].flgclm EQ 1 OR $
                ism[nn].ion[i].state[j,k].flgclm EQ 9 then begin
                  ;; Good measurement
                  ionflg = 1
                  dla_calcclm, ism, nn, i, ism[nn].ion[i].state[j,k].clm, $
                    ism[nn].ion[i].state[j,k].sigclm
                  ism[nn].elm[i].flgclm = 1
                  ism[nn].elm[i].flginst = ism_updfd(ism[nn].elm[i].flginst, $
                                                     ism[nn].ion[i].state[j,k].flginst)
              endif
          endfor

          ;; Saturated
          if ionflg NE 1 then begin
              for k=1,tidx do begin
                  if ism[nn].ion[i].state[j,k].flgclm EQ 3 OR $
                    ism[nn].ion[i].state[j,k].flgclm EQ 11 then begin
                      ionflg = 2
                      ism[nn].elm[i].flgclm = 2
                      dla_calcclm, ism, nn, i, ism[nn].ion[i].state[j,k].clm, $
                        ism[nn].ion[i].state[j,k].sigclm
                      ism[nn].elm[i].flginst = ism_updfd(ism[nn].elm[i].flginst, $
                                                         ism[nn].ion[i].state[j,k].flginst)
                  endif
              endfor

              ;; Upperlimit
              if ionflg NE 2 then begin
                  for k=1,tidx do begin
                      if ism[nn].ion[i].state[j,k].flgclm EQ 5 OR $
                        ism[nn].ion[i].state[j,k].flgclm EQ 13 then begin
                          ionflg = 3
                          if ism[nn].elm[i].clm LE 0. then $
                            ism[nn].elm[i].clm = ism[nn].ion[i].state[j,k].clm $
                          else $
                            ism[nn].elm[i].clm = $
                            ism_sig*ism[nn].ion[i].state[j,k].sigclm < $
                            ism[nn].ion[i].state[j,k].clm
                          ;;
                          ism[nn].elm[i].flgclm = 3
                          ism[nn].elm[i].flginst = $
                            ism_updfd(ism[nn].elm[i].flginst, $
                                      ism[nn].ion[i].state[j,k].flginst)
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
          tidx = ism[nn].ion[i].indx[j]
          ;; Normal
          for k=1,tidx do begin
              if ism[nn].ion[i].state[j,k].flgclm LT 2 OR $
                ism[nn].ion[i].state[j,k].flgclm EQ 9 OR $
                ism[nn].ion[i].state[j,k].flgclm EQ 8 then begin
                  ionflg = 1
                  dla_allclm, ism, nn, i, j, ism[nn].ion[i].state[j,k].clm, $
                    ism[nn].ion[i].state[j,k].sigclm
                  ism[nn].ion[i].state[j].flgclm  = 1
                  ism[nn].ion[i].state[j].flginst = $
                    ism_updfd(ism[nn].ion[i].state[j].flginst, $
                              ism[nn].ion[i].state[j,k].flginst)
              endif
          endfor
          ;; Saturated
          if ionflg NE 1 then begin
              for k=1,tidx do begin
                  if ism[nn].ion[i].state[j,k].flgclm EQ 3 OR $
                    ism[nn].ion[i].state[j,k].flgclm EQ 2 OR $
                    ism[nn].ion[i].state[j,k].flgclm EQ 10 OR $
                    ism[nn].ion[i].state[j,k].flgclm EQ 11 then begin
                      ionflg = 2
                      ism[nn].ion[i].state[j].flgclm  = 2
                      dla_allclm, ism, nn, i, j, ism[nn].ion[i].state[j,k].clm, $
                        ism[nn].ion[i].state[j,k].sigclm
                      ism[nn].ion[i].state[j].flginst = $
                        ism_updfd(ism[nn].ion[i].state[j].flginst, $
                                  ism[nn].ion[i].state[j,k].flginst)
                  endif
              endfor
              ;; Upperlimit
              if ionflg NE 2 then begin
                  for k=1,tidx do begin
                      if ism[nn].ion[i].state[j,k].flgclm EQ 4 OR $
                        ism[nn].ion[i].state[j,k].flgclm EQ 12 OR $
                        ism[nn].ion[i].state[j,k].flgclm EQ 13 then begin
                          ionflg = 3
                          ;;  Allow for blends for upper limits
                          if ism[nn].ion[i].state[j].clm LE 0. then $
                            ism[nn].ion[i].state[j].clm = $
                            ism[nn].ion[i].state[j,k].clm $
                          else $
                            ism[nn].ion[i].state[j].clm = $
                            ism[nn].ion[i].state[j,k].clm < $
                            ism[nn].ion[i].state[j].clm 
                          ;;
                          ism[nn].ion[i].state[j].flgclm  = 3
                          ism[nn].ion[i].state[j].flginst = $
                            ism_updfd(ism[nn].ion[i].state[j].flginst, $
                                      ism[nn].ion[i].state[j,k].flginst)
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
