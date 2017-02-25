;+ 
; NAME:
; dla_allabd
;  V1.2
;
; PURPOSE:
;    Given a DLA structure, calculate abundances for various
;    transitions using the AODM.  Also fill up elemental abundances
;    using dominant ions.
;
; CALLING SEQUENCE:
;   dla_allabd, dla, nn
;
; INPUTS:
;   dla -- DLA structure
;   nn  -- Index of the structure
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;  /FLG_PLT --  Plot the transitions
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_allabd, dla
;
;
; PROCEDURES CALLED:
;  upd_fd
;
; REVISION HISTORY:
;   01-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function dla_updfd, val, flg

  if val EQ 0 then return, flg
  
  if (val MOD (2*flg)) GT (flg-1) then return, val $
  else return, val + flg
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_allabnd, dla, nn, FLG_PLT=flg_plt

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'dla_allabd, dla, nn, [v1.2]'
    return
  endif 

  if not keyword_set( FINEDAT ) then $
    finedat = getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'

  ;; Fine structure
  readcol, finedat, fin_Z, fin_ion, fin_j, fin_wav, fin_E, fin_A, $
           FORMAT='I,I,F,D,F,F', /silent

  if not keyword_set( DLA_SIG ) then dla_sig = 3.

  ;; Initialize
  dla[nn].ion.state.flgclm = -1
  dla[nn].ion.state[*,0].flgclm = 0

  dumc = ''
  dumf = 0.0
  dumf2 = 0.0
  dumd = 0.d
  dumd2 = 0.d
  dumi = 0

  ;; Open data file
  openr, lun, strtrim(dla[nn].abndfil,2), /get_lun
  readf, lun, dumc
  dlanm = dumc
  ;; Data file(s)
  readf, lun, dumi
  dla[nn].ndfil = dumi

  ;; HIRES [1]
  if dla[nn].ndfil MOD 2 GT 0 then begin
      readf, lun, dumc
      dla[nn].Hfil = dumc
  endif
  ;; ESI  [2]
  if dla[nn].ndfil MOD 4 GT 1 then begin
      readf, lun, dumc
      dla[nn].Efil = dumc
  endif
  ;; UVES [4]
  if dla[nn].ndfil MOD 8 GT 3 then begin
      readf, lun, dumc
      dla[nn].Ufil = dumc
  endif
  ;; X? [8]
  if dla[nn].ndfil MOD 16 GT 7 or dla[nn].ndfil MOD 32 GT 15 then begin
      readf, lun, dumc
      dla[nn].Xfil = dumc
  endif
  ;; MIKE b
  if dla[nn].ndfil MOD 32 GT 15 then begin
      readf, lun, dumc
      ism[nn].MBfil = dumc
  endif
  ;; MIKE r
  if dla[nn].ndfil MOD 64 GT 31 then begin
      readf, lun, dumc
      ism[nn].MRfil = dumc
   endif

  ;; FIRE [16]
  if dla[nn].ndfil MOD 128 GT 63 then begin
      readf, lun, dumc
      dla[nn].Ffil = dumc
   endif

  ;; Redshift
  readf, lun, dumd
;  if dumd NE dla[nn].zabs then begin
  if dumd-dla[nn].zabs gt 1e-10 then begin
      print, 'dla_allabnd: Update zabs fool!', dumd, dla[nn].zabs
      stop
  endif

  ;; Tables
  readf, lun, dumc
  dla[nn].tab_fil = strtrim(dumc,2)

  ;; HI
  readf, lun, dumd, dumd2
  if dumd NE dla[nn].NHI then begin
      print, 'dla_allabnd: Update NHI fool!', dumd, dla[nn].NHI
      stop
  endif
  dla[nn].ion[1].indx[1] = 1
  dla[nn].ion[1].state[1].flgclm = 1
  dla[nn].ion[1].state[1].clm = dumd
  dla[nn].ion[1].state[1].lambda = 1215.6701d
  dla[nn].ion[1].state[1].sigclm = dumd2
  dla[nn].ion[1].state[1,1].flgclm = 1
  dla[nn].ion[1].state[1,1].clm = dumd
  dla[nn].ion[1].state[1,1].lambda = 1215.6701d
  dla[nn].ion[1].state[1,1].sigclm = dumd2
  dla[nn].elm[1].clm = dumd
  dla[nn].elm[1].sigclm = dumd2

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
      dla[nn].elm[Znmb].flginst = fins
      ein = 10^(cin+ein) - 10^cin
      ;; Calc the column (and fill)
      dla_calcclm, dla, nn, ZNmb, 10^cin, ein
      dla[nn].elm[Znmb].flgclm = 1
  endfor

  while not eof(lun) do begin
      ;; Ion flag
      readf, lun, dumi
      ionflg = dumi
      ;; Ion info
      readf, lun, dumd, dumf, dumf2, dumi
      lambda = dumd
      vmin = dumf
      vmax = dumf2
      flg_fil = dumi

      ;; Get atomic info
      getion, lambda, ion, elm, Z=atnmb

      ;; CII*
      if lambda EQ 1335.7077d then ion = 6

      ;; Increment + basic info
;      idx = dla[nn].ion[atnmb].indx[ion]++ + 1
      idx = dla[nn].ion[atnmb].indx[ion] + 1
      dla[nn].ion[atnmb].indx[ion] = dla[nn].ion[atnmb].indx[ion] + 1
      dla[nn].ion[atnmb].state[ion,idx].lambda = lambda
      dla[nn].ion[atnmb].state[ion,idx].flginst = flg_fil
      

      ;; Deal with by hand columns
      if ionflg GE 8 then begin
          if ionflg GE 16 then continue  ;; EW by hand
          dla[nn].ion[atnmb].state[ion,idx].clm = 10.^vmin
          dla[nn].ion[atnmb].state[ion,idx].sigclm = $
            (10^(vmin+vmax) - 10^(vmin-vmax))/2	      
          dla[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
          continue
      endif

      ;; Read in data
      if not keyword_set( SV_FIL ) then sv_fil = -1
      ;mdflg = flg_fil mod 64
      ;if flg_fil GT 128 then iflg=2 else iflg=0
      ;;; MR had to modify the above to make
      ;;; it work with the higher numbered value
       mdflg = flg_fil mod 128
      if flg_fil GT 128 then begin
         if flg_fil LT 256 then iflg=2 else iflg = 13
      endif else iflg=0
     if mdflg NE SV_FIL then begin
          case mdflg of
              1: dat = x_readspec(dla[nn].Hfil, /struct, /autofsig)
              2: dat = x_readspec(dla[nn].Efil, /struct, /autofsig)
              4: dat = x_readspec(dla[nn].Ufil, /struct, /autofsig, INFLG=iflg)
              8: dat = x_readspec(dla[nn].Xfil, /struct, /autofsig, INFLG=iflg)
              64: dat = x_readspec(dla[nn].Ffil, /struct, /autofsig) ;; FIRE
              else:stop
          endcase
      endif else sv_fil = mdflg

      ;if flg_fil EQ 136 then stop
      ;; Velocity array
      dumv = (x_allvelo(dat.wv, dla[nn].zabs, lambda, $
                            [vmin, vmax], ALL_PMNX=pmnx))[*]
      dat.velo[pmnx[0]:pmnx[1]] = dumv[0:pmnx[2,0]]

      ;; AODM
      x_aodm, dat.wv, dat.fx[pmnx[0]:pmnx[1]], $
        dat.sig[pmnx[0]:pmnx[1]], lambda, clm, sigc, $
        VELO=dumv[0:pmnx[2,0]], FLG_SAT=flg_sat
      dla[nn].ion[atnmb].state[ion,idx].clm = clm
      dla[nn].ion[atnmb].state[ion,idx].sigclm = sigc

      ;;;;;;;
      ;; Flag

      ;; CII*
      if lambda EQ 1335.7077d then begin
          ;; Saturation
          if flg_sat GT 0 then begin
              dla[nn].ion[atnmb].state[ion,idx].flgclm = 2
              continue
          endif
          if clm LT dla_sig*sigc then begin
              dla[nn].ion[atnmb].state[ion,idx].flgclm = 4
              dla[nn].ion[atnmb].state[ion,idx].clm = dla_sig*sigc
              continue
          endif
          if ionflg GE 4  then begin
              dla[nn].ion[atnmb].state[ion,idx].clm = clm > dla_sig*sigc
              dla[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
              continue
          endif
      endif

      ;; Standard
      if flg_sat GT 0 OR ionflg EQ 2 OR ionflg EQ 3  then begin
          ;; Lower limits
          if (ionflg MOD 8) gt 1 then begin ; Keep flg given by hand
              dla[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
              print, 'dla_allabnd: Taking value in the file'
          endif else dla[nn].ion[atnmb].state[ion,idx].flgclm = ionflg + 2
      endif else begin

          ;; Upper limit?  ion_flg = 4,5
          if clm GT dla_sig*sigc then begin
              dla[nn].ion[atnmb].state[ion,idx].flgclm = ionflg
              ;; Blends
              if ionflg EQ 4 OR ionflg EQ 5 then $
                dla[nn].ion[atnmb].state[ion,idx].clm = dla_sig*sigc > clm 
          endif else begin
              ;; Upper limit!!
              if ionflg LT 4 then begin
                  ;; Adopt dla_sig (3) sigma limit
                  dla[nn].ion[atnmb].state[ion,idx].flgclm = ionflg + 4
                  dla[nn].ion[atnmb].state[ion,idx].clm = dla_sig*sigc 
              endif else begin
                  dla[nn].ion[atnmb].state[ion,idx].flgclm = ionflg 
                  dla[nn].ion[atnmb].state[ion,idx].clm = $
                    dla_sig*sigc > clm
              endelse
          endelse
      endelse
  endwhile

  ;;ccccccccccccccccccccccccccccccccccccc
  ;;  Update elemental column if appropriate
  ;;ccccccccccccccccccccccccccccccccccccc

  ;; Loop on Z
  for i=2,100 do begin
      ;; Loop on ion
      for j=1,6 do begin
          ionflg = 0
          ;; Deal with by hand!
          if dla[nn].elm[i].flgclm EQ 1 then ionflg=1
          ;; CII*
          if i EQ 6 AND j EQ 6 then continue
          tidx = dla[nn].ion[i].indx[j]
          ;; Count for elemental
          ;; Loop on the index
          ;; Normal: dla[nn].ion[atnmb].indx
          for k=1,tidx do begin
              if dla[nn].ion[i].state[j,k].flgclm EQ 1 OR $
                dla[nn].ion[i].state[j,k].flgclm EQ 9 then begin
                  ;; Good measurement
                  ionflg = 1
                  dla_calcclm, dla, nn, i, dla[nn].ion[i].state[j,k].clm, $
                               dla[nn].ion[i].state[j,k].sigclm
                  dla[nn].elm[i].flgclm = 1
                  dla[nn].elm[i].flginst = $
                    dla_updfd(dla[nn].elm[i].flginst, $
                              dla[nn].ion[i].state[j,k].flginst)
              endif
          endfor

          ;; Saturated
          if ionflg NE 1 then begin
              for k=1,tidx do begin
                  if dla[nn].ion[i].state[j,k].flgclm EQ 3 OR $
                    dla[nn].ion[i].state[j,k].flgclm EQ 11 then begin
                      ionflg = 2
                      dla[nn].elm[i].flgclm = 2
                      dla_calcclm, dla, nn, i, dla[nn].ion[i].state[j,k].clm, $
                        dla[nn].ion[i].state[j,k].sigclm
                      dla[nn].elm[i].flginst = dla_updfd(dla[nn].elm[i].flginst, $
                                                         dla[nn].ion[i].state[j,k].flginst)
                  endif
              endfor

              ;; Upperlimit
              if ionflg NE 2 then begin
                  for k=1,tidx do begin
                      if dla[nn].ion[i].state[j,k].flgclm EQ 5 OR $
                        dla[nn].ion[i].state[j,k].flgclm EQ 13 then begin
                          ionflg = 3
                          if dla[nn].elm[i].clm LE 0. then $
                            dla[nn].elm[i].clm = dla[nn].ion[i].state[j,k].clm $
                          else $
                            dla[nn].elm[i].clm = $
                            dla[nn].ion[i].state[j,k].clm < dla[nn].elm[i].clm
                          ;;
                          dla[nn].elm[i].flgclm = 3
                          dla[nn].elm[i].flginst = $
                            dla_updfd(dla[nn].elm[i].flginst, $
                                      dla[nn].ion[i].state[j,k].flginst)
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
  for i=2,100 do begin
      ;; Loop on ion
      for j=1,6 do begin
          ionflg = 0
          ;; CII*
          if i EQ 6 AND j EQ 6 then continue
          ;; Loop on the index
          tidx = dla[nn].ion[i].indx[j]
          ;; Normal
          for k=1,tidx do begin
              if min(abs(dla[nn].ion[i].state[j,k].lambda-fin_wav)) LT 0.005 $
                then continue
              if dla[nn].ion[i].state[j,k].flgclm LT 2 OR $
                dla[nn].ion[i].state[j,k].flgclm EQ 9 OR $
                dla[nn].ion[i].state[j,k].flgclm EQ 8 then begin
                  ionflg = 1
                  dla_allclm, dla, nn, i, j, dla[nn].ion[i].state[j,k].clm, $
                    dla[nn].ion[i].state[j,k].sigclm
                  dla[nn].ion[i].state[j].flgclm  = 1
                  dla[nn].ion[i].state[j].flginst = $
                    dla_updfd(dla[nn].ion[i].state[j].flginst, $
                              dla[nn].ion[i].state[j,k].flginst)
              endif
          endfor
          ;; Saturated
          if ionflg NE 1 then begin
              for k=1,tidx do begin
                  if min(abs(dla[nn].ion[i].state[j,k].lambda-fin_wav)) LT 0.005 $
                    then continue
                  if dla[nn].ion[i].state[j,k].flgclm EQ 3 OR $
                    dla[nn].ion[i].state[j,k].flgclm EQ 2 OR $
                    dla[nn].ion[i].state[j,k].flgclm EQ 10 OR $
                    dla[nn].ion[i].state[j,k].flgclm EQ 11 then begin
                      ionflg = 2
                      dla[nn].ion[i].state[j].flgclm  = 2
                      dla_allclm, dla, nn, i, j, dla[nn].ion[i].state[j,k].clm, $
                        dla[nn].ion[i].state[j,k].sigclm
                      dla[nn].ion[i].state[j].flginst = $
                        dla_updfd(dla[nn].ion[i].state[j].flginst, $
                                  dla[nn].ion[i].state[j,k].flginst)
                  endif
              endfor
              ;; Upperlimit
              if ionflg NE 2 then begin
                  for k=1,tidx do begin
;                      if i EQ 6 and j EQ 2 then stop
                  if min(abs(dla[nn].ion[i].state[j,k].lambda-fin_wav)) LT 0.005 $
                        then continue
                      if dla[nn].ion[i].state[j,k].flgclm EQ 4 OR $
                        dla[nn].ion[i].state[j,k].flgclm EQ 5 OR $
                        dla[nn].ion[i].state[j,k].flgclm EQ 12 OR $
                        dla[nn].ion[i].state[j,k].flgclm EQ 13 then begin
                          ionflg = 3
                          ;;  Allow for blends for upper limits
                          if dla[nn].ion[i].state[j].clm LE 0. then $
                            dla[nn].ion[i].state[j].clm = $
                            dla[nn].ion[i].state[j,k].clm $
                          else $
                            dla[nn].ion[i].state[j].clm = $
                            dla[nn].ion[i].state[j,k].clm < $
                            dla[nn].ion[i].state[j].clm 
                          ;;
                          dla[nn].ion[i].state[j].flgclm  = 3
                          dla[nn].ion[i].state[j].flginst = $
                            dla_updfd(dla[nn].ion[i].state[j].flginst, $
                                      dla[nn].ion[i].state[j,k].flginst)
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
