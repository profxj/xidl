;+ 
; NAME:
; fuse_calcxh
;  V1.1
;
; PURPOSE:
;    Calculate X/H values for FUSE observations.  It grabs the info
;    from the ion_fil and an input XH_fil and then outputs the XH
;    values into a XH data file.  The Input file must contain
;    ionization information (Cloudy file, U values) to do the final
;    calcalculations of X/H.  You need to see an example file to get
;    this set correctly.
;
; CALLING SEQUENCE:
;   fuse_calcxh, xhfil
;
; INPUTS:
;  xhfil -- Input file to run the XH calculations and output
;
; RETURNS:
;
; OUTPUTS:
;  outfil -- Output file of X/H values
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calcxh, 'Input/z0495_XH.inp'
;
; PROCEDURES CALLED:
;  fuse_calcxh_parse
;
; REVISION HISTORY:
;   21-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_calcxh_parse, xhfil, ion_fil, cldy, pht_flg, $
                       Uflg, Ubest, NHI, mtl, outfil
                       

  ion_fil= ''
  cfil= ''
  outfil= ''
  pht_flg = 0
  uflg = 0
  NHI = 0.d
  mtl = 0.d
  Ubest = fltarr(3)
  ;; Open
  close, 1
  openr, 1, xhfil
  readf, 1, ion_fil, FORMAT='(a)'
  readf, 1, pht_flg, FORMAT='(i)'
  case pht_flg of
      0: readf, 1, outfil, FORMAT='(a)'
      1: begin
          ;; Open cloudy
          readf, 1, cfil, FORMAT='(a)'
          prs_cldygrid, cldy, cfil
          ;; Grab Cloudy info
          readf, 1, Uflg, FORMAT='(i)'
          readf, 1, Ubest, FORMAT='(3f)'
          readf, 1, NHI, mtl, FORMAT='(2f)'
          readf, 1, outfil, FORMAT='(a)'
      end
      else: stop
  endcase
  close, 1

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fuse_calcxh, xhfil

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'fuse_calcxh, xhfil (v1.1)' 
    return
  endif 

  if not keyword_set( nH ) then nH = -1.d


  ;; Parse xhfil
  close, /all
  fuse_calcxh_parse, xhfil, ion_fil, cldy, pht_flg, Uflg, $
    Ubest, NHI, mtl, outfil
  

  ;; Outfil
  openw, 2, outfil
          
  ;; Read ion_fil
  readcol, ion_fil, atmnum, ionnum, Nval, Nsig, Nflg, FORMAT='I,I,F,F,I', /silent
  nion = n_elements(ionnum)

  ;; Grep out 
  gd = where(cldy.NHI EQ NHI AND cldy.FeH EQ mtl AND cldy.nH EQ nH, ngd)
  if ngd EQ 0 then stop
  scldy = cldy[gd]
  
  ;; HI
  printf, 2, atmnum[0], ionnum[0], Nval[0], Nsig[0], '1', $
    FORMAT='(i2,1x,i2,1x,f6.3,1x,f5.3,1x,a1)'
  ;; HII
  idum = 2
  if Nflg[0] NE 1 then stop
  case uflg of 
      1: begin
          ;; best
          pmod = scldy.X[1,1]
          IC = interpol(pmod, scldy.U, Ubest)
          XH = Nval[0] - IC
          flg = 1
          NHII = XH[0]
          sigHII = total( abs(XH-XH[0]) ) /2.
      end
      2: begin ;; Lower limit on U
          uindx = where(scldy.U GE Ubest[0], nui)
          ucldy = scldy[uindx]
          if nui EQ 0 then stop
          IC = ucldy.X[1,1]
          mxIC = max(IC)
          NHII = Nval[0] - mxIC
          flg = 2
          sigHII = 0.
      end
      else: stop
  endcase
  printf, 2, atmnum[0], idum, NHII, sigHII, flg, $
    FORMAT='(i2,1x,i2,1x,f6.3,1x,f5.3,1x,i1)'

  ;; Loop on metals
  for qq=1L, nion-1 do begin

      ;; Abund
      getabnd, elm, atmnum[qq], abnd, flag=1

      case uflg of 
          1: begin
              ;; best
              pmod = scldy.X[atmnum[qq],ionnum[qq]]-scldy.X[1,1]
              IC = interpol(pmod, scldy.U, Ubest)
              XH = Nval[qq] - Nval[0] + 12. - abnd - IC
              sig = 0.
              if Nflg[qq] EQ 2 or Nflg[qq] EQ 3 then begin
                  flg = 2
                  XH[0] = min(XH)
              endif
              if Nflg[qq] EQ 4 or Nflg[qq] EQ 5 then begin
                  flg = 3
                  XH[0] = max(XH)
              endif
              if Nflg[qq] EQ 0 or Nflg[qq] EQ 1 then begin
                  flg = 1
                  sig = total( abs(XH-XH[0]) ) /2.
              endif
          end
          2: begin ;; Lower limit on U
              flg = -1
              if Nflg[qq] EQ 1 OR Nflg[qq] EQ 3 then begin
                  uindx = where(scldy.U GE Ubest[0], nui)
                  ucldy = scldy[uindx]
                  if nui EQ 0 then stop
                  IC = ucldy.X[atmnum[qq],ionnum[qq]]-ucldy.X[1,1]
                  mxIC = max(IC)
                  XH = Nval[qq] - Nval[0] + 12. - abnd - mxIC
                  flg = 2
                  sig = 0.
              endif
          end
          else: stop
      endcase
              
              

      ;; Print
      if flg NE -1 then $
        printf, 2, atmnum[qq], ionnum[qq], XH[0], sig, flg, $
        FORMAT='(i2,1x,i2,1x,f6.3,1x,f5.3,1x,i1)'
  endfor

  close, /all

  return
end
