;+ 
; NAME:
; ism_calcew
;  V1.1
;
; PURPOSE:
;    Calculates the EW values for a series of transitions given the
;    an input .clm file.  Output is put into a COG structure
;
; CALLING SEQUENCE:
;   ism_calcew, clmfil, all_cog, LLIST=, /NOCONTI, ISM=
;
; INPUTS:
;   clmfil -- ASCII .clm file
;
; RETURNS:
;
; OUTPUTS:
;  all_cog  -- COG structure containing the EW values
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-Apr-2006 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ism_calcew, clmfil, all_cog, LLIST=llist, NOCONTI=noconti, ISM=ism, $
                ROOT=root

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'ism_calcew, clmfil, cog, LLIST= [v1.1]'
    return
  endif 

  if not keyword_set( ROOT ) then root = ''
  if not keyword_set( FINEDAT ) then $
    finedat = getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'

  ;; Fine structure
  readcol, finedat, fine_Z, fine_ion, fine_j, fine_wav, fine_E, fine_A, $
           FORMAT='I,I,F,D,F,F', /silent

  if not keyword_set( ISM_SIG ) then ism_sig = 3.
  ism = {ismstruct}
  nn = 0

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
  a = findfile(strtrim(root+clmfil,2),count=na)
  if na EQ 0 then begin
      if keyword_set(all_cog) then delvarx, all_cog
      print, 'ism_calcew: No clm file!', strtrim(root+clmfil,2)
      return
  endif
  openr, lun, strtrim(root+clmfil,2), /get_lun
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
  ;; MIKE b
  if ism[nn].ndfil MOD 32 GT 15 then begin
      readf, lun, dumc
      ism[nn].MBfil = dumc
  endif
  ;; MIKE r
  if ism[nn].ndfil MOD 64 GT 31 then begin
      readf, lun, dumc
      ism[nn].MRfil = dumc
  endif
  ;; Unusual
  if ism[nn].ndfil MOD 128 GT 63 then begin
      readf, lun, dumc
      ism[nn].Xfil = dumc
  endif

  ;; Redshift
  readf, lun, dumd
  ism[nn].zabs = dumd

  ;; Tables
  readf, lun, dumc
  ism[nn].tab_fil = strtrim(dumc,2)

  ;; HI
  readf, lun, dumd, dumd2
  ism[nn].NHI  = dumd

  ;; Stuff
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

  svstate = lonarr(100,3)  ; Z, state, J level [0,1,2,etc.]

  ncog = 0
  all_cog = replicate({cogstrct},1)
  while not eof(lun) do begin
      ;; Ion flag
      readf, lun, dumi
      ionflg = dumi
      ;; Ion info
      if not keyword_set(NOCONTI) then begin
          readf, lun, dumd, dumf, dumf2, dumi, dumf3
          contierr = dumf3
      endif else readf, lun, dumd, dumf, dumf2, dumi
      if ionflg GE 8 then begin
          if ionflg LT 16 then continue  ;; Colm measurement by hand
          flg_hand = 1
      endif else flg_hand = 0

      lambda = dumd
      vmin = dumf
      vmax = dumf2
      flg_fil = dumi

      ;; Get atomic info
      getion, lambda, ion, elm, Z=atnmb

      ;; Skip?
      if keyword_set(LLIST) then begin
          mn = min(abs(lambda - LLIST),imn)
          if mn GT 0.002 then continue
      endif

      ;; Fine-structure?
      mt = (where(abs(lambda-fine_wav) LT 0.01, nmt))[0]
      if nmt EQ 0 then Jval=0 else begin
          fE = fine_E[where(fine_Z EQ atnmb and fine_ion EQ ion)]
          uE = [0., fE[ uniq( fE, sort(fE) ) ]]
          mtE = where(abs(fine_E[mt] - uE) LT 0.01)
          Jval = mtE[0]
      endelse
      Jval = Jval[0]

      ;; New level?
      mtlvl = where(svstate[*,0] EQ atnmb AND $
                    svstate[*,1] EQ ion AND $
                    svstate[*,2] EQ Jval, nmtlv)
      if nmtlv EQ 0 then begin
          ncog = ncog+1
          all_cog = [all_cog, {cogstrct}]
          svstate[ncog-1,0] = atnmb
          svstate[ncog-1,1] = ion
          svstate[ncog-1,2] = Jval
          icog = ncog
          all_cog[icog].ionnm = elm+'_'+strtrim(ion,2)+'_'+strtrim(Jval,2)
          all_cog[icog].ion = ion
          all_cog[icog].Z = atnmb
      endif else icog = mtlvl[0]+1  ; Funny offset
      all_cog[icog].nlin = all_cog[icog].nlin + 1
          
      ;; Increment + basic info
      idx = ism[nn].ion[atnmb].indx[ion]++ + 1
      ism[nn].ion[atnmb].state[ion,idx].lambda = lambda
      ism[nn].ion[atnmb].state[ion,idx].flginst = flg_fil
      
      ;; Deal with by hand columns (NG for EW)
      if ionflg GE 16 then begin
          all_cog[icog].ew[all_cog[icog].nlin-1] = vmin
          all_cog[icog].sigew[all_cog[icog].nlin-1] = vmax
          all_cog[icog].wrest[all_cog[icog].nlin-1] = lambda
          continue
      endif


      ;; Read in data
      if not keyword_set( SV_FIL ) then sv_fil = -1
      mdflg = flg_fil mod 64
      if flg_fil GT 128 then begin
         if flg_fil LT 256 then iflg=2 else iflg = 13
      endif else iflg=0
      if mdflg NE SV_FIL then begin
          case mdflg of
              1: dat = x_readspec(ism[nn].Hfil, /struct, /autofsig)
              2: dat = x_readspec(ism[nn].Efil, /struct, /autofsig)
              4: dat = x_readspec(ism[nn].Ufil, /struct, /autofsig,INFLG=iflg)
              8: dat = x_readspec(ism[nn].Xfil, /struct, /autofsig, inflg=iflg)
              13: dat = x_readspec(ism[nn].Ufil, /struct, /autofsig, INFLG=13) ;; Murphy UVES
              16: dat = x_readspec(ism[nn].MBfil, /struct, /autofsig)
              32: dat = x_readspec(ism[nn].MRfil, /struct, /autofsig)
              64: dat = x_readspec(ism[nn].Xfil, /struct, /autofsig,inflg=2)
          endcase
      endif else sv_fil = mdflg

      ;; Measure EW
      x_pixminmax, dat.wv, lambda, ism[nn].zabs, vmin, vmax, $
                   PIXMIN=pmin, PIXMAX=pmax
      ew = x_calcew(dat.wv, dat.fx, [pmin,pmax], $
                    dat.sig, sigew, /fpix)

      ;if flg_fil EQ 136 then stop
      ;; Save
      all_cog[icog].dv[all_cog[icog].nlin-1,*] = [vmin,vmax]
      all_cog[icog].ew[all_cog[icog].nlin-1] = ew / (1. + ism[nn].zabs)
      all_cog[icog].sigew[all_cog[icog].nlin-1] = sigew / (1. + ism[nn].zabs)
      all_cog[icog].wrest[all_cog[icog].nlin-1] = lambda
      getfnam, lambda, fval
      all_cog[icog].f[all_cog[icog].nlin-1] = fval
      all_cog[icog].inst[all_cog[icog].nlin-1] = flg_fil

  endwhile

  ;; Close shop
  free_lun, lun

  ;; Chop off the padding
  if n_elements(all_cog) EQ 1 then begin
      delvarx, all_cog
      return
  endif
  all_cog = all_cog[1:*]

  ;; Sort
  bs = all_cog.z*100 + all_cog.ion
  srt = sort(bs)
  all_cog = all_cog[srt]

      
  
  return
end
