;+ 
; NAME:
; lls_struct
;  V1.2
;
; PURPOSE:
;    Given a list of LLS base files, fill up a LLS structure.  
;
; CALLING SEQUENCE:
;   
;   lls_struct, lls, [list]
;
; INPUTS:
;  [list] -- List of DLA base files.  [default:
;            '/u/xavier/LLS/Lists/tot_lls.lst']
;
; RETURNS:
;
; OUTPUTS:
;  lls  -- IDL DLA structure
;
; OPTIONAL KEYWORDS:
;  /ION - Input ionic column densities
;  /NOELM - Supress elemental and ion structures (saves memory)
;  /NORELM - Supress inputting Elemental [X/H] values
;  /NOHIS -- Suppress HI error in [X/H] values
;  /EW    -- Fill up the ion arrays with EW values instead of column
;            densities.  This reads the .EW files instead of .ion
;  /FINE  -- This fills up arrays for fine-structure states
;            (e.g. SiII*).  By default CII* is already considered 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_dlalst, dla, '/u/xavier/DLA/Lists/tot_dla.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   08-Mar-2005 Written by JXP (based on parse_dlalst)
;-
;------------------------------------------------------------------------------
pro lls_struct, supstrc, list, NOELM=noelm, ION=ion, NOHIS=nohis, FILE=file,$
                NOSDSS=nosdss, NORELM=norelm, ROOT=root, NOSYS=nosys, EW=ew, $
                FINE=fine, VPFIT=vpfit, NOFILLELM=NOFILLELM, DEBUG=debug

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lls_struct, struct, [list], /NOELM, /ION, /NOHIS ROOT=, /EW (v1.2)' 
    return
  endif 

  if keyword_set(NOELM) and (not keyword_set(ION)) then llss = {llsnoestruct} else llss = {llsstruct} 
  if keyword_set(NOSYS) then llss = {llsnosys}
  
  ; Parse the Base File
  if not keyword_set( FILE ) then begin
      readcol, list, listnms, format='A', /silent
      nlls = n_elements(listnms)
  endif else begin
      listnms = list
      nlls = n_elements(list)
   endelse

  ;; Skip those with ##
  keep = where( strmatch( strmid(listnms,0,1), '#') EQ 0, nlls)
  listnms = listnms[keep]

  ;; Add ROOT
  if keyword_set(ROOT) then listnms = root+listnms else root = ''

  ;; Create Big struct
  supstrc = replicate(llss,nlls)

  dumc = ''
 
  close, 1
  for i=0,nlls-1 do begin
      fil = listnms[i]
      print, 'lls_struct: Reading ', fil
      nline = file_lines(fil)
      donelin = 0L
      openr, 1, fil
      supstrc[i].llsfil = strtrim(fil,2)
      readf, 1, dumc
      ;; QSO stuff
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].qso = prs[0]
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].qso_ra = prs[0]
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].qso_dec = prs[0]
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].qso_zem = float(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].flg_mag = long(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].qso_mag = float(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].srvy = long(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].srvy_mag = float(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      ndumc = n_elements(prs)
      for j=0,ndumc-3 do begin
          supstrc[i].ref = supstrc[i].ref+' '+prs[j]
      endfor
      ;; SDSS stuff
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      ndumc = n_elements(prs)
      supstrc[i].sdss_plate = long(prs[0])
      supstrc[i].sdss_mjd = long(prs[2])
      supstrc[i].sdss_fibid = long(prs[1])
      readf, 1, dumc
      ;; Redshift, NHI, NH
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].zabs = double(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].nhi = float(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].signhi[0] = float(prs[0])
      supstrc[i].signhi[1] = float(prs[1])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].nh = float(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].nhsig[0] = float(prs[0])
      supstrc[i].nhsig[1] = float(prs[1])
      ;; Kinematics
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].vmn = double(prs[0])
      supstrc[i].vmx = double(prs[1])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].fdelv = double(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].fmm = double(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].fedg = double(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].ftpk = double(prs[0])
      ;; [M/H] (NH weighted)
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].flg_mh = long(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].mhave = double(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      if n_elements(prs) EQ 4 then $
         supstrc[i].mhsig = double(prs[0:1]) $
      else $
         supstrc[i].mhsig[*] = double(prs[0])

      ;; D/H
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].flg_dh = long(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].dh = double(prs[0])
      ;; Abund (in components)
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      supstrc[i].nsys = long(prs[0])
      readf, 1, dumc
      prs = strsplit(dumc,' ', /extract)
      ;; Cloudy
      if prs[0] EQ '!' then $
         supstrc[i].cldyfil = 'dum.fil' $
      else supstrc[i].cldyfil = prs[0]
      if keyword_set(NOSYS) then begin
          close, 1
          continue
      endif
      ;; Subsystems
      for j=0,supstrc[i].nsys-1 do begin
          ;; Subsystem Name
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].name = prs[0]
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].zabs = double(prs[0])
          ;; NHI
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].nhi = double(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].nhisig = double(prs[0:1])
          ;; NH
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].nh = double(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].nhsig = double(prs[0:1])
          ;; x
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].logx = double(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].sig_logx = double(prs[0:1])
          ;; b
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].b = double(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].bsig = double(prs[0])
          ;;
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].abndfil = prs[0]
          ;; U
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].U = float(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].Usig = double(prs[0:1])
          ;; Kin
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].flg_low = long(prs[0])
          ;; [a/H]
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].flg_alpha = long(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].alphah = double(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].sig_alphah = double(prs[0])
          if n_elements(prs) EQ 5 then $
             supstrc[i].systems[j].sig_alphah = double(prs[0:1]) $
          else $
             supstrc[i].systems[j].sig_alphah = double(prs[0])
          ;; [Fe/H]
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].flgFe = long(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          supstrc[i].systems[j].feh = double(prs[0])
          readf, 1, dumc
          prs = strsplit(dumc,' ', /extract)
          if n_elements(prs) EQ 5 then $
             supstrc[i].systems[j].sig_feh = double(prs[0:1]) $
          else $
             supstrc[i].systems[j].sig_feh = double(prs[0])
          ;; VPFIT file
          readf, 1,  format='(a60)', dumc
          supstrc[i].systems[j].vpfil = strtrim(dumc,2)
      endfor
      close, 1
  endfor

  ;;  Read in Elm info
  if not keyword_set ( NOELM ) and not keyword_set(NOFILLELM) and $ 
    not keyword_set(NOSYS) then begin
      for nn=0L,n_elements(supstrc)-1 do begin
          for qq=0L,supstrc[nn].nsys-1 do begin
             if keyword_set(DEBUG) then print, supstrc[nn].llsfil
              lls_fillelm, supstrc, nn, qq, ROOT=root
          endfor
      endfor
  endif

  ;;  Read in Ion info
  if keyword_set ( ION ) then begin
     for nn=0L,n_elements(supstrc)-1 do begin
        ;; Fill up the sub-systems
        for qq=0L,supstrc[nn].nsys-1 do begin
           lls_fillion, supstrc, nn, qq, ROOT=root
        endfor
        ;; Total for the full LLS
        lls_sumionsys, supstrc, nn 
     endfor
  endif

  ;;  Read in EW info
  if keyword_set ( EW ) then begin
      for nn=0L,n_elements(supstrc)-1 do begin
          for qq=0L,supstrc[nn].nsys-1 do begin
              lls_fillew, supstrc, nn, qq, ROOT=root
          endfor
      endfor
  endif

  ;; Fine structure
  if keyword_set ( FINE ) then begin
      for nn=0L,n_elements(supstrc)-1 do begin
          for qq=0L,supstrc[nn].nsys-1 do begin
              tmp = supstrc[nn].systems[qq]
              x_fineabnd, tmp
              ;; Save
              supstrc[nn].systems[qq] = tmp
          endfor
      endfor
  endif

  ;; Vpfit
  if keyword_set(VPFIT) then begin
      print, 'lls_struct:  Using VPFIT values if they exist'
      for nn=0L,n_elements(supstrc)-1 do begin
          for qq=0L,supstrc[nn].nsys-1 do begin
              ;; Read in
              if strlen(supstrc[nn].systems[qq].vpfil) GT 0 then $
                prs_vpfit, supstrc[nn].systems[qq].vpfil, vpstr $
              else continue ; vpstr = {newabslinstrct}
              ;; Find ions
              uniion = vpstr[ uniq(vpstr.ion, sort(vpstr.ion)) ].ion
              nion = n_elements(uniion)
              ;; Sum em up
              for jj=0L,nion-1 do begin
                  ipos = strpos(uniion[jj],'*')
                  if ipos NE (-1) then begin
                      ;; CII*
                      if strmid(uniion[jj],0,1) EQ 'C' then jion=6
                  endif else jion = 0
                  ;; Sum up vpfit
                  vp_clm = 0.
                  vp_sig = 0.
                  vpmt = where(strmatch(vpstr.ion,uniion[jj]),nvp)
                  for tt=0L,nvp-1 do begin
                      x_logclm, Ni, sig, vpstr[vpmt[tt]].N, $
                                vpstr[vpmt[tt]].nsig, /rev
                      vp_clm = vp_clm + Ni
                      vp_sig = sqrt(vp_sig^2 + sig^2)
                  endfor
                  x_logclm, vp_clm, vp_sig, logn, logs
                  ;; Fill in
                  getion, uniion[jj], ion, /INM, Z=Z
                  supstrc[nn].systems[qq].ion[Z].state[ion,jion].flgclm = 1
                  supstrc[nn].systems[qq].ion[Z].state[ion,jion].clm = logn
                  supstrc[nn].systems[qq].ion[Z].state[ion,jion].sigclm = logs
              endfor
          endfor
      endfor
   endif

  ;; Cut bad ones (mainly those skipped with a # sign)
  gd = where(supstrc.zabs GT 0.)
  supstrc = supstrc[gd]


end
