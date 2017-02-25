;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_flag.pro                
; Author: Kathy Cooksey                      Date: 18 Feb 2008
; Project: Metal-line System Survey with Jason X. Prochaska
; Description: Evaluate likelihood of lines
; Input: 
;   strct_fil -- 
; Optional:
;
; Output: 
;   strct_fil
; Calls:
;   dblt_retrieve
;   civ_aodm_mtch
; Example:
; History:
;   18 Feb 2008  created by KLC
;    4 Mar 2008  AODM flag calc wrong
;    5 Mar 2008  Also enable subtraction of flag
;   12 Mar 2008  Be sure flag correct; fix EW ratio bug
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_aodm

pro civ_flag_summ,flgsumm,siiv=siiv

  flgsumm = strarr(9)

  ;; Flag legend:
  if keyword_set(siiv) then begin
     flgsumm[8] = '256 - detect SiIV 1393 at 3sigma EW'
     flgsumm[7] = '128 - detect SiIV 1402 at 3sigma EW '
     flgsumm[6] = '64 - EW ratio of 1393/1402 (1-sR <= R <= 2+sR) when both lines >= 1sigma'
  endif else begin
     flgsumm[8] = '256 - detect CIV 1548 at 3sigma EW'
     flgsumm[7] = '128 - detect CIV 1550 at 3sigma EW '
     flgsumm[6] = '64 - EW ratio of 1548/1550 (1-sR <= R <= 2+sR) when both lines >= 1sigma'
  endelse 
  flgsumm[5] = '32 - dv <= 10 km/s'
  flgsumm[4] = '16 - Feature found at location of Lya w/EW >= 3sigma'
  flgsumm[3] = '8 - 1548 outside Lya forest'
  flgsumm[2] = '4 - 1548 outside H2 forest (1138.867 H2 B0-0P(7))'
  flgsumm[1] = '2 - tau_AOD per element in agreement for >=68.3%'
  flgsumm[0] = '1 - Features found at location of other (non_Lya) lines w/EW >= 3sigma'
end                             ;civ_flag_summ


function civ_flag_prsflag,flg,nflg
allflg = intarr(30) ;large enough to contain number of flags
if flg eq 0 then begin
   nflg = 0
   return,allflg
endif 

  ;; Seperate binary flags
  lgv = alog(flg)/alog(2)
  qq = fix(lgv+0.00001)
  nn = flg-2^qq
  lft = where(nn ne 0,nlft,complement=rgt,ncomplement=nrgt)
  if nlft eq 0 then done = 1 else done = 0
  if nrgt ne 0 then nn[rgt] = 0
  allflg[0] = 2^qq
  cnt = 1

  while not done do begin
     lgv = alog(nn[lft])/alog(2)
     qq = fix(lgv+0.00001)
     allflg[cnt] = 2^qq
     cnt = cnt + 1
     nn[lft] = nn[lft]-2^qq
     lft = where(nn ne 0,nlft,complement=rgt,ncomplement=nrgt)
     if nlft eq 0 then done = 1 
     if nrgt ne 0 then nn[rgt] = 0
  endwhile

  allflg = allflg[0:cnt-1]
  nflg = n_elements(allflg)
  return,allflg
end                             ;civ_flag_prsflag


function civ_flag_strflag,flg,bybyte=bybyte

flgarr = civ_flag_prsflag(flg,nflg)
flg_str = ''

if nflg eq 0 then return,flg_str ;no information

for ii=0,nflg-1 do begin
    if keyword_set(bybyte) then $
      flg_str = flg_str + strtrim(flgarr[ii],2)+'; ' $
    else case flgarr[ii] of
        256: flg_str = flg_str + '3sig 1548; '
        128: flg_str = flg_str + '3sig 1550; '
        64: flg_str = flg_str + 'EW ratio agree; '
        32: flg_str = flg_str + 'dv<=10km/s; '
        16: flg_str = flg_str + 'Lya too; '
        8: flg_str = flg_str + 'Outside Lya forest; '
        4: flg_str = flg_str + 'Outside H2 forest; '
        2: flg_str = flg_str + '68.3% profile agree; '
        1: flg_str = flg_str + 'Other lines too; '
        else: stop,'civ_flag_strflag: problematic flag'
    endcase
endfor

flg_str = strmid(flg_str,0,strlen(flg_str)-2)
return,flg_str
end                             ;civ_flag_strflag

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro civ_flag,strct_fil,dvlim=dvlim,flgdescript=flgdescript,$
             savfil=savfil,hdf5_fil=hdf5_fil,dblt_name=dblt_name

  if not keyword_set(dvlim) then dvlim = 10. ;km/s  

  ;;;;;;;;;;;;;;;;;;;;
  ;; Flag legend:
  ;; 256 - detect CIV 1548 at 3sigma EW 
  flg1548 = 256
  ;; 128 - detect CIV 1550 at 3sigma EW 
  flg1550 = 128
  ;; 64 - EW ratio of 1548/1550 (1-sR <= R <= 2+sR) 
  ;;      when both lines >= 1sigma
  rtoflg = 64
  ;; 32 - dv <= 10 km/s (incorp zsig and/or add dvflg2)
  dvflg = 32
  ;; 16 - Feature found at location of Lya w/EW >= 3sigma
  lyaflg = 16
  ;; 8 - 1548 outside Lya forest 
  lyaforflg = 8
  ;; 4 - 1548 outside H2 forest (1138.867 H2 B0-0P(7))
  h2forflg = 4
  ;; 2 - tau_AOD per element in agreement for >=68.3%
  tauflg = 2
  ;; 1 - Features found at location of other (non_Lya) lines w/EW >= 3sigma
  sysflg = 1
  ;;;;;;;;;;;;;;;;;;;;
  flgtot = total(flg1548 or flg1550 or rtoflg or dvflg or lyaflg or $
                 lyaforflg or h2forflg or tauflg or sysflg)

  if keyword_set(flgdescript) then begin
     flgdescript = ['Flags -- '+strtrim(flg1548,2)+$
                    ': EW(1548)/sigEW(1548) >= 3; '+$
                    strtrim(flg1550,2)+': ibid for 1550',$
                    '      '+strtrim(rtoflg,2)+': 1-sR <= R <= 2+sR; ', $
                    '      '+strtrim(dvflg,2)+': dv <= '+strtrim(dvlim,2)+$
                    ' km/s; '+$
                    strtrim(lyaflg,2)+': potential Lya w/significance >= 3; ', $
                    '      '+strtrim(lyaforflg,2)+': outside Lya forest; '+$
                    strtrim(h2forflg,2)+': ibid for H2; ',$
                    '      '+strtrim(tauflg,2)+': >=68.3% AOD profile similar', $
                    '      '+strtrim(sysflg,2)+$
                    ': potential non-Lya lines, sig>=3; ' $
                   ]
     return
  endif                         ;/flgdescript

  ;;;;;;;;;;;;;;;;
  ;; Read in/process structure
  ;;;;;;;;;;;;;;;;
  if size(strct_fil,/type) eq 7 then begin ;must read in
     if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
           stop,'civ_flag: restored structure not of expected name'
        strct = civcand
     endif else $
        strct = xmrdfits(strct_fil, 1, /silent) 
  endif else strct = strct_fil
  nstrct = n_elements(strct)

  if not keyword_set(dblt_name) then begin
     ;; Assume all first filled element define structure focus
     gd = where(strtrim(strct[0].ion,2) ne '')
     prs = strsplit(strct[0].ion[gd[0]],/extract,count=nprs) 
     dblt_name = prs[0]
  endif 

  civinfo = dblt_retrieve(dblt_name)
  lyainfo = dblt_retrieve('Lya') 


  for ii=0,nstrct-1 do begin
     nwflg = 0                  ; restart every time

     ;;;;;;;;;;;;;;;;
     ;; Find CIV 
     ;;;;;;;;;;;;;;;;
     civ = where(stregex(strct[ii].ion,dblt_name,/boolean),nciv,complement=sys)
     if nciv ne 2 then stop,'civ_flag: no '+dblt_name+' doublet to evaluate'
     srt = sort(strct[ii].wrest[civ])
     civ = civ[srt]
     
     ;;;;;;;;;;;;;;;;
     ;; Verify whether even EW error use
     ;;;;;;;;;;;;;;;;
     test = where(strct[ii].sigew[civ] eq 0.,ntest)

     ;;;;;;;;;;;;;;;;
     ;; Line significance
     ;;;;;;;;;;;;;;;;
     if ntest eq nciv and keyword_set(hdf5_fil) then begin
        nsig1548 = 1.e4         ; just large
        nsig1550 = 1.e4 
     endif else begin
        nsig1548 = strct[ii].ew[civ[0]]/strct[ii].sigew[civ[0]]
        nsig1550 = strct[ii].ew[civ[1]]/strct[ii].sigew[civ[1]]
     endelse 
     if nsig1548 ge 3. then nwflg = nwflg + flg1548
     if nsig1550 ge 3. then nwflg = nwflg + flg1550

     ;;;;;;;;;;;;;;;;
     ;; EW Ratio
     ;;;;;;;;;;;;;;;;
     ewrto = strct[ii].ew[civ[0]]/strct[ii].ew[civ[1]]
     if ntest eq nciv and keyword_set(hdf5_fil) then $
        sigewrto = 0. $
     else sigewrto = abs(ewrto)*sqrt(1./nsig1548^2+1./nsig1550^2)
     if nsig1548 ge 1. and nsig1550 ge 1. and $ ; so can measure EW
        ewrto ge 1.-sigewrto and ewrto le 2.+sigewrto then $
           nwflg = nwflg + rtoflg

     ;;;;;;;;;;;;;;;;
     ;; dv limit
     ;;;;;;;;;;;;;;;;
     dv = 2.998e5*(strct[ii].zabs[civ[1]]-strct[ii].zabs[civ[0]])/$
          (1+strct[ii].zabs[civ[0]])
     if abs(dv) le dvlim then nwflg = nwflg + dvflg 

     ;;;;;;;;;;;;;;;;
     ;; Lya > 3sigma EW
     ;;;;;;;;;;;;;;;;
     if ntest eq nciv and keyword_set(hdf5_fil) then $
        gd = where(abs(strct[ii].wrest[sys]-lyainfo.wvI) lt 1e-4,ngd) $
     else $
        gd = where(abs(strct[ii].wrest[sys]-lyainfo.wvI) lt 1e-4 and $
                   strct[ii].ew[sys]/strct[ii].sigew[sys] ge 3.,ngd)
     if ngd ne 0 then nwflg = nwflg + lyaflg

     ;;;;;;;;;;;;;;;;
     ;; Lya Forest
     ;;;;;;;;;;;;;;;;
     if strct[ii].wrest[civ[0]]*(1+strct[ii].zabs[civ[0]]) gt $
        lyainfo.wvI*(1+strct[ii].zqso) then $
           nwflg = nwflg + lyaforflg 

     ;;;;;;;;;;;;;;;;
     ;; H2 Forest (H2 B0-0P(7))
     ;;;;;;;;;;;;;;;;
     if strct[ii].wrest[civ[0]]*(1+strct[ii].zabs[civ[0]]) gt $
        1138.867*(1+strct[ii].zqso) then nwflg = nwflg + h2forflg
     
     ;;;;;;;;;;;;;;;;
     ;; Profile shape
     ;;;;;;;;;;;;;;;;
     ngd = civ_aodm_mtch(strct[ii])
     if ngd ge 0.683 then nwflg = nwflg + tauflg

     ;;;;;;;;;;;;;;;;
     ;; Other (non-Lya) lines found with > 3sigma EW
     ;;;;;;;;;;;;;;;;
     if ntest eq nciv and keyword_set(hdf5_fil) then $
        gd = where(strct[ii].wrest[sys] gt 0. and $
                   abs(strct[ii].wrest[sys]-lyainfo.wvI) gt 1e-4,ngd) $
     else $
        gd = where(strct[ii].wrest[sys] gt 0. and $
                   abs(strct[ii].wrest[sys]-lyainfo.wvI) gt 1e-4 and $
                   strct[ii].ew[sys]/strct[ii].sigew[sys] ge 3.,ngd)
     if ngd ne 0 then nwflg = nwflg + sysflg

     ;; Finish
     strct[ii].flg_sys[0] = nwflg
     tmp = civ_flag_prsflag(strct[ii].flg_sys[0],nflg)
     strct[ii].flg_sys[1] = nflg

  endfor                        ;loop nstrct

  ;;;;;;;;;;;;;;;;
  ;; Save
  ;;;;;;;;;;;;;;;;
  if size(strct_fil,/type) eq 7 then begin
     spawn,'cp '+strct_fil+' '+strtrim(strct_fil,2)+'.bkp'
     if keyword_set(savfil) then begin
        civcand = strct
        save,civcand,filename=strct_fil
     endif else mwrfits, strct, strct_fil, /create,/silent
     print,'civ_flag: overwrote ',strct_fil
  endif else strct_fil = strct

  test = max(strct.flg_sys[0],min=mn)
  if test gt flgtot or mn lt 0 then stop,'civ_flag: error assigning flags'
end
