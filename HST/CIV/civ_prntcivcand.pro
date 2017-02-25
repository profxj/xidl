;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_prntcivcand.pro               
; Author: Kathy Cooksey                      Date: 20 Feb 2008
; Project: HST Metal-line System Survey (CIV Focus) with 
;          Xavier Prochaska
; Description: Print CIV structure
; Input: 
;   strct_fil - name of structure file (FITS or IDL sav) or
;               structure itself
;   outfil - name of output table
; Optional:
;   latex - outfil will be latex formatted
;   savfil - strct_fil is IDL save file of civcand
; Output: 
;   outfil
; Example:
; History:
;   20 Feb 2008  created by KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function civ_prntcivcand_srtwrest,civstrct,num
ncivstrct = n_elements(civstrct.wrest)
mask = replicate(0,ncivstrct)

;; Assume first element defines structure focus
gd = where(strtrim(civstrct[0].ion,2) ne '')
prs = strsplit(civstrct[0].ion[gd[0]],/extract,count=nprs) 
dblt_name = prs[0]

;; CIV or whichever doublet
gd = where(stregex(civstrct.ion,dblt_name,/boolean),ngd)
if ngd ne 2 then $
  stop,'civ_prntcivcand_srtwrest: '+dblt_name+' doublet not present'
srt = sort(civstrct.wrest[gd])
indx = gd[srt]
mask[gd] = 1

;; HI
gd = where(stregex(civstrct.ion,'HI',/boolean),ngd)
if ngd ne 0 then begin
    srt = reverse(sort(civstrct.wrest[gd]))
    indx = [indx,gd[srt]]
    mask[gd] = 1
endif

;; Rest
gd = where(civstrct.wrest gt 0 and mask eq 0,ngd)
if ngd ne 0 then begin
    srt = sort(civstrct.wrest[gd])
    indx = [indx,gd[srt]]
    mask[gd] = 1
endif 

num = n_elements(indx)
return,indx

end                             ;civ_prntcivcand_srtwrest

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro civ_prntcivcand,strct_fil,outfil,latex=latex,wiki=wiki,savfil=savfil,$
                    asis=asis,_extra=extra

if size(strct_fil,/type) eq 7 then begin
    if keyword_set(savfil) then begin
        restore,strct_fil
        if not keyword_set(civcand) then $
          stop,'civ_prntcivcand: save file structure not named civcand'
    endif else civcand = xmrdfits(strct_fil,1,/silent)
endif else civcand = strct_fil

dvdr = '-'
for ii=0,86 do dvdr = dvdr+'-'

openw,1,outfil

iqso = uniq(civcand.qso,sort(civcand.qso))
nqso = n_elements(iqso)

if not keyword_set(asis) then begin
    for qq=0,nqso-1 do begin

        ;;;;;;;;;;;;;;;;
        ;; Normal ASCII
        ;;;;;;;;;;;;;;;;
        if not keyword_set(latex) and not keyword_set(wiki) then begin
            ;; QSO (zqso = )
            printf,1,dvdr,format='(a0)' ;horizontal line
            printf,1,civcand[iqso[qq]].qso,civcand[iqso[qq]].zqso,$
              format='(8x,a'+string(strlen(strtrim(civcand[iqso[qq]].qso,2)))+$
              ',1x,"(zqso = ",f7.5,")")'
            printf,1,dvdr,format='(a0)' ;horizontal line

            ;; Ion, wobs, wr, zabs, zsig, EWr, sigEWr, logN, siglogN, flgcolm, instr
            printf,1,'Ion','Wobs','Wrest','zabs','zsig','EWr','sigEWr','logN','siglogN',$
              'flgN','Instr',$
              format='(a10,2x,a9,1x,a9,1x,a9,1x,a9,2x,a5,1x,a5,2x,a5,1x,a4,1x,a4,2x,a5)'
            printf,1,dvdr,format='(a0)' ;horizontal line

            ;; System (based on CIV)
            sys = where(civcand[iqso[qq]].qso eq civcand.qso,nsys)
            if nsys eq 0 then stop,'civ_prntcivcand: probelmatic lack of matching QSO'
            sys = civcand[sys]
            for ss=0,nsys-1 do begin
                ;; Lines associated with system (CIV)
                lin = civ_prntcivcand_srtwrest(sys[ss],nlin)
                if nlin eq 0 then $
                  stop,'civ_prntcivcand: problematic lack of matching system'

                ;; zsys = #.#; flg_sys = # (totnum)
                printf,1,sys[ss].zabs[lin[0]],sys[ss].flg_sys[0],sys[ss].flg_sys[1],$
                  format='(8x,"zsys = ",f8.5,"; flg_sys = ",i4," (",i2,")")'
                printf,1,dvdr,format='(a0)' ;horizontal line
                writecol,outfil,strtrim(sys[ss].ion[lin],2),$
                  sys[ss].wrest[lin]*(1+sys[ss].zabs[lin]),sys[ss].wrest[lin],$
                  sys[ss].zabs[lin],sys[ss].zsig[lin],$
                  round(sys[ss].ew[lin]),round(sys[ss].sigew[lin]),$
                  sys[ss].ncolm[lin],sys[ss].signcolm[lin],sys[ss].flg_colm[lin],$
                  sys[ss].instr[lin],filnum=1,$
                  fmt='(a10,2x,f9.4,1x,f9.4,1x,f9.5,1x,f9.5,2x,i5,1x,i5,2x,f5.2,1x,f4.2,1x,i4,2x,i5)'
                printf,1,dvdr,format='(a0)' ;horizontal line
                
            endfor              ;loop nsys
        endif                   ;normal ASCII table


        ;;;;;;;;;;;;
        ;; WIKI
        ;;;;;;;;;;;;
        if keyword_set(wiki) then begin
            npglim = 30
            svsys = where(civcand[iqso[qq]].qso eq civcand.qso,nsys)
            if nsys eq 0 then stop,'civ_prntcivcand: probelmatic lack of matching QSO'

            npage = ceil(nsys/double(npglim)) ;cap number of tables
            for pp=0,npage-1 do begin
                ;; QSO (zqso = )
                printf,1,civcand[iqso[qq]].qso,civcand[iqso[qq]].zqso,$
                  format='("======  ",8x,a'+string(strlen(strtrim(civcand[iqso[qq]].qso,2)))+$
                  ',1x,"(zqso = ",f7.5,")  ======")'
                printf,1,''


                ;; System (based on CIV)
                imx = pp+npglim < nsys
                sys = civcand[svsys[pp:imx-1]]
                for ss=0,n_elements(sys)-1 do begin
                    ;; Lines associated with system (CIV)
                    lin = civ_prntcivcand_srtwrest(sys[ss],nlin)
                    if nlin eq 0 then $
                      stop,'civ_prntcivcand: problematic lack of matching system'

                    printf,1,''
                    ;; zsys = #.#; flg_sys = # (totnum)  CENTER
                    printf,1,sys[ss].zabs[lin[0]],sys[ss].flg_sys[0],sys[ss].flg_sys[1],$
                      format='("====",8x,"zsys = ",f8.5,"; flg_sys = ",i4," (",i2,")  ====")'

                    ;; Ion, wobs, wr, zabs, zsig, EWr, sigEWr, logN, siglogN, flgcolm, instr
                    printf,1,'Ion','W<sub>obs</sub>','W<sub>rest</sub>',$
                      'z<sub>abs</sub>','EW<sub>r</sub>','sig(EW<sub>r</sub>)','logN','sig(logN)',$
                      'flg<sub>N</sub>','Instr',$
                      format='("^  ",a,"  ^  ",2x,a,"  ^  ",1x,a,"  ^  ",1x,a,"  ^  ",2x,a,"  ^  ",1x,a,"  ^  ",2x,a,"  ^  ",1x,a,"  ^  ",1x,a,"  ^  ",2x,a,"  ^  ")'
                    writecol,outfil,strtrim(sys[ss].ion[lin],2),$
                      sys[ss].wrest[lin]*(1+sys[ss].zabs[lin]),sys[ss].wrest[lin],$
                      sys[ss].zabs[lin],$
                      round(sys[ss].ew[lin]),round(sys[ss].sigew[lin]),$
                      sys[ss].ncolm[lin],sys[ss].signcolm[lin],sys[ss].flg_colm[lin],$
                      sys[ss].instr[lin],filnum=1,$
                      fmt='("|",a10,"|",2x,f9.4,"|",1x,f9.4,"|",1x,f9.5,"|",2x,i5,"|",1x,i5,"|",2x,f5.2,"|",1x,f4.2,"|",1x,i4,"|",2x,i5,"|")'
                    
                endfor          ;loop nsys
                printf,1,''
                printf,1,''
                printf,1,''
            endfor              ;loop npage
        endif                   ;wiki table

    endfor                      ;loop nqso


endif else begin                ;not asis
;;;;;;;;;;;
;; As is
;;;;;;;;;;;
    nciv = n_elements(civcand)
    for ii=0,nciv-1 do begin
        svqso = ''              ;comparison
        if not keyword_set(latex) and not keyword_set(wiki) then begin
            if strtrim(civcand[ii].qso,2) ne svqso then begin
                ;; QSO (zqso = )
                printf,1,dvdr,format='(a0)' ;horizontal line
                printf,1,civcand[ii].qso,civcand[ii].zqso,$
                  format='(8x,a'+string(strlen(strtrim(civcand[ii].qso,2)))+$
                  ',1x,"(zqso = ",f7.5,")")'
                printf,1,dvdr,format='(a0)' ;horizontal line
                svqso = strtrim(civcand[ii].qso,2) ;save
            endif               ;compare qso

            ;; Ion, wobs, wr, zabs, zsig, EWr, sigEWr, logN, siglogN, flgcolm, instr
            printf,1,'Ion','Wobs','Wrest','zabs','zsig','EWr','sigEWr',$
              'logN','siglogN','flgN','Instr',$
              format='(a10,2x,a9,1x,a9,1x,a9,1x,a9,2x,a5,1x,a5,2x,a5,1x,a4,1x,a4,2x,a5)'
            printf,1,dvdr,format='(a0)' ;horizontal line

            ;; Lines associated with system (CIV)
            lin = civ_prntcivcand_srtwrest(civcand[ii],nlin)
            if nlin eq 0 then $
              stop,'civ_prntcivcand: problematic lack of matching system'

            ;; zsys = #.#; flg_sys = # (totnum)
            printf,1,civcand[ii].zabs[lin[0]],civcand[ii].flg_sys[0],$
              civcand[ii].flg_sys[1],$
              format='(8x,"zsys = ",f8.5,"; flg_sys = ",i4," (",i2,")")'
            printf,1,dvdr,format='(a0)' ;horizontal line
            writecol,outfil,strtrim(civcand[ii].ion[lin],2),$
              civcand[ii].wrest[lin]*(1+civcand[ii].zabs[lin]),civcand[ii].wrest[lin],$
              civcand[ii].zabs[lin],civcand[ii].zsig[lin],$
              round(civcand[ii].ew[lin]),round(civcand[ii].sigew[lin]),$
              civcand[ii].ncolm[lin],civcand[ii].signcolm[lin],$
              civcand[ii].flg_colm[lin],civcand[ii].instr[lin],filnum=1,$
              fmt='(a10,2x,f9.4,1x,f9.4,1x,f9.5,1x,f9.5,2x,i5,1x,i5,2x,f5.2,1x,f4.2,1x,i4,2x,i5)'
            printf,1,dvdr,format='(a0)' ;horizontal line
            
        endif                   ;normal ASCII table
    endfor                      ;loop nciv


    ;;;;;;;;;;;;
    ;; WIKI
    ;;;;;;;;;;;;
    if keyword_set(wiki) then begin
        npglim = 30
        npage = ceil(nciv/double(npglim)) ;cap number of tables
        for pp=0,npage-1 do begin
            svqso = ''
            ;; System (based on CIV)
            imx = (pp+1)*npglim < nciv
            sys = civcand[pp*npglim:imx-1]
            for ss=0,n_elements(sys)-1 do begin
                if strtrim(sys[ss].qso,2) ne svqso then begin
                    ;; QSO (zqso = )
                    printf,1,sys[ss].qso,sys[ss].zqso,$
                      format='("======  ",8x,a'+$
                      string(strlen(strtrim(sys[ss].qso,2)))+$
                      ',1x,"(zqso = ",f7.5,")  ======")'
                    printf,1,''
                    svqso = strtrim(sys[ss].qso,2)
                endif           ;compare svqso


                ;; Lines associated with system (CIV)
                lin = civ_prntcivcand_srtwrest(sys[ss],nlin)
                if nlin eq 0 then $
                  stop,'civ_prntcivcand: problematic lack of matching system'

                printf,1,''
                ;; zsys = #.#; flg_sys = # (totnum)  CENTER
                printf,1,sys[ss].zabs[lin[0]],sys[ss].flg_sys[0],$
                  sys[ss].flg_sys[1],$
                  format='("====",8x,"zsys = ",f8.5,"; flg_sys = ",i4," (",i2,")  ====")'

                ;; Ion, wobs, wr, zabs, zsig, EWr, sigEWr, logN,
                ;; siglogN, flgcolm, instr
                printf,1,'Ion','W<sub>obs</sub>','W<sub>rest</sub>',$
                  'z<sub>abs</sub>','EW<sub>r</sub>','sig(EW<sub>r</sub>)',$
                  'logN','sig(logN)','flg<sub>N</sub>','Instr',$
                  format='("^  ",a,"  ^  ",2x,a,"  ^  ",1x,a,"  ^  ",1x,a,"  ^  ",2x,a,"  ^  ",1x,a,"  ^  ",2x,a,"  ^  ",1x,a,"  ^  ",1x,a,"  ^  ",2x,a,"  ^  ")'
                writecol,outfil,strtrim(sys[ss].ion[lin],2),$
                  sys[ss].wrest[lin]*(1+sys[ss].zabs[lin]),sys[ss].wrest[lin],$
                  sys[ss].zabs[lin],$
                  round(sys[ss].ew[lin]),round(sys[ss].sigew[lin]),$
                  sys[ss].ncolm[lin],sys[ss].signcolm[lin],sys[ss].flg_colm[lin],$
                  sys[ss].instr[lin],filnum=1,$
                  fmt='("|",a10,"|",2x,f9.4,"|",1x,f9.4,"|",1x,f9.5,"|",2x,i5,"|",1x,i5,"|",2x,f5.2,"|",1x,f4.2,"|",1x,i4,"|",2x,i5,"|")'
                
            endfor              ;loop nsys
            printf,1,''
            printf,1,dvdr
            printf,1,'     Page ',strtrim(pp+2,2)
            printf,1,dvdr
            printf,1,''
            printf,1,''
        endfor                  ;loop npage
    endif                       ;wiki table
endelse                         ;end asis
close,1
print,'civ_prntcivcand: created ',outfil

end
