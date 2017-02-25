pwd
pwd;; civ_search.pro
;; 21 Feb 2008
@measureSN

pro civ_search,siiv=siiv
mlss_dir = strtrim(getenv('MLSS_DIR'),2)
if mlss_dir eq '' then mlss_dir = strtrim(getenv('HOME'),2)
mlss_dir = mlss_dir + '/'
cd,mlss_dir

if keyword_set(siiv) then dblt_name = 'SiIV' $
else dblt_name = 'CIV'

;;;;;;;;;;;;;;;;;
;; HOME DIRECTORY
;;;;;;;;;;;;;;;;;
;homedir = mlss_dir+'analysis/'
;homedir = mlss_dir+'analysis/lowsnrspec/'
homedir = mlss_dir+'analysis/nocivspec/'
cd,homedir

;;;;;;;;;;;;;;;;;
;; LIST QUASARS
;;;;;;;;;;;;;;;;;
;;exclude CVS and random directories (e.g. gallin)
spawn,'ls -d [3ABD-Z]* CS*',qso 
nqso = n_elements(qso)

readcol,mlss_dir+'summary/targets_radec_fmt.lst',targ,ra,dec,zqso,$
  delimiter=',',format='a,a,a,f'
targ = strtrim(targ,2)
rad = 1.
decd = 1.
x_radec,ra,dec,rad,decd         ;convert to degrees


;;;;;;;;;;;;;;;;;
;; LOOP QUASARS
;;;;;;;;;;;;;;;;;
for ii=0,nqso-1 do begin 
    ;;;;;;;;;;;;;;;;;
    ;; CHANGE DIRECTORY; CREATE NAMES
    ;;;;;;;;;;;;;;;;;
    curdir = homedir+qso[ii]+'/' 
    cd,curdir 	

    instr = curdir+'lists/'+qso[ii]+'instr.lst' 
    test = file_search(instr,count=ntest) 
    if ntest eq 0 then continue 

    auto = curdir+qso[ii]+'abslin_auto.fits' 
    root = curdir+'search/'+qso[ii] 

    ;;;;;;;;;;;;;;;;;
    ;; Search for features
    ;;;;;;;;;;;;;;;;;
;;    search_all,instr,auto,3.,/calcb,name=qso[ii],/linlst,/table,root=root 


    mtch = where(targ eq qso[ii],nmtch) 
    
    if nmtch ne 1 then stop,'civ_search: no zqso for ',qso[ii] 
    ;;;;;;;;;;;;;;;;;
    ;; FIND, PLOT ALL DOUBLETS/PAIRS
    ;;;;;;;;;;;;;;;;;
    cd,curdir+'search/dblts/' 
    civroot = curdir+'search/dblts/'+qso[ii]+strlowcase(dblt_name)
    ;; Search 5000 km/s beyond zqso 
    civ_find,auto,siiv=siiv,root=civroot,$
      zlim=[-1000./3e5,5000./3e5*(1+zqso[mtch[0]])+zqso[mtch[0]]],$
      instrfil=instr 
    civfits = civroot+'.fits' 

    ;; Verify something found
    test = file_search(civfits,count=ntest) 
    if ntest eq 0 then continue 

    ;; Load QSO info
    civstrct = xmrdfits(civfits,1,/silent) 
    civstrct.qso = targ[mtch[0]]  
    civstrct.zqso = zqso[mtch[0]] 
    civstrct.ra = rad[mtch[0]] 
    civstrct.dec = decd[mtch[0]] 
    mwrfits,civstrct,civfits,/create,/silent 

    civ_calcewn,civfits 
    
    flg_sat = 1.
    flg_mtch = 1.
    civ_aodm,civfits,/log,nbin=3 ,flg_sat=flg_sat,flg_mtch=flg_mtch ;	,/view 
    
    civ_flag,civfits,dvlim=10. 
    
    civ_velplt,civfits,/outline,/label,concat=civroot+'_all.ps'

    civ_prntcivcand,civfits,civroot+'.tab'
;    stop 

    ;;;;;;;;;;;;;;;;;
    ;; HOME DIRECTORY
    ;;;;;;;;;;;;;;;;;
    cd,homedir 

endfor                          ;loop nqso

end
