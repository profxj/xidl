;+
; NAME:
;   lls_mkclm_idlin
;
; PURPOSE:
;   Generates a .clm file from a ID lines .fits file
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;  NEWDAT=  -- Use this data file instead of original
;  ZEM= -- QSO emission redshift.  If set, only include lines outside
;          the forest
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   ;
; REVISION HISTORY:
;  Generated on 10 July 2013 by JXP
;  Some tweaks (clobber, changed paths) by MF, August 2013 
;-
;------------------------------------------------------------------------------
pro lls_mkclm_idlin, idfil, outfil, DATFIL=datfil, ROOT=root, DATFLG=datflg, $
                     NHI=nhi, ZEM=zem, clobber=clobber 

  ;if not keyword_set(ROOT) then root = getenv('LLSTREE')+'/Abund/'
  if not keyword_set(ROOT) then root = ''
  if not keyword_set(DATFLG) then datflg = 1 ;; Default to HIRES
  if not keyword_set(NHI) then NHI = [17.0, 0.2] ;; Default to HIRES
  
  if keyword_set(ZEM) then wv_forest = (1+zem)*1215.6701 else wv_forest = 1.

  ;; Check to see if outfil exists
  fil = file_search(root+outfil, count=nfil)
  if nfil GT 0 and ~keyword_set(clobber) then begin
     print, 'lls_mkclm_idlin: Will not over-write ', outfil
     print, 'lls_mkclm_idlin: Use /CLOBBER if you wish'
     return
  endif
     
  ;; Read
  idlin = xmrdfits(idfil,1,/silen)
  ntrans = n_elements(idlin)

  close, /all

  ;; ;;;;;;;;;;;;;;;;;;;;
  ;; Open for writing
  openw, 11, root+outfil
  
  ;;deal with subfolders if present 
  pfold=strpos(outfil,'/')
  if(pfold ge 0) then trimmedout=strmid(outfil,pfold+1,strlen(outfil)-pfold) else trimmedout=outfil

  ;; Data files
  prso = strsplit(trimmedout,'.', /extrac)
  

  printf, 11, prso[0]
  printf, 11, strtrim(datflg,2)  
  if not keyword_set(DATFIL) then datfil = idlin[0].datfil
  printf, 11, strtrim(datfil,2)  
  
  ;; zabs
  printf, 11, string(idlin[0].zabs, format='(f7.5)')

  ;; Output
  printf, 11, 'Abund/Tables/'+prso[0]+'.'+prso[1]+'.ion'

  ;; NHI
  printf, 11, string(NHI[0], format='(f5.2)')+', '+$
          string(NHI[1], format='(f4.2)')

  ;; Dummy line
  printf, 11, '0'

  ;; Loop on transitions
  for qq=0L,ntrans-1 do begin
     ;; Forest?
     if idlin[qq].wrest*(1+idlin[qq].zabs) LT wv_forest then continue

     ;;  NG?
     if idlin[qq].flg_anly EQ 0 then continue
     
     ;; Limit flag
     case idlin[qq].flg_limit of
        1: flg = 0
        2: flg = 2
        3: flg = 4
        else: stop
     endcase
     printf, 11, strtrim(flg,2)

     ;; Wave, Integration, and instrument
     printf, 11, idlin[qq].wrest, ',', round(idlin[qq].dv[0]), ',', round(idlin[qq].dv[1]), ',', datflg, $
             format='(f9.4,a1,1x,i5,a1,1x,i5,a1,1x,i3)'
  endfor

  close, /all

  print, 'lls_mkclm_idlin:  Remember to add your .clm files to the SVN!!'
  return
end
