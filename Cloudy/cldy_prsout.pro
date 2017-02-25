;+
; NAME:
; cldy_prsout
;  V1.1
;
; PURPOSE:
;    parse a single CLOUDY output file
;
; CALLING SEQUENCE:
;   cldy_prsout, outfil, cldystrct=cldystrct
;
; INPUT:
;    outfil: name of CLOUDY output file
;
; OUTPUT:
;    cldystrct: output cloudy strct
;
; REVISION HISTORY:
;   12-Sep-2005 Written by GEP
;   27-Apr-2006 Modified to work with Cloudy v06.02 syntax
;-
;----------------------------------------------------------------------

pro cldy_prsout, outfil, cldystrct=cldystrct, zabs=zabs, SPEC=spec

    cldystrct={strctcldy}
    nlin = file_lines(outfil)
    close, 1
    openr, 1, outfil

    if not keyword_set(SPEC) then spec = 'HM'

    dumc=''

    cldystrct.flg = 1
    cldystrct.spec = spec

    if keyword_set(zabs) then cldystrct.z = zabs

    tagsdone = 1
    ntags = 0

    indata = 0
    datadone = 0
    flg_uval = 0
    
    for i=0L, nlin-1 do begin
        readf, 1, dumc
        prs = strsplit(dumc,' ', /extract)
        nprs = n_elements(prs)
        if nprs EQ 1 then continue
        if tagsdone then begin
            case prs[1] of
                'CMB': begin
                   if not keyword_set(zabs) then begin
                      cldystrct.z = float(prs[3])
                      ntags = ntags + 1
                   endif 
                end
                'hextra': begin
                    cldystrct.heat = float(prs[2])
                    ntags = ntags+1
                end
                'hden': begin
                    cldystrct.nh = double(prs[2])
                    ntags = ntags+1
                end
                'metals': begin
                    cldystrct.feh = double(prs[2])
                    ntags = ntags+1
                end
                'globule': begin
                    ;; Find the = sign
                    pos1 = strpos(dumc,'=')
                    pos2 = strpos(dumc,',')
                    cldystrct.nh = double(strmid(dumc,pos1+1,pos2-pos1-1))
                end
                'stop': begin
                    case prs[2] of
                        'neutral': begin
                            cldystrct.nhi = double(prs[5])
                            ntags = ntags+1
                        end
                        'column' : begin
                            cldystrct.nhi = double(prs[4])
                            ntags = ntags+1
                        end
                        'temperature': begin
                            ;;Modified by KLC for Cloudy v06.02 syntax
                            if stregex(dumc,'=',/boolean) then $
                              cldystrct.Tval = double(prs[4]) $
                            else begin
                                if stregex(prs[3],'K',/boolean) then $
                                  prs[3] = strmid(prs[3],0,strpos(prs[3],'K'))
                                cldystrct.Tval = double(prs[3])
                            endelse        
                            ;ntags = ntags+1
                        end                             
                        else: 
                    endcase
                end 
                'ionization': begin
                    if prs[2] EQ 'parameter' then begin
                        cldystrct.u = double(prs[4])
                        ntags = ntags+1
                        flg_uval = 1
                    endif
                end 
                else:
             endcase
;            if ntags EQ 5 then tagsdone = 0
        endif
        for j=0L,nprs-1 do begin
            if prs[j] EQ 'Log(U):' and flg_uval NE 1 then begin
               cldystrct.u = double(prs[j+1])
               flg_uval = 1
            endif
            if prs[j] EQ 'Ionisation' then begin
                if prs[j+1] EQ '(over' and prs[j+2] EQ 'radius)' then begin
                    ;cldystrct.X[0,0] = double(prs[1])
                    ;cldystrct.X[1,0] = double(prs[2])
                    cldystrct.X[1,1] = double(prs[1])
                    cldystrct.X[1,2] = double(prs[2])
                    if prs[3] NE '(H2)' then $
                        cldystrct.X[1,3] = double(prs[3])
                        ;cldystrct.X[2,0] = double(prs[3])
                    for k=2,30 do begin
                        readf, 1, dumc
                        ;; Extra line? -- JXP on Oct 28, 2013.  Works
                        ;;                for v10.0
                        if strmatch(strmid(dumc,0,1),'-') OR $ 
                           strmatch(strmid(dumc,1,1),' ') then begin
                           k=k-1
                           continue
                        endif
                        prs = strsplit(dumc,' ', /extract)
                        nprs = n_elements(prs)
                        jnk =''
                        for nn=1,nprs-1 do jnk = jnk+prs[nn]
                        nums = strsplit(jnk,'-',/extract)
                        nnums = n_elements(nums)
                                ;if k GT 23 then begin -- This might
                                ;                         have worked
                                ;                         for old
                                ;                         veresions
                        ;   stop
                        ;    if nnums LT 5 then begin
                        ;        k = k - 1
                        ;        continue
                        ;    endif
                        ;endif
                        for nn=1,nnums do $
                               cldystrct.X[k,nn] = -1.*double(nums[nn-1])
                    endfor
;                    datadone = 1
;                    stop
                    break
                endif
            endif

            ;; H temperature
            if prs[j] EQ 'Temperature' then begin
                if prs[j+1] EQ '(over' and prs[j+2] EQ 'radius)' then begin
                    ;cldystrct.X[0,0] = double(prs[1])
                    ;cldystrct.X[1,0] = double(prs[2])
                    cldystrct.Tgas[0] = double(prs[1])
                    cldystrct.Tgas[1] = double(prs[2])
                    datadone = 1
                    break
                 endif
             endif
         endfor
        if datadone then break
     endfor

    close, 1

end
