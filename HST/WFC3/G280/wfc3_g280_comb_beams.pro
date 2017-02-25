;+ 
; NAME:
; wfc3_g280_comb_beams
;
; PURPOSE:
;  Program that combines the different beams using long_combspec 
;
; CALLING SEQUENCE:
;   wfc3_g280_comb_beams, wfc3_g280_strct, ii
;
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;   ii -- the object in the structure to combine
;
; RETURNS:
;
; OUTPUTS:
;   Update structure with combined beam's flux, wave and sig
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  wfc3_g280_center_direct, wfc3_g280_strct, ii
;
; PROCEDURES CALLED:
;  long_combspec
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

pro wfc3_g280_comb_beams, wfc3_g280_strct, ii

  ;; loop over image + sky

  for jj=0,1 do begin

     ;; calculate the maximum counts
     cnta=wfc3_g280_strct[ii].cnta
     cntc=wfc3_g280_strct[ii].cntc
     cnt=cnta>cntc
     wfc3_g280_strct[ii].cnt=cnt
     
     ;; deal with how the expected input for long_combspec looks like
     influx=dblarr(cnt,2)
     case jj of
        0: begin
           influx(cnt-cnta:cnt-1,0)=reverse(wfc3_g280_strct[ii].skya(0:cnta-1)*1d17)
           influx(cnt-cntc:cnt-1,1)=wfc3_g280_strct[ii].skyc(0:cntc-1)*1d17
        end
        1: begin
           influx(cnt-cnta:cnt-1,0)=reverse(wfc3_g280_strct[ii].fluxa(0:cnta-1)*1d17)
           influx(cnt-cntc:cnt-1,1)=wfc3_g280_strct[ii].fluxc(0:cntc-1)*1d17
        end
        else: stop
     endcase
     
     inivar=dblarr(cnt,2)
     inivar(cnt-cnta:cnt-1,0)=1/reverse((wfc3_g280_strct[ii].flux_siga(0:cnta-1)*1d17)^2)
     inivar(cnt-cntc:cnt-1,1)=1/(wfc3_g280_strct[ii].flux_sigc(0:cntc-1)*1d17)^2
     
     loglam=dblarr(cnt,2)
     loglam(cnt-cnta:cnt-1,0)=reverse(alog10(wfc3_g280_strct[ii].wavea(0:cnta-1)))
     loglam(cnt-cntc:cnt-1,1)=alog10(wfc3_g280_strct[ii].wavec(0:cntc-1))
     
     newloglam=loglam(*,0)
     LONG_COMBSPEC, influx, inivar, loglam, newloglam=newloglam, $
                    newflux=newflux, newivar=newivar ;, sigrej=2.5
     
     case jj of
        0: begin
           wfc3_g280_strct[ii].sky=newflux*1d-17
        end
        1: begin
           wfc3_g280_strct[ii].wave=10d^newloglam
           wfc3_g280_strct[ii].flux=newflux*1d-17
           wfc3_g280_strct[ii].sig=1/sqrt(newivar)*1d-17
        end
        else: stop
     endcase
  endfor
  
end
