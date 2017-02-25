;+ 
; NAME:
; sdss_compare
;    Version 1.1
;
; PURPOSE:
;  Given the structure containing DLA info, combine the values to
; create a final quality value.
;   
; CALLING SEQUENCE:
;   sdss_compare, strct, FSTRCT=
;
; INPUTS:
;  strct -- QALstrct with the relevant info
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  FSTRCT=  -- Struct containing the combined values [required?]
;
; COMMENTS:
;
; EXAMPLES:
;   sdss_chkbal, 'sdss_DR1_QAL.fits'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   2003 Written by SHF
;-
;------------------------------------------------------------------------------
pro sdss_compare, tmp, FSTRCT=fstrct

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'sdss_compare, strct, FSTRCT= [v1.1]'
    return
  endif 

ntot = tmp.ndla1+tmp.ndla2
if ntot EQ 0 then return
qtot=fltarr(ntot)
z=fltarr(tmp.ndla1+tmp.ndla2)
if tmp.ndla1 NE 0 then msk=bytarr(tmp.ndla1)
if tmp.ndla2 NE 0 then total=fltarr(tmp.ndla2) 


;;if no metals detected
if tmp.ndla2 EQ 0 then begin 
    for i=0, tmp.ndla1-1 do begin
        
        if tmp.dla_score1[i] GT .6 AND tmp.dla_score1[i] LT .65 then qtot[i]=1.       
        if tmp.dla_score1[i] GE .65 AND tmp.dla_score1[i] LT .675 then qtot[i]=1.5
        if tmp.dla_score1[i] GE .675 AND tmp.dla_score1[i] LT .7 then qtot[i]=3.
        if tmp.dla_score1[i] GE .7 AND tmp.dla_score1[i] LT .725 then qtot[i]=4.
        if tmp.dla_score1[i] GE .725 AND tmp.dla_score1[i] LT .75 then qtot[i]=5.
        if tmp.dla_score1[i] GE .75 AND tmp.dla_score1[i] LT .775 then qtot[i]=6.
        if tmp.dla_score1[i] GE .775 AND tmp.dla_score1[i] LT .8 then qtot[i]=6.5
        if tmp.dla_score1[i] GE .8 AND tmp.dla_score1[i] LT .85 then qtot[i]=7.
        if tmp.dla_score1[i] GE .85 AND tmp.dla_score1[i] LT .9 then qtot[i]=7.5
        if tmp.dla_score1[i] GE .9 AND tmp.dla_score1[i] LT .95 then qtot[i]=7.75
        if tmp.dla_score1[i] GE .95 AND tmp.dla_score1[i] LE 1. then qtot[i]=8.
        
        z[i]=tmp.dla_zabs1[i]
        
    endfor
endif else begin
    
;;loop over metal lines
    for j=0, tmp.nDLA2-1 do begin
        
;;calc. quality of each z2
    
        total[j]=round(tmp.dla_hits[j]/tmp.dla_score2[j])
        
        case round(tmp.dla_hits[j]) of	
            0:stop
            1:stop
            2:begin
                case total[j] of 
                    0:stop
                    1:qtot[j]=3.6
                    2:qtot[j]=3.5
                    3:qtot[j]=3.4
                    4:qtot[j]=3.3
                    5:qtot[j]=3.1
                    6:qtot[j]=2.9
                    7:qtot[j]=2.7
                    8:qtot[j]=2.5
                    9:qtot[j]=2.3
                    10:qtot[j]=2.1
                    11:qtot[j]=1.9
                    12:qtot[j]=1.7
                    else: stop
                endcase
            end
            3:begin
                case total[j] of 
                    1:stop
                    2:qtot[j]=3.1
                    3:qtot[j]=3.0
                    4:qtot[j]=2.9
                    5:qtot[j]=2.8
                    6:qtot[j]=2.7
                    7:qtot[j]=2.4
                    8:qtot[j]=2.2
                    9:qtot[j]=2.0
                    10:qtot[j]=2.8
                    11:qtot[j]=2.6
                    12:qtot[j]=2.4
                    else:stop
                endcase
            end
            4:begin
                case total[j] of
                    1:stop
                    2:stop
                    3:qtot[j]=7.6
                    4:qtot[j]=7.5
                    5:qtot[j]=7.4
                    6:qtot[j]=7.3
                    7:qtot[j]=7.1
                    8:qtot[j]=6.9
                    9:qtot[j]=6.7
                    10:qtot[j]=6.5
                    11:qtot[j]=6.3
                    12:qtot[j]=6.1
                    else:stop
                endcase
            end
            5:begin
                case total[j] of 
                    1:stop
                    2:stop
                    3:stop
                    4:qtot[j]=8.6
                    5:qtot[j]=8.5
                    6:qtot[j]=8.4
                    7:qtot[j]=8.3
                    8:qtot[j]=8.2
                    9:qtot[j]=8.1
                    10:qtot[j]=8.0
                    11:qtot[j]=7.9
                    12:qtot[j]=7.8
                    else:stop
                endcase
            end
            6:begin
                case total[j] of
                    1:stop
                    2:stop
                    3:stop
                    4:stop
                    5:qtot[j]=9.1
                    6:qtot[j]=9.0
                    7:qtot[j]=8.9
                    8:qtot[j]=8.8
                    9:qtot[j]=8.7
                    10:qtot[j]=8.6
                    11:qtot[j]=8.5
                    12:qtot[j]=8.4
                    else:stop
                endcase
            end
            7:begin
                case total[j] of
                    1:stop
                    2:stop
                    3:stop
                    4:stop
                    5:qtot[j]=9.6
                    6:qtot[j]=9.6
                    7:qtot[j]=9.5
                    8:qtot[j]=9.4
                    9:qtot[j]=9.3
                    10:qtot[j]=9.2
                    11:qtot[j]=9.1
                    12:qtot[j]=9.0
                    else:stop
                endcase
            end
            else: qtot[j] = 10.
        endcase
;;stop
        
;;pick out which z1's are very close to this z2
        iclose=where(abs(tmp.dla_zabs2[j]-tmp.dla_zabs1) LT .015)
        
;stop
;;if no close matches are found for this z2
        if iclose[0] EQ -1 then begin
            
;;check for coverage, if no, let be, if yes (but no detection), take 3
;;off from quality of this z2
            if (tmp.DLA_zabs2[j]+1)*1215.6701d $
              GE tmp.start_wave AND tmp.start_wave GT 0 then qtot[j]=(qtot[j]-3.)>0. 
;stop        
        endif else begin
            
;;if close match is found to this z2, calc. total quality of this z         
            
            if tmp.dla_score1[iclose] GT .6 AND tmp.dla_score1[iclose] LT .65 then qtot[j]=qtot[j]+1.       
            if tmp.dla_score1[iclose] GE .65 AND tmp.dla_score1[iclose] LT .675 then qtot[j]=qtot[j]+1.5
            if tmp.dla_score1[iclose] GE .675 AND tmp.dla_score1[iclose] LT .7 then qtot[j]=qtot[j]+3.
            if tmp.dla_score1[iclose] GE .7 AND tmp.dla_score1[iclose] LT .725 then qtot[j]=qtot[j]+4.
            if tmp.dla_score1[iclose] GE .725 AND tmp.dla_score1[iclose] LT .75 then qtot[j]=qtot[j]+5.
            if tmp.dla_score1[iclose] GE .75 AND tmp.dla_score1[iclose] LT .775 then qtot[j]=qtot[j]+6.
            if tmp.dla_score1[iclose] GE .775 AND tmp.dla_score1[iclose] LT .8 then qtot[j]=qtot[j]+6.5
            if tmp.dla_score1[iclose] GE .8 AND tmp.dla_score1[iclose] LT .85 then qtot[j]=qtot[j]+7.
            if tmp.dla_score1[iclose] GE .85 AND tmp.dla_score1[iclose] LT .9 then qtot[j]=qtot[j]+7.5
            if tmp.dla_score1[iclose] GE .9 AND tmp.dla_score1[iclose] LT .95 then qtot[j]=qtot[j]+7.75
            if tmp.dla_score1[iclose] GE .95 AND tmp.dla_score1[iclose] LE 1. then qtot[j]=qtot[j]+8.
            
;;assign high score to z with both
            if tmp.dla_score1[iclose] GT .7 AND tmp.dla_hits[j] GE 2. then qtot[j]=qtot[j]>10.
            
            msk(iclose)=1b
                  
        endelse
        
;;take z from metal test
        z[j]=tmp.dla_zabs2[j]
        
    endfor
    
;stop
    
;;loop over dla's that don't match
    if tmp.ndla1 NE 0 then begin
        badmatch=where(msk EQ 0b, nbad)
        for i=0, nbad-1 do begin
        
;;calc. total quality of these badmatch z's    
        
            if tmp.dla_score1[badmatch[i]] GT .6 AND tmp.dla_score1[badmatch[i]] LT .65 then qtot[tmp.ndla2+i]=1.       
            if tmp.dla_score1[badmatch[i]] GE .65 AND tmp.dla_score1[badmatch[i]] LT .675 then qtot[tmp.ndla2+i]=1.5
            if tmp.dla_score1[badmatch[i]] GE .675 AND tmp.dla_score1[badmatch[i]] LT .7 then qtot[tmp.ndla2+i]=3.
            if tmp.dla_score1[badmatch[i]] GE .7 AND tmp.dla_score1[badmatch[i]] LT .725 then qtot[tmp.ndla2+i]=4.
            if tmp.dla_score1[badmatch[i]] GE .725 AND tmp.dla_score1[badmatch[i]] LT .75 then qtot[tmp.ndla2+i]=5.
            if tmp.dla_score1[badmatch[i]] GE .75 AND tmp.dla_score1[badmatch[i]] LT .775 then qtot[tmp.ndla2+i]=6.
            if tmp.dla_score1[badmatch[i]] GE .775 AND tmp.dla_score1[badmatch[i]] LT .8 then qtot[tmp.ndla2+i]=6.5
            if tmp.dla_score1[badmatch[i]] GE .8 AND tmp.dla_score1[badmatch[i]] LT .85 then qtot[tmp.ndla2+i]=7.
            if tmp.dla_score1[badmatch[i]] GE .85 AND tmp.dla_score1[badmatch[i]] LT .9 then qtot[tmp.ndla2+i]=7.5
            if tmp.dla_score1[badmatch[i]] GE .9 AND tmp.dla_score1[badmatch[i]] LT .95 then qtot[tmp.ndla2+i]=7.75
            if tmp.dla_score1[badmatch[i]] GE .95 AND tmp.dla_score1[badmatch[i]] LE 1. then qtot[tmp.ndla2+i]=8.
            
;;take z from dla test
            z[tmp.ndla2+i]=tmp.dla_zabs1[badmatch[i]]
        
        endfor
    endif
endelse

print, 'z, quality, hits, score2, score1'
;printcol, z, qtot, tmp.dla_hits, tmp.dla_score2, tmp.dla_score1

zsub = where(z GT 0.,nz)
fstrct.dla_z[0:nz-1]=z[zsub]
fstrct.dla_quality[0:nz-1]=qtot[zsub]

;stop
end



    
    
    
