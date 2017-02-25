;\tiny
;\begin{verbatim}
function esopar, header, parameter

; ============================================================
; Since sxpar is unable to handle the ESO hierarchical FITS
; keywords, a routine was needed for this purpose.
;  
;

; Daniel Apai
; 2001., AIU Jena
; ============================================================


idx=where(strpos(header,parameter) ne -1,cnt)


idx=idx(0)
   if cnt eq 1 then begin	
        tp=strpos(header[idx],"'")
;	tp=long(tp(0))
	if (tp eq -1) then begin
	   p1=strpos(header[idx],"=")
	   p2=strpos(header[idx],"/", p1+1)
           val_str = strmid(header[idx], p1+2, p2-p1-2)
           IF strcompress(val_str, /rem) EQ 'T' THEN value = 'T' $
           ELSE IF strcompress(val_str, /rem) EQ 'F' THEN value = 'F' $
           ELSE value = double(val_str)
        endif else begin
 	
 	p1=strpos(header[idx],"'")
	;print, "String begins at ", p1
	p2=strpos(header[idx],"'", p1+1)
	value=strmid(header[idx],p1+1, p2-p1-1)   
   	;print, value
	endelse
   endif  else value=-1    

;value=(float(strmid(head(idx),31,13)))[0] $

return, value
end
;\end{verbatim}
;\normalsize
