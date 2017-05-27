;+
; NAME:
;
; COMBINE_ARCS
;
; PURPOSE:
;
;   Combine multiple arcs together, dealing with saturated regions
;   semi-intelligently.  This routine might absorb the arc reading as well.
;
; CATEGORY:
;
;  Calibration
;
; CALLING SEQUENCE:
;
;  cleanarc=combine_arcs(arcarray,arcsatarray,dirtyarc=dirtyarc)
;
; INPUTS:
;
;  arcarray - (NX x NY x #of arcs) array, containing each arc as a 2d
;             image
;  arcsatarray - (NX x NY x #of arcs) array, 1 on bad pixels, 0 otherwise
;
; OUTPUTS:
;
;  cleanarc - combined arc, with saturated areas in only one arc masked out
;  dirtyarc - combined arc, with no masking
; PROCEDURE:
;
;  Take multiple arcs, some of which may be badly saturated.  Mask out
;  the regions with nonzero flux on the saturated lines, and replace
;  with median value in that arc.  Then, sum the arcs together.
;
; MODIFICATION HISTORY:
;   created jan 2003mar12
;-


function combine_arcs,arcarray,arcsatarray,dirtyarc=dirtyarc

  s = size(arcarray, /dim)
  nx=s[0]
  ny=s[1]
  narcs=s[2]

; total # of saturated pixels at a given point
  totsat=total(arcsatarray,3)
  totsat=smooth(totsat,7)

  dirtyarc=total(arcarray,3)

; to get the median quickly
  sample=randomu(seed,20000)*float(nx)*ny

  for arc=0,narcs-1 do begin
      medarc=median((arcarray[*,*,arc])[sample])
      mask=bytarr(nx,ny)

; need the dilation - sometimes satmask is spotty
      dilsat=dilate(arcsatarray[*,*,arc],bytarr(21,9)+1)
      dilsat=(dilsat AND arcarray[*,*,arc] gt 4.5E4) OR arcsatarray[*,*,arc]

; work one column at a time
      for i=0,nx-1 do begin
;          whbad=where(dilsat[i,*] AND totsat[i,*] LT narcs-0.5,badct)
          whbad=where(dilsat[i,*],badct)
          if badct gt 0 then begin
              fluxderiv=deriv(arcarray[i,*,arc])
              ; choose at most 2 points on each saturated line
              whuse=whbad[where(shift(whbad,-1)-whbad NE 1 $
                                OR shift(whbad,1)-whbad NE (-1),nuse)]
              for j=0,nuse-1 do begin
; bottom end of saturated line
                  edge1=max(where(fluxderiv le 0. AND $
                                  lindgen(ny) lt whuse[j] AND $
                                  dilsat[i,*] eq 0))
                  if edge1 eq -1 then edge1 = 0
; top end of saturated line
                  edge2=min(where(fluxderiv ge 0. AND $
                                  lindgen(ny) gt whuse[j] AND $
                                  dilsat[i,*] eq 0))
                  if edge2 eq -1 then edge2 =(ny-1)
                  mask[i,edge1:edge2]=1

              endfor
         endif     
     endfor     

     arcarray[*,*,arc]=arcarray[*,*,arc]*(mask eq 0)+medarc*mask
;     arcarray[*,*,arc]=djs_maskinterp(arcarray[*,*,arc],mask,iaxis=1)


  endfor

; sum arcs.  could make this a mean, if we wanted.
  combarc=total(arcarray,3)


return,combarc
end
