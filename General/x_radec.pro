;+ 
; NAME:
; x_radec
;
; PURPOSE:
;    Turns RA, DEC in XX:XX:XX.X -DD:DD:DD.D format to decimal deg
;      or VICE VERSA
;
; CALLING SEQUENCE:
;   
;   x_radec, ra, dec, rad, decd, /ARCS, /FLIP
;
; INPUTS:
;   ra, dec    - RA and DEC in in RR:RR:RR.R -DD:DD:DD.D format 
;                 Colons are required as separators
;   rad, decd  - RA and DEC in decimal degrees (double)
;
; RETURNS:
;
; OUTPUTS:
;   ra, dec    - RA and DEC in in RR:RR:RR.R -DD:DD:DD.D format 
;                 Colons are required as separators
;   rad, decd  - RA and DEC in decimal degrees (double)
;
; OPTIONAL KEYWORDS:
;   ARCS - Outputs in arcseconds
;   FLIP - Gives RA and DEC from decimal RA,DEC
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_radec, '21:12:23.1', '-13:13:22.2', rad, decd
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   04-Aug-2001 Written by JXP
;   20-Nov-2003 Added + sign for /flip
;-
;------------------------------------------------------------------------------

pro x_radec, ra, dec, rad, decd, ARCS=arcs, FLIP=flip

  if (N_params() LT 4) then begin 
    print,'Syntax - ' + $
             'x_radec, ra, dec, rad, decd, /ARCS, /FLIP (v1.0)'
    return
  endif 

  if not keyword_set( FLIP ) then begin

;  Remove blanks
      
      rastr = strcompress(ra, /REMOVE_ALL)
      decstr = strcompress(dec, /REMOVE_ALL)
      
      
;  Parse on colons and convert to float
      
      rav = strsplit(rastr, ':', /extract)
      decv = strsplit(decstr, ':', /extract)
      
                                ; Deal with negative DEC
      flg_neg = 0
      if strmid(decv[0],0,1) EQ '-' then flg_neg = 1
      
      rah = double(rav[0])
      ram = double(rav[1])
      ras = double(rav[2])
      
      decdeg = double(decv[0])
      decm = double(decv[1])
      decs = double(decv[2])
      
;  Convert to deg
      
      rad = (360./24.)*(rah + ram/60. + ras/3600.)
      decd = abs(decdeg) + decm/60. + decs/3600.
      
      if flg_neg EQ 1 then decd = (-1)*decd
      
      if keyword_set( ARCS ) then begin
          rad = rad * 3600
          decd = decd * 3600
      endif
;;;;;;;;;; PASS BACK STRING (FLIP)  ;;;;;;;;;;;;;;;; 
  endif else begin

;  RA
      rah = fix(rad*24./360.)
      ram = fix((rad - double(rah)*360./24.)*24.*60./360.)
      ras = (rad - double(rah)*360./24. - double(ram)*360/(24.*60.))*24*60*60/360.
      if rah GE 10 then srah = strtrim(rah,2) $
        else srah = '0'+strtrim(rah,2)
      if ram GE 10 then sram = strtrim(ram,2) $
        else sram = '0'+strtrim(ram,2)
      if ras GE 10 then sras = strmid(strtrim(ras,2),0,5) $
        else sras = '0'+strmid(strtrim(ras,2),0,4)

      ra = strjoin([srah,':',sram,':',sras])

; DEC

      if decd LT 0. then begin
          if abs(decd) LT 10 then sdec = '-0'+strtrim(fix(abs(decd)),2) $
          else sdec = strtrim(fix(decd),2) 
      endif else begin
          if abs(decd) LT 10 then sdec = '+0'+strtrim(fix(decd),2) $
          else sdec = '+'+strtrim(fix(decd),2) 
      endelse

      adecd = abs(decd)
      decm = fix( (adecd - fix(adecd))*60.)
      decs = (adecd - fix(adecd) - double(decm)/60.)*60*60.

      if decm GE 10 then sdecm = strtrim(decm,2) $
        else sdecm = '0'+strtrim(decm,2)
;      if decs GE 10 then sdecs = strtrim(decs,2) $
;        else sdecs = '0'+strtrim(decs,2)
      if decs GE 10 then sdecs = strmid(strtrim(decs,2),0,4) $
        else sdecs = '0'+strmid(strtrim(decs,2),0,3)

      dec = strjoin([sdec,':',sdecm,':',sdecs])

  endelse

end

