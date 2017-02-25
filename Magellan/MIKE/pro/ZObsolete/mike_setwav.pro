;+ 
; NAME:
; mike_setwav   
;     Version 1.1
;
; PURPOSE:
;    Process and combine arc files  
;
; CALLING SEQUENCE:
;   
;  mike_setwav, mike, slit, /CLOBBER
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Integer defining setup
;
; RETURNS:
;
; OUTPUTS:
;  One processed, combined Arc Image (e.g. Arcs/Arc_ECH##.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_setwav, mike, 0.5
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_setwav, side, ordr_str, WV_ENDS=wv_ends, CDELT=cdelt, $
                 BIN=bin, ALL_CRVAL=all_crval, NPIX=npix

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_mkarc, side, ordrstr, WV_ENDS=, CDELT= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( BIN ) then bin = 1.d

; Loop on orders
  nordr = n_elements(ordr_str)
  if arg_present( ALL_CRVAL ) then all_crval = dblarr(nordr)
  wv_ends = dblarr(nordr,2)
  for ii=0L,nordr-1 do begin

      if side EQ 2 then begin
          case ordr_str[ii].order of
              76: wv_ends[ii,*] = [4430.d, 4590.]
              75: wv_ends[ii,*] = [4490.d, 4650.]
              74: wv_ends[ii,*] = [4555.d, 4710.]
              73: wv_ends[ii,*] = [4620.d, 4770.]
              72: wv_ends[ii,*] = [4690.d, 4850.]
              71: wv_ends[ii,*] = [4755.d, 4915.]
              70: wv_ends[ii,*] = [4820.d, 4990.]
              69: wv_ends[ii,*] = [4895.d, 5060.]
              68: wv_ends[ii,*] = [4970.d, 5130.]
              67: wv_ends[ii,*] = [5040.d, 5210.]
              66: wv_ends[ii,*] = [5110.d, 5290.]
              65: wv_ends[ii,*] = [5195.d, 5370.]
              64: wv_ends[ii,*] = [5280.d, 5450.]
              63: wv_ends[ii,*] = [5360.d, 5540.]
              62: wv_ends[ii,*] = [5450.d, 5630.]
              61: wv_ends[ii,*] = [5540.d, 5720.]
              60: wv_ends[ii,*] = [5630.d, 5820.]
              59: wv_ends[ii,*] = [5730.d, 5915.]
              58: wv_ends[ii,*] = [5830.d, 6015.]
              57: wv_ends[ii,*] = [5930.d, 6125.]
              56: wv_ends[ii,*] = [6040.d, 6235.]
              55: wv_ends[ii,*] = [6150.d, 6355.]
              54: wv_ends[ii,*] = [6260.d, 6465.]
              53: wv_ends[ii,*] = [6380.d, 6580.]
              52: wv_ends[ii,*] = [6505.d, 6710.]
              51: wv_ends[ii,*] = [6630.d, 6850.]
              50: wv_ends[ii,*] = [6760.d, 6980.]
              49: wv_ends[ii,*] = [6910.d, 7125.]
              48: wv_ends[ii,*] = [7050.d, 7275.]
              47: wv_ends[ii,*] = [7200.d, 7600.]
              46: wv_ends[ii,*] = [7350.d, 7730.]
              else: stop
          endcase
          cdelt = 2.234d / 2.9979E5 / alog(10.d) * bin  
          npix = 5000L / bin
          if arg_present( ALL_CRVAL ) then begin
              tot_wave = 10^(alog10(3000.d) + dindgen(200000L)*cdelt)
              mn = min(abs(wv_ends[ii,0]-tot_wave),imn)
              all_crval[ii] = alog10(tot_wave[imn])
          endif
      endif else begin ; blue side
          case ordr_str[ii].order of
              75: wv_ends[ii,*] = [4698.d, 4800.]
              76: wv_ends[ii,*] = [4635.d, 4740.]
              77: wv_ends[ii,*] = [4575.d, 4675.]
              78: wv_ends[ii,*] = [4515.d, 4620.]
              79: wv_ends[ii,*] = [4460.d, 4560.]
              80: wv_ends[ii,*] = [4400.d, 4505.]
              81: wv_ends[ii,*] = [4348.d, 4450.]
              82: wv_ends[ii,*] = [4295.d, 4395.]
              83: wv_ends[ii,*] = [4240.d, 4340.]
              84: wv_ends[ii,*] = [4190.d, 4290.]
              85: wv_ends[ii,*] = [4140.d, 4236.]
              86: wv_ends[ii,*] = [4095.d, 4190.]
              87: wv_ends[ii,*] = [4048.d, 4140.]
              88: wv_ends[ii,*] = [4000.d, 4095.]
              89: wv_ends[ii,*] = [3955.d, 4050.]
              90: wv_ends[ii,*] = [3910.d, 4005.]
              91: wv_ends[ii,*] = [3868.d, 3960.]
              92: wv_ends[ii,*] = [3825.d, 3915.]
              93: wv_ends[ii,*] = [3785.d, 3875.]
              94: wv_ends[ii,*] = [3745.d, 3830.]
              95: wv_ends[ii,*] = [3705.d, 3790.]
              96: wv_ends[ii,*] = [3665.d, 3750.]
              97: wv_ends[ii,*] = [3628.d, 3712.]
              98: wv_ends[ii,*] = [3590.d, 3675.]
              99: wv_ends[ii,*] = [3555.d, 3638.]
             100: wv_ends[ii,*] = [3520.d, 3600.]
             101: wv_ends[ii,*] = [3485.d, 3565.]
             102: wv_ends[ii,*] = [3450.d, 3530.]
             103: wv_ends[ii,*] = [3415.d, 3495.]
             104: wv_ends[ii,*] = [3380.d, 3462.]
             105: wv_ends[ii,*] = [3350.d, 3430.]
             106: wv_ends[ii,*] = [3318.d, 3400.]
             107: wv_ends[ii,*] = [3285.d, 3366.]
             108: wv_ends[ii,*] = [3255.d, 3340.]
             109: wv_ends[ii,*] = [3225.d, 3310.]
              else: stop
          endcase
          cdelt = 2.234d / 2.9979E5 / alog(10.d) * bin  
          npix = 5000L / bin
          if arg_present( ALL_CRVAL ) then begin
              tot_wave = 10^(alog10(3000.d) + dindgen(200000L)*cdelt)
              mn = min(abs(wv_ends[ii,0]-tot_wave),imn)
              all_crval[ii] = alog10(tot_wave[imn])
          endif
      endelse
  endfor
;  ALL DONE
  print, 'mike_setwav: All Done! '
  return


end

