;+ 
; NAME:
; mtch_wcs   
;   Version 1.1
;
; PURPOSE:
;    I dont recall what this program is for at all.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Jan-2004 Written by JXP (revised x_register)
;-
;------------------------------------------------------------------------------

pro mtch_wcs, lst1, lst2
;mtch_wcs, 'Lists/wcsU.list', 'Lists/wcsB.list'

 ;; Input
 readcol, lst1, fil1, FORMAT='A', /silent
 readcol, lst2, fil2, FORMAT='A', /silent

 ;; Loop on files
 nfil = n_elements(fil1)

 for q=0L,nfil-1 do begin
     ;; Print
     print, 'Reading ', fil1[q]
     ;; Open the img
;     h0 = headfits(fil1[q],EXTEN=0)
     i0 = xmrdfits(fil1[q], 0, h0, /silent)

     i1 = xmrdfits(fil1[q], 1, h1, /silent)
     i2 = xmrdfits(fil1[q], 2, h2, /silent)
     i3 = xmrdfits(fil1[q], 3, h3, /silent)
     i4 = xmrdfits(fil1[q], 4, h4, /silent)
     i5 = xmrdfits(fil1[q], 5, h5, /silent)
     i6 = xmrdfits(fil1[q], 6, h6, /silent)
     i7 = xmrdfits(fil1[q], 7, h7, /silent)
     i8 = xmrdfits(fil1[q], 8, h8, /silent)

     ;; Loop on the new headers
     for k=1,8 do begin
         case k of 
             1: hd = h1
             2: hd = h2
             3: hd = h3
             4: hd = h4
             5: hd = h5
             6: hd = h6
             7: hd = h7
             8: hd = h8
         endcase
         hnew = headfits(fil2[q], EXTEN=k)

         ;; CRVAL1
         crval1 = sxpar(hnew, 'CRVAL1')
         sxaddpar, hd, 'CRVAL1', crval1
         ;; CRVAL2
         crval2 = sxpar(hnew, 'CRVAL2')
         sxaddpar, hd, 'CRVAL2', crval2
         ;; CRPIX1
         crpix1 = sxpar(hnew, 'CRPIX1')
         sxaddpar, hd, 'CRPIX1', crpix1
         ;; CRPIX2
         crpix2 = sxpar(hnew, 'CRPIX2')
         sxaddpar, hd, 'CRPIX2', crpix2
         ;; CD1_1
         cd1_1 = sxpar(hnew, 'CD1_1')
         sxaddpar, hd, 'CD1_1', cd1_1
         ;; CD2_1
         cd2_1 = sxpar(hnew, 'CD2_1')
         sxaddpar, hd, 'CD2_1', cd2_1
         ;; CD1_2
         cd1_2 = sxpar(hnew, 'CD1_2')
         sxaddpar, hd, 'CD1_2', cd1_2
         ;; CD2_2
         cd2_2 = sxpar(hnew, 'CD2_2')
         sxaddpar, hd, 'CD2_2', cd2_2

         ;; WAT0_001
         wat0_001 = sxpar(hnew, 'WAT0_001')
         sxaddpar, hd, 'WAT0_001', wat0_001
         ;; WAT1_001
         wat1_001 = sxpar(hnew, 'WAT1_001')
         sxaddpar, hd, 'WAT1_001', wat1_001
         ;; WAT1_002
         wat1_002 = sxpar(hnew, 'WAT1_002')
         sxaddpar, hd, 'WAT1_002', wat1_002
         ;; WAT1_003
         wat1_003 = sxpar(hnew, 'WAT1_003')
         sxaddpar, hd, 'WAT1_003', wat1_003
         ;; WAT1_004
         wat1_004 = sxpar(hnew, 'WAT1_004')
         sxaddpar, hd, 'WAT1_004', wat1_004
         ;; WAT1_005
         wat1_005 = sxpar(hnew, 'WAT1_005')
         sxaddpar, hd, 'WAT1_005', wat1_005
         ;; WAT2_001
         wat2_001 = sxpar(hnew, 'WAT2_001')
         sxaddpar, hd, 'WAT2_001', wat2_001
         ;; WAT2_002
         wat2_002 = sxpar(hnew, 'WAT2_002')
         sxaddpar, hd, 'WAT2_002', wat2_002
         ;; WAT2_003
         wat2_003 = sxpar(hnew, 'WAT2_003')
         sxaddpar, hd, 'WAT2_003', wat2_003
         ;; WAT2_004
         wat2_004 = sxpar(hnew, 'WAT2_004')
         sxaddpar, hd, 'WAT2_004', wat2_004
         ;; WAT2_005
         wat2_005 = sxpar(hnew, 'WAT2_005')
         sxaddpar, hd, 'WAT2_005', wat2_005

         case k of 
             1: h1 = hd
             2: h2 = hd
             3: h3 = hd
             4: h4 = hd
             5: h5 = hd
             6: h6 = hd
             7: h7 = hd
             8: h8 = hd
         endcase
     endfor

     print, 'Writing ', fil1[q]
     ;; Output
     mwrfits, i0, fil1[q], h0, /create
     mwrfits, i1, fil1[q], h1
     mwrfits, i2, fil1[q], h2
     mwrfits, i3, fil1[q], h3
     mwrfits, i4, fil1[q], h4
     mwrfits, i5, fil1[q], h5
     mwrfits, i6, fil1[q], h6
     mwrfits, i7, fil1[q], h7
     mwrfits, i8, fil1[q], h8
     
 endfor

 return
end	
