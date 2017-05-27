;+
; NAME:
;   deimos_badchip
;
; PURPOSE:
;   Generate bad pixel mask (bad columns, vignetting also)
;
; CALLING SEQUENCE:
;   deimos_badchip, path
; 
; INPUTS:
;   path   - path for output file
;
; OUTPUTS:
;    invvar(2k,4k)  = 0 for good regions; otherwise:
;                     +1 for bad pixels/columns
;                     +2 for >50% vignetted pixels
;                     +4 for >5% vignetted pixels (used to indicate
;                     suspect regions; if a vignetting vigflat is
;                     available, all pixels >0% vignetted are marked); written
;                     to file deimos_badmask.fits.Z (in 8 separate
;                     HDU's, with proper WCS mosaic keywords)
;
;
; MODIFICATION HISTORY:
;    18-Jan-2002 MD
;    26-May-2002 MD - flag bad regions in DEIMOS red mosaic
;    05-Jun-2002 DPF added path 
;    11-Jun-2002 DPF - significant changes
;    20-Jul-2002 DPF - another careful look, using derivative images
;    12-Aug-2002 JAN - added vignetted
;    17-Oct-2002 JAN - used quantitative vignetting determinations
;    03-Feb-2003 MD  - extended a few bad columns
;-
pro deimos_badchip, path

  invvar = bytarr(2048, 4096)
  
  r2 = fltarr(1024, 1024)
  x = indgen(1024)
 ; for i=0, 1023 do r2[i, x] = (4096.-x)^2 + (4096.-i)^2 ;radius of circle
 ; r2limit = 4096.^2 + (4096*3./4.)^2 ;vignetted region
 ; vignette = r2 GE r2limit  ;define good region
  invvar1 = invvar
 ; invvar1[0:1023, 0:1023] = vignette*2B

  ; create x,y arrays for vignetted sides
  xind = lindgen(2048, 4096)
  yind =  xind / 2048
  xind =  xind MOD 2048

  ;if NOT keyword_set(path) then message, 'please call with a path!'
  if NOT keyword_set(path) then path = './'

  maskfilename = path+'deimos_badmask.fits'
  vigflatname = getenv('CALIB_DATA')+'/vignettingflat.fits.gz'

	if file_test(vigflatname) eq 0 then novigflat=1 else novigflat=0

  mwrfits, junk, maskfilename, /create ;dummy primary


; NOTE: corner vignetting in 1/4 and 5/8 is pretty symmetric 
;(but NOT 1/5 or 4/8!) 
; eventually, we might want to use the 'average' vignetting of each
; pair (with proper flippings) rather than the below whvig2/whpartvig

  for i=1, 8 do begin
    case i of
      1: begin
            chip1 = invvar
            chip1[241, 3904:4095] = 1B ;bad column
            chip1[1052:1053, 835:4095] = 1B ;bad columns
	    chip1[1049:1056,835:843]=1B
            chip1[1016, 3190:4095] = 1B ;CTE
            chip1[1081, 3027:4095] = 1B ;CTE
            chip1[1171:1172, 482:565] = 1B
	    chip1[1170:1173,482:485]=1B
            chip1[370:377, 542:546] = 1B ;dead region
            chip1[1361:1369, 1788:1793] = 1B ;dead region
; bad pix from superdarks:
            chip1[0:1, *] = 1B
            chip1[*, 4094:4095] = 1B
            chip1[2046:2047, *] = 1B
            chip1[1172, *] = 1B
; column 1922 is slightly (0.25 cts) hot
             chip1[608, 3543] = 1B
            chip1[711:712, 3592] =  1B
              chip1[2043, 3609] =  1B
              chip1[1933, 3550] = 1B
              chip1[1547, 2691] =  1B                
              chip1[276, 2529] =  1B
              chip1[344, 2512] =  1B
              chip1[903, 2438] =  1B
               chip1[1085, 2207] =  1B
              chip1[253, 1931] =  1B
               chip1[1450, 1483] =  1B
                chip1[775, 1649] =  1B
               chip1[795, 1490] =  1B
               chip1[505, 1445] =  1B
               chip1[566, 1304] =  1B
               chip1[169, 1139] =  1B
              chip1[681, 1131] =  1B
             chip1[635, 1039] =  1B
             chip1[1423, 1276] =  1B
             chip1[1505, 1255] =  1B
            chip1[1190, 925] =  1B
            chip1[419, 353] =  1B
            chip1[766, 275] =  1B
            chip1[1615, 86] =  1B
            chip1[1715, 5] =  1B
            chip1[1607, 1415] =  1B
            chip1[778, 1291 ] = 1B     
            chip1[2045, 2504 ] = 1B   
            chip1[1669, 548  ] = 1B 
            chip1[1810, 519] = 1B 
            chip1[1960, 423] = 1B 
            chip1[1860, 334] = 1B
  
;cold pixels
            chip1[1054:1056, 857:858] = 1B
            ;chip1[1409, 275] =  1B
            ;chip1[1454, 108] =  1B
 ;  1415,1098 maybe            
          ; vignetted side
            whvig = where(xind lt ( (( yind/100)*0.8 +13) > 15) < 40)
            ; 50% vignetted corner
            if novigflat then whvig2 = where(yind lt 1224.8654 $
                           -1.2717381*xind $
                           +0.00017275134*xind^2) else begin
                vigflat=mrdfits(vigflatname,i)
                whvig2=where(vigflat lt 0.5)
            endelse
                        
            chip1[whvig] = chip1[whvig] OR 2B
            chip1[whvig2] = chip1[whvig2] OR 2B

            ; define >5% vignetted region
            if novigflat then whpartvig =where(yind  lt  1563.8401 $
                             -1.3760047*xind $
                             +0.00021877510*xind^2) else $
                  whpartvig=where(vigflat lt 1.0)

            chip1[whpartvig] = chip1[whpartvig] OR 4B

            mkhdr, header, chip1

            mosaic_keywords, header, 1
            mwrfits, chip1, maskfilename, header; chip1
            delvarx, chip1
         end
      2: begin
            chip2 = invvar
            chip2[376:387, *] = 1B ;hot columns
        
                                ; 388-390 still hot by 0.5-1 cts/pixels
                                ; (~4 more columns affected at lower
                                ; levels); 383-hot by 4.5; 384 - 2.8;
                                ;  385-2;386-1.3; 387-1
            chip2[330:331, 3323:4095] = 1B ; bad col
	    chip2[327:333,3329:3333]=1B
            chip2[463, *] = 1B ; trap
            chip2[464, 1960:2000] = 1B
            chip2[488:489, *] = 1B ;hot columns
            chip2[232:235, 3838:3840] = 1B ;cold spot
            chip2[1333:1334, *] = 1B ; hot 
            chip2[1373, *] = 1B
            chip2[1549:1551, 424:426] = 1B
            chip2[1783:1787, 3179:3182] = 1B
            chip2[1772:1783, 3644:3648] = 1B ; cold spot
            chip2[1709:1711, 2017:2019] = 1B

; bad pix from superdarks:
            chip2[0:1, *] = 1B
            chip2[*, 4094:4095] = 1B
            chip2[2046:2047, *] = 1B
           
            chip2[1110, 1600:*] = 1B
            chip2[1374, *] = 1B
            chip2[471, 162] = 1B
            chip2[1234, 25] = 1B
           chip2[1603, 139] = 1B
           chip2[1495, 924] = 1B
           chip2[1335, 992:993] = 1B
           chip2[1110:1111, 918:*] = 1B
           chip2[70, 1115] = 1B
           chip2[461, 1108] = 1B
           chip2[1914, 1171] = 1B
           chip2[1332, 1558:1562]= 1B
           chip2[1853, 1937] = 1B
           chip2[1976:1981, 2193:2195] = 1B
           chip2[1437, 2428] = 1B
           chip2[462, 1966] = 1B
           chip2[97:98, 2448:2460] = 1B
           chip2[373, 2351] = 1B
           chip2[1345, 2698] = 1B
           chip2[1916, 3281] = 1B
           chip2[364:375, 3149:3151] = 1B
           chip2[490:491, 3788:3792] = 1B
           chip2[629, 3766:3767] = 1B
           chip2[1645, 3749] = 1B
           chip2[447, 1952] = 1B
           chip2[461:462, 2014:2018] = 1B
           chip2[1194, 3362] = 1B

            mkhdr, header, chip2
            mosaic_keywords, header, 2
            mwrfits, chip2,  maskfilename, header; chip2
            delvarx, chip2
         end     
      3: begin
            chip3 = invvar
            chip3[221, *] = 1B ;hot column
            chip3[366, 2960:4095] = 1B ; CTE
            chip3[940, 1583:4095] = 1B ; bad column
            chip3[1166:1168, 751:780] = 1B ; hot spot
            chip3[1167, *] = 1B ; CTE
            chip3[1169, 3480:3488] = 1B ; hot spot
            chip3[1280, 2169:4095] = 1B ; CTE
            chip3[1301:1302, 2186:4095] = 1B ; CTE
            chip3[1299:1300, 2186:2200] = 1B ; hot spot
            chip3[817:819, *] = 1B ;hot column
            chip3[851, 1200:4095] = 1B ;blocked colum
            chip3[1744:1747, *] = 1B ;bad CTE
            chip3[651:656, 1324:1332] = 1B
            chip3[649:650, 1330:1334] = 1B
            chip3[768:772, 1062:1064] = 1B
            chip3[799:802, 1044:1045] = 1B
            chip3[1537:1542, 782:784] = 1B
            chip3[1858, *] = 1B
            chip3[2005, 827:4095] = 1B
            chip3[2026:2030, 1558:1563] = 1B
            chip3[2039, 1538] = 1B
            chip3[2037, 1588] = 1B
            
; bad pix from superdarks:
            chip3[0:1, *] = 1B
            chip3[*, 4094:4095] = 1B
            chip3[2046:2047, *] = 1B
            chip3[222:228, 3636:3642] = 1B
            chip3[793, 1600:2050] = 1B
            chip3[251, 247] = 1B
            chip3[899, 158:159] = 1B
            chip3[1747:1748, 156:165] = 1B
            chip3[1747:1752, 157:162] = 1B
            chip3[1255:1257, 615] = 1B
            chip3[1421, 543] = 1B
            chip3[1172, 558]= 1B
            chip3[202, 584]= 1B
            chip3[1486, 1030] = 1B
           chip3[478, 1343] = 1B
           chip3[791, 1642] = 1B
           chip3[191, 1505] = 1B
           chip3[151, 1656] = 1B
           chip3[148, 2263] = 1B
           chip3[914, 2231] = 1B
           chip3[1277, 2167] = 1B
           chip3[1293, 2189] = 1B
           chip3[1307, 2189] = 1B
           chip3[504, 2434:2449] = 1B
           chip3[1314, 2908] = 1B
           chip3[1558, 2542] = 1B
           chip3[1864, 2587:2588] = 1B
           chip3[1979:1983, 2623] = 1B
           chip3[477, 3297:3323] = 1B
           chip3[221, 3191] = 1B
           chip3[687, 3364] = 1B
           chip3[634, 3434] = 1B
           chip3[146, 3871] = 1B
           chip3[594, 25] = 1B
;check   
            chip3[819:823, 3219 ] = 1B
            chip3[819:823, 3223] = 1B
            chip3[815:816,2698:2699] = 1B

            chip3[1897:1898, 2575:2576] = 1B
            chip3[1900:1903, 2572:2573] = 1B
; from schapman's data
            chip3[260,1400:*]=1B
            
;check 819:827,1004:1006 - bad pixflat crossing column
            mkhdr, header, chip3
            mosaic_keywords, header, 3
            mwrfits, chip3,  maskfilename, header; chip3
            delvarx, chip3
         end
      4: begin
            chip4 = reverse(invvar1, 1)
            chip4[47, *] = 1B ; hot
            chip4[571, 253:301] = 1B ; trap
            chip4[743, 2925:4095] = 1B ; hot
            chip4[744, *] = 1B ; hot
            chip4[791, *] = 1B ; slightly hot?
            chip4[48, 1878:*] = 1B
            chip4[46, 1881:1887] = 1B
            chip4[188:192, 3775:3779] = 1B
            chip4[345:347, 3895:3900] = 1B
            chip4[356:358, 3923] = 1B
            chip4[1162:1165, 3540:3543] = 1B
            chip4[787:794, 519:522] = 1B
            chip4[689:696, 2211:2215] = 1B
            chip4[1056:1068, 3933:3940] = 1B
            chip4[997, 2134:4095] = 1B
            chip4[998, *] = 1B ; hot col
            chip4[999, *] = 1B
            chip4[1000:1001, 2133:2144] = 1B
; bad pix from superdarks:
            chip4[0:1, *] = 1B
            chip4[*, 4094:4095] = 1B
            chip4[2046:2047, *] = 1B
            chip4[133, 419] = 1B
            chip4[436, 542] = 1B
            chip4[664, 255] = 1B
            chip4[712, 286:306]= 1B
            chip4[1133, 340:364]= 1B
            chip4[1279, 578] = 1B
            chip4[1933, 478:480] = 1B
            chip4[1464:1469, 239:242] = 1B
            chip4[701,2516:2580] = 1B
            chip4[702, 2516:2535] = 1B   
            chip4[790, 517] = 1B
            chip4[2027, 862] = 1B
            chip4[504, 1088] = 1B
            chip4[1541, 1083:1085] = 1B
            chip4[1681, 1389] = 1B
            chip4[590, 1307] = 1B
            chip4[702, 1524] = 1B
            chip4[218, 1880] = 1B
            chip4[368, 1994:1995] = 1B
            chip4[951, 2081] = 1B
            chip4[547, 2448:2454] = 1B
           chip4[546, 2448:2466] = 1B
           chip4[382, 2443] = 1B
           chip4[576, 2714] = 1B
           chip4[418, 2786] = 1B
           chip4[246, 2735] = 1B
           chip4[547, 3016] = 1B
           chip4[745, 2936:*] = 1B
           chip4[746, 2941] = 1B
           chip4[1021, 3052:3066] = 1B
           chip4[1022, 3052:3054] = 1B
           chip4[790, 3252:3300] = 1B
           chip4[792, 3252:3280] = 1B
           chip4[757, 3501] = 1B
           chip4[462, 3515:3530] = 1B
           chip4[463, 3515:*] = 1B
           chip4[464, 3515:3519] = 1B
           chip4[462:467, 3515] = 1B
           chip4[1510:1513,3967:3973] = 1B
           chip4[127, 4005:4006] = 1B
           chip4[251, 4085] = 1B
           chip4[47:49, *] = 1B ; hot
           chip4[743:746, *] = 1B
           chip4[547, 2447:2454] = 1B
           chip4[546, 2447:2466] = 1B
           chip4[548, 2447:2452] = 1B
           chip4[462, 3515:3524] = 1B
           chip4[463, 3515:3556] = 1B
           chip4[464, 3515:3519] = 1B
           chip4[997:1001, 2133:2150] = 1B
           chip4[1513, 3974:*] = 1B
           chip4[1510, 3974:4000] = 1B
           chip4[1510:1513, 3974:3979] = 1B
           chip4[1541, 1083:1085] = 1B
           chip4[1681, 1385] = 1B
           chip4[1681, 1389] = 1B
           chip4[1667, 1371] = 1B
           chip4[1279, 578] = 1B
           chip4[1933, 478:480] = 1B

           chip4[571,250:*] = 1B

        ; vignetted edge
            whvig = where(xind gt 1968.+(41.-yind/100.)^2/39)
            chip4[whvig] = chip4[whvig] OR 2B


          ; 50% vignetted corner
            if novigflat then whvig2 = where(yind lt -367.51087 $
                           +0.25707244*xind $
                           +0.00024931130*xind^2) else begin
                vigflat=mrdfits(vigflatname,i)
                whvig2=where(vigflat lt 0.5)
            endelse
                        
            chip4[whvig2] = chip4[whvig2] OR 2B

            ; define >5% vignetted region
            if novigflat then whpartvig =where(yind lt -65.497610 $
                             -0.048433023*xind $
                             +0.00044756528*xind^2) else $
                  whpartvig=where(vigflat lt 1.)
            chip4[whpartvig] = chip4[whpartvig] OR 4B

            mkhdr, header, chip4
            mosaic_keywords, header, 4
            mwrfits, chip4,  maskfilename, header; chip4
            delvarx, chip4
         end

      5: begin
            chip5 = reverse(invvar1, 2)
            chip5[22:23, 1180:1480] = 1B ; hot region
            chip5[25:26, *] = 1B ; hot
            chip5[127:130, *] = 1B ; hot
	    chip5[131,0:800] = 1B ;hot		
            chip5[255, 0:3126] = 1B ; CTE
            chip5[1155, 0:135] = 1B ; bad CTE
            chip5[1445, 0:315] = 1B ; bad CTE
            chip5[1535:1538, 0:473] = 1B ; bad CTE - messy!
	    chip5[1532:1540,466:474] =1B ; the head of the bad region
            chip5[1618, *] = 1B ; bad CTE
            chip5[1637, *] = 1B ; slightly hot at 3715

; bad pix from superdarks:
            chip5[0:1, *] = 1B
            chip5[*, 4094:4095] = 1B
            chip5[2046:2047, *] = 1B
            chip5[*, 0:1] = 1B
            chip5[124:126, 425] = 1B
            chip5[120:126, 429] = 1B
            chip5[1442, 314:315] = 1B
            chip5[16, 1173] = 1B
            chip5[91, 1103:1104] = 1B
            chip5[1129, 1142:1143] = 1B
            chip5[1961, 1777] = 1B
            chip5[835, 1728] = 1B
            chip5[772, 1767] = 1B
            chip5[478, 1857] = 1B
            chip5[172, 2264] = 1B
            chip5[1780, 2115] = 1B
            chip5[1870, 2433] = 1B
           chip5[1985, 2498] = 1B
           chip5[1211, 2460] = 1B
           chip5[1004, 2891] = 1B
           chip5[1604, 2786:2787] = 1B
           chip5[1603, 2788:2789] = 1B
           chip5[2034, 2801] = 1B
          chip5[1764, 3287] = 1B
          chip5[1646, 3267] = 1B
          chip5[1073, 3385] = 1B
          chip5[873, 3351] = 1B
          chip5[1360, 3724:3725] = 1B
          chip5[1416, 3513] = 1B
          chip5[1430:1435, 3539:3541] = 1B
          chip5[1617:1621, 3715] = 1B
          chip5[1620:1621, 3716:3717] = 1B
          chip5[1637, 3704:3715] = 1B
          chip5[1616:1637, 3715] = 1B
          chip5[375:376, 3837] = 1B
          chip5[460, 3954] = 1B
          chip5[217, 3589] = 1B
          chip5[403, 464] = 1B
          chip5[772, 1766] = 1B
          chip5[126, 2338] = 1B
          chip5[19, 3227] = 1B
          chip5[240, 3744] = 1B
          chip5[118:126, *] = 1B
          chip5[132, 4092] = 1B

          chip5[24, *] = 1B
          chip5[22:23, 0:200] = 1B
          chip5[14:25, 1060:1600] = 1B
          chip5[22:23, 1600:*] = 1B
          chip5[17:19, 3227] = 1B
; bad serial register region          
          chip5[*, 4073:*] = 1B 
          whreg = where(yind gt 4164-0.0434783*xind-10 AND $
                        xind ge 1869 AND xind le 1961)
          chip5[whreg] = 1B
          whreg = where(yind gt 4429-0.178571*xind-10 AND $
                        xind ge 1961 AND xind le 2045)
          chip5[whreg] = 1B
          
          whreg = where(yind gt 4117-0.33333333*xind-10 AND $
                        xind ge 103 AND xind le 115)
          chip5[whreg] = 1B

          whreg = where(yind gt 4256-1.54545*xind-10 AND $
                        xind ge 115 AND xind le 126)
          chip5[whreg] = 1B
   
            whvig = where(xind lt 40)
            chip5[whvig] = chip5[whvig] OR 2B

           ; 50% vignetted corner
            if novigflat then whvig2 = where(yind gt 2995.0006 $
                           +1.4685326*xind $
                           -0.00041152133*xind^2) else begin
                vigflat=mrdfits(vigflatname,i)
                whvig2=where(vigflat lt 0.5)
            endelse

            chip5[whvig2] = chip5[whvig2] OR 2B

            ; define >5% vignetted region
            if novigflat then whpartvig = where(yind gt 2697.6157 $
                              +1.2399439*xind $
                              -0.00016719548*xind^2) else $
              whpartvig=where(vigflat lt 1.)

            chip5[whpartvig] = chip5[whpartvig] OR 4B


            mkhdr, header, chip5
            mosaic_keywords, header, 5
            mwrfits, chip5,  maskfilename, header; chip5
            delvarx, chip5
         end
      6: begin
            chip6 = invvar
            ; this is the cleanest CCD I have ever seen! - DPF


; bad pix from superdarks:
            chip6[0:1, *] = 1B
            chip6[2046:2047, *] = 1B
            chip6[*, 0:1] = 1B
 
            chip6[1189, 455] = 1B
            chip6[797, 771] = 1B
            chip6[407, 809] = 1B
            chip6[1844:1845, 1102] = 1B
            chip6[1433, 1014] = 1B
            chip6[1006, 980] = 1B
            chip6[40:41, 1069] = 1B
            chip6[820, 1532] = 1B
            chip6[1224, 1531] = 1B
            chip6[1705:1706, 1813] = 1B
            chip6[639, 1951] = 1B
            chip6[115, 3353] = 1B
            chip6[1077, 3889:3890] = 1B
            chip6[1286, 3761] = 1B
            chip6[125, 150:2300] = 1B
            mkhdr, header, chip6
            mosaic_keywords, header, 6

            mwrfits, chip6,  maskfilename, header; chip6
            delvarx, chip6
         end
      7: begin
            chip7 = invvar

            chip7[425:427, 0:3415] = 1B ;blocked columns
	    chip7[422:431,3407:3416]=1B
            chip7[555:561, 176:179] = 1B
            chip7[676, 0:2101] = 1B ; CTE
            chip7[1176:1177, 0:3507] = 1B ; bad CTE
	    chip7[1172:1180,3502:3509]=1B
            chip7[1710:1723, 487:497] = 1B



; bad pix from superdarks:
            chip7[0:1, *] = 1B
            chip7[2046:2047, *] = 1B
            chip7[*, 0:1] = 1B
            chip7[776, 67] = 1B
            chip7[1776, 204] = 1B
            chip7[294, 1117] = 1B
            chip7[1927, 1636] = 1B
            chip7[1722, 1994] = 1B
            chip7[83, 2333] = 1B
            chip7[574, 2420] = 1B
            chip7[1011, 2613:2628] = 1B
            chip7[256, 2595] = 1B
            chip7[1156, 3227] = 1B
            chip7[1173, 2678] = 1B
            chip7[155, 3765] = 1B
            chip7[559, 3908:3909] = 1B
            chip7[839, 4040:4045] = 1B
            chip7[960, 1333] = 1B
            chip7[1638, 3339] = 1B
            chip7[1168, 2911] = 1B
            chip7[608, 3786] = 1B
            chip7[608, 3783] = 1B

            chip7[266, 1773] = 1B
            chip7[106, 1556] = 1B
            chip7[707, 1859] = 1B
            chip7[549, 1523] = 1B
            chip7[1721, 1994] = 1B
            chip7[1962, 1748] = 1B
            chip7[1334, 1391] = 1B
            chip7[831, 1215] = 1B
            chip7[764, 1123] = 1B
            chip7[111, 637] = 1B
            chip7[1005, 604] = 1B
            chip7[1788, 180] = 1B
            chip7[1964, 2] = 1B
            mkhdr, header, chip7
            mosaic_keywords, header, 7
            mwrfits, chip7,  maskfilename, header; chip7
            delvarx, chip7
         end
      8: begin
            chip8 = rotate(invvar1, 2)
            chip8[439:440, *] = 1B ; hot col
            chip8[508:512, 0:1021] = 1B ; really whacked
	    chip8[505:516,1014:1023]=1B
            chip8[532, 0:2031] = 1B ; bad col
	    chip8[783,3992:3999]=1B ; JAN
            chip8[805:806, 0:3932] = 1B ; bad CTE
            chip8[804:808,3929:3934]=1B
            chip8[925:934, *] = 1B ; hot cols 
; hot by:      920    0.0625000
;     921     0.187500
;     922     0.312500
;     923     0.375000
;     924     0.625000
;     925     0.875000
;     926      1.37500
;     927      2.00000
;     928      2.12500
;     929      5.62500
;     930      18.2500

            chip8[978:984, 1298:1302] = 1B
            chip8[789:794, 1236:1239] = 1B

; bad pix from superdarks:
            chip8[0:1, *] = 1B
            chip8[2046:2047, *] = 1B
            chip8[*, 0:1] = 1B
            chip8[167, 146] = 1B
            chip8[1763, 128] = 1B
            chip8[1554, 1374] = 1B
            chip8[1091, 1414:1415] = 1B
            chip8[716, 1726] = 1B
            chip8[1683, 1412] = 1B
            chip8[2035, 2000] = 1B
            chip8[1346, 2136:2163] = 1B
            chip8[445, 2398] = 1B
            chip8[1637:1638, 2436] = 1B
            chip8[330, 2733] = 1B
            chip8[456, 2745] = 1B
            chip8[730, 2690] = 1B
            chip8[165, 2994] = 1B
            chip8[844, 3554] = 1B
            chip8[203, 3564] = 1B
            chip8[807, 3928] = 1B
            chip8[ 532, 4040:*] = 1B ; maybe
            chip8[507:513, 1013] = 1B
            chip8[ 1679, 3980] = 1B
            chip8[1998, 3903] = 1B
            chip8[532, 4094] = 1B
            chip8[513, 1012] = 1B
            chip8[508, 1013] = 1B
            chip8[2040, 2214:2215 ] = 1B        
            chip8[1863, 3097] = 1B
            chip8[844, 3354] = 1B
            chip8[1799, 3879] = 1B
            chip8[1758, 3819] = 1B
            chip8[1765, 176] = 1B
            chip8[2027, 84] = 1B
            chip8[1992, 1089] = 1B
            chip8[1897, 1263] = 1B
            chip8[1091, 1409:1424] = 1B
            chip8[615, 1306] = 1B
            chip8[312, 1776] = 1B
            chip8[325, 2045] = 1B
            chip8[1636, 2436] = 1B
            chip8[1135, 2730] = 1B
            chip8[120, 2589] = 1B
            chip8[1520, 3311] = 1B
            chip8[923:930, 3992:3995] = 1B
            chip8[434:440, 2409:2410] = 1B
            chip8[805, *] = 1B
            chip8[2042:2045, 1248:1251] = 1B ; cold spot
            chip8[2021:2022, 1717:1718] = 1B ; cold spot

            whvig = where(xind gt (1969.+yind^2./26.E4))
            chip8 = chip8 OR rotate(invvar1, 2)
            chip8[whvig] = chip8[whvig] OR 2B

          ; 50% vignetted corner
            if novigflat then whvig2 = where(yind gt 4789.6676 $
                           -0.64562068*xind $
                           -0.00010495970*xind^2) else begin
                vigflat=mrdfits(vigflatname,i)
                whvig2=where(vigflat lt 0.5)
            endelse

            chip8[whvig2] = chip8[whvig2] OR 2B

            ; define >5% vignetted region
            if novigflat then whpartvig = where(yind gt 4807.9611 $
                             -1.0267016*xind $
                              +1.7633758e-05*xind^2) else $
              whpartvig=where(vigflat lt 1.)

            chip8[whpartvig] = chip8[whpartvig] OR 4B



            mkhdr, header, chip8
            mosaic_keywords, header, 8
            mwrfits, chip8,  maskfilename, header; chip8
            delvarx, chip8
         end
    endcase
  endfor  
  spawn, 'compress -f '+maskfilename
  return 
end







