;+ 
; NAME:
; hires_arcformm
;     Version 1.1
;
; PURPOSE:
;  Processes a bunch of Arc fits into a file format for Michael Murphy
;  to analyze.
;
; CALLING SEQUENCE:
;   
;  hires_arcformm, 
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
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_fitarc
;
; REVISION HISTORY:
;   01-Jul-2007 Created by JXP 
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; .com x_fitarc
; hires_arcformm, 1, '/raid/Keck/HIRES/ThAr/HIRESb/'
; hires_arcformm, 2, '/raid/Keck/HIRES/ThAr/HIRESb_2/'
; hires_arcformm, 3, '/raid/Keck/HIRES/ThAr/HIRESb_3/', nchip=1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_arcformm, flg, outdir, NCHIP=nchip


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_arcformm, flg, outdir [v1.0]'
      return
  endif 
;  cd, getenv('XIDL_DIR')+'/Spec/Arcs/', CURRENT=current
;  RESOLVE_ROUTINE, 'x_fitarc'
;  cd, current

  if not keyword_set(NCHIP) then nchip = 3
  hires_rslvall

  tmp = { $
        direct: '', $
        setup: 0L, $
        xdangl: 0., $
        echangl: 0., $
        decker: '' $
        }
        
  case flg of 
      1: begin  ;; HIRESb (Aug 2006)
          pth = '/Volumes/data/Keck/HIRES/RawData/'
          dirs = ['2006-Aug-17', $
                  '2006-Aug-18', $
                  '2006-Aug-19', $
                  '2006-Aug-20']
          gdarcs = [ [0,1], $  ; 17Aug, setup 1
                     [0,2], $
                     [0,3], $
                     [1,1], $
                     [1,2], $
                     [1,3], $
                     [2,2], $
                     [2,3], $
                     [3,2] $
                   ]
      end
      2: begin  ;; HIRESb (Aug 2006) -- Select two setups
          pth = '/Volumes/data/Keck/HIRES/RawData/'
          dirs = ['2006-Aug-17', $
                  '2006-Aug-18']
          gdarcs = [ [0,3], $  ; 17Aug, setup 3
                     [1,1]  $  ; 18Aug, setup 1
                   ]
      end
      3: begin  ;; HIRESb (Aug 2006) -- Select two setups
          pth = '/Volumes/data/Keck/HIRES/RawData/'
          dirs = ['2006-Aug-18']
          gdarcs = [ [0,1]  $  ; 18Aug, setup 1
                   ]
      end
      4: begin  ;; HIRESb J2123 (Aug 2006) -- E3 only
          pth = '/Volumes/data/Keck/HIRES/RawData/'
          dirs = ['2006-Aug-19', $
                  '2006-Aug-20']
          gdarcs = [ [0,4], $  ; 19Aug, setup 4
                     [1,6]  $  ; 20Aug, setup 6
                   ]
      end
      else: stop
  endcase

  sz = size(gdarcs,/dimen)
  if n_elements(sz) GT 1 then nset = sz[1] else nset = 1
  strct= replicate(tmp, nset)
  for qq=0L,nset-1 do begin
      ;; Fill up
      strct[qq].direct = dirs[gdarcs[0,qq]]
      strct[qq].setup = gdarcs[1,qq]
      ;; Change directory
      cd, pth+strct[qq].direct

      ;; Open hires file
      hires = hires_ar()

      ;; Parse on setup
      gd = where(hires.setup EQ strct[qq].setup)
      strct[qq].xdangl = hires[gd[0]].xdangl
      strct[qq].echangl = hires[gd[0]].echangl
      strct[qq].decker = strtrim(hires[gd[0]].decker,2)

      ;;Create the directory
      newdir = '/'+strct[qq].direct+'_'+ $
               string(strct[qq].xdangl,format='(f6.4)')+$
               '_'+strtrim(string(strct[qq].echangl,format='(f7.4)'),2)+'_'+$
               strct[qq].decker

      a = file_search(outdir+newdir, count=na)
      if na EQ 0 then spawn, '\mkdir '+outdir+newdir

      ;; Find the Arcs
      arcs = where(hires.setup EQ strct[qq].setup and $
                   hires.type EQ 'ARC' and hires.flg_anly EQ 1, narc)
      if narc EQ 0 then stop

      ;; Loop
      for ss=0L,narc-1 do begin
          for ii=1L,nchip do begin
              fil = hires_getfil('arc_fit', strct[qq].setup, CHIP=ii, $
                                 FRAME=hires[arcs[ss]].frame, /NAME)
              ;; Process
              i1 = strpos(fil,'/',/reverse_search) > 0
              i2 = strpos(fil,'.idl')
              fitsfile = strmid(fil,i1,i2-i1+1)+'fits'

              x_fitarc_fitsout, fil, FITSFILE=outdir+newdir+fitsfile, $
                                MJD=hires[arcs[ss]].date, $
                                EXPTIME=hires[arcs[ss]].exp, GDO=gdordr

              ;; Read in 2D data
              out_str = hires_getfil('arc_2Dfit', strct[qq].setup, CHIP=ii, $
                                 FRAME=hires[arcs[ss]].frame)
              mwrfits, out_str, outdir+newdir+fitsfile

              ;; Order structure
              ordr_str = hires_getfil('ordr_str', strct[qq].setup, CHIP=ii, $
                                      FRAME=hires[arcs[ss]].frame)
              mwrfits, ordr_str[gdordr], outdir+newdir+fitsfile
              ;; Compress
              spawn, 'gzip -f '+outdir+newdir+fitsfile
          endfor
      endfor
  endfor
         
          
  
; All done
  print, 'hires_arcformm: All done!'

  return
end
