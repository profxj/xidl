;; This script provides an example of how to quickly calibrate new HIRES-red
;; (or blue) setups

;; ECHANGL = 0.7058, XDANGL = 1.4180
;; 
;; 0) Plug your setup into the HIRES EFS at:
;;    https://www2.keck.hawaii.edu/inst/hires/HIRES-efs-master/efs.html
;;    to get an idea of what orders are covered. Keep this up as a reference.
;; --------- 
;; 1) Read in hires structure for this reduction
hires = hires_ar('hires-C2-XD_1.41.fits')
setup = 3
;; --------- 
;; 2) Run the automated wavelength solution routine
;;    hires_allarc.pro through with the few_lines keyword. This will
;;    reidentify using archival arcs when it work, but it will not crash
;;    on orders for which auto-id shits the bed.
hires_allarc, hires, setup,/clobber,/few_lines
stop
;; ---------
;; 3) 
;;    Now inspect the QA file gv QA/Arcs03/qa_arcfit_B0015.ps
;;    Decide if you just need to tweak a few orders or whether all the
;;    orders need to be re-fit
;; 
;;    Fitting the Blue Arc [Orders:53-66]
;; 
orders = [65]
;;
;; --------- 
;; 4) Scroll up on the on screen output from hires_allarc.pro
;; for the line which states what the typical shift is between your
;; arc and the template. Record the template filenamename and the
;; shift. Note that shift and typical stretch you get (below) here.
;; for reference.
;;
;; Code will write: 
;; IDL> hires_fitarc:  Using fil /Users/joe/IDL/xidl/Keck/HIRES/CALIBS/ARCS/hires_tmplR2x1B76.idl for our guess.
templatefile='hires_tmplR2x1B76.idl'
;;
;;
;;  x_fitarc: Order   54 Shifting -26 
;;
;; Shift: -26
;; Stretch: 50
;; --------- 
;; 5) Some parameters that could improve the functionality/stability
;;    of x_identify
;;  
toler = 2.0  ;; tolerance between centroid and line
pkwdth= 6.0  ;; set according to spectrqal binning pkwdth = (6L / rbin) > 2L
mxoff=5.0    ;; maximum offset
pksig=7.0    ;; significance
thin=1       ;; use thin lines
;;
;; --------- 
;; 6) Before we identify/re-identify the arc lines launch a x_specplot
;; gui with a calibrated UVES ThAr spectrum for reference (if you are
;; used to id-ing arclines from line plots you will thank me for this
;; step!)
;; 
linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst'
;; You may find it handy to have this linelist open in a text editor
;; emacs linlist &
thar = hires_murphy_thar(wave_thar=wave_thar)
x_arclist, linlist, lines
yrange=[-1.0,1e4]
x_specplot,thar,wav=wave_thar,xsize=1400,ysize=250 $
           ,ytwo=yrange[1]*0.7 + 0*lines.wave $
           ,two_wave=lines.wave,psym2=2
;; 
;; This launches an interactive spectrum of reference ThAr spectrum
;; with the lines in the linelist marked by the purple
;; stars. After hand-IDing a few lines in x_reidentify, mouse over the
;; lines in this reference spectrum to get their wavelength and use the 'M'
;; feature of x_reidentify to quickly add more lines. 
;;
;; --------- 
;; 7) Tweak or re-fit the arc using hires_tweakarc and x_identify 
;; 
hires_tweakarc, 'Arcs/Fits/Arc_B0015_fit.idl', orders $
                ,templatefile $
                ,xsize=1400,ysize=400,pkwdth=pkwdth,toler=toler $
                ,mxoff=mxoff,pksig=pksig,thin=thin
;;
;; --------- 
;; 8) Option A):
;;    Once you are happy with your fits, you can update the 1-d
;;    wavelength solutions and the 2-d wavelength fits and other steps
;;    that build on them by running:
;; 
;;    hires_allarc_sngl('Raw/hires0015.fits', hires, setup, 1, /chk,
;;       /nooned, /clobber)
;;
;;    Note that you need to run the /nooned option, otherwise you will
;;    overwrite and *LOSE* your previous fit.
;;
;;    This option makes the most sense if you just made a small tweak
;;    to a few orders. I'm not a fan of this option, both because of the
;;    possibility of overwriting the fit, and beacuse the 1-d QA files
;;    i.e. QA/Arcs03/qa_arcfit_B0029.ps do not get updated (although the
;;    2-d QA files do get updated). 
;;
;; --------- 
;; 8) Option B)
;;    Archive solution as a new hires template file:
;;
fitfil= 'Arcs/Fits/Arc_B0015_fit.idl'
orders= [53,66]
outfil= getenv('XIDL_DIR') + $
        '/Keck/HIRES/CALIBS/ARCS/hires_tmplR3x1B3.idl'
hires_mktempl,fitfil,orders,outfil
;; 
;;     You also need to add a line like below for your template file
;;
;;     Name           chip side    ECH    XD      Dck  specbin  spatbin  min_o max_o
;; hires_tmplR3x1B2.idl  1    RED  0.  1.3135   C2    1     3    53  67
;; 
;;     To the hires_templ file at:
;; 
;;     getenv('IDL_DIR') + '/Keck/HIRES/Redux/pro/Arcs/Templates'
;;
;;     If you have svn commit priveleges commit this arc and the
;;     template file (contact JXP or JFH about that). 
;;

