;+ 
; NAME:
; x_obss2n
;    Version 1.1
;
; PURPOSE:
;     This program computes count rates and expected S/N for
;     the HIRES spectrograph.  It mimics the 'old' FORTRAN code.
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
;     HISTORY
;     ??/??/?? Written by G. Donald Penrod 
;     11/13/89 Modified by S. Vogt       - made inputs less confusing
;     06/08/92 Modified by M. Keane      - changed from Hamilton to HIRES 
;     02/07/96 Modified by C. Churchill  - structured queries by function
;                                        - set defaults for Decker C1
;                                        - added comments
;     20-Oct-2005 Ported to IDL by JXP
;-
;------------------------------------------------------------------------------
pro x_input_obss2n, infil, str_instr, str_tel

 ;; Parse infil
 readcol, infil, card, val, FORMAT='A,A'

 ;; Instrument
 mtch = where(strtrim(card,2) EQ 'INSTR',nmt)

 ;; HIRES cards
 mtch = where(strtrim(card,2) EQ 'HIRES_DECKER',nmt)
 if nmt NE 0 then hires_decker = val[mtch[0]]

 if nmt EQ 0 then begin
     ;; OLD HIRES
     x_inithires, str_instr, KECKTEL=str_tel, DECKER=hires_decker
 endif

 return
end
 

pro x_obss2n, infil, GUI=gui
;
  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_obss2n, infil, GUI=, [v1.0]'
      return
  endif 

  ;; 
  if not keyword_set(NWV) then nwv = 5L

  ;; Initialize
  x_init_obss2n, infil, str_instr, str_tel

  ;; Wave val
  wv_val = str_instr.wvmnx[0] + $
    (str_instr.wvmnx[1]-str_instr.wvmnx[0])*findgen(nwv)/float(nwv)
c
c  compute order, communicate blaze center, query
       m      = nint(MLAMBDA/wave)
       center = real_get('Blaze center (A)', MLAMBDA/m, 3000., 10000.)
c
c  communicate free spectral range, query
       fsr = real_get('Free spectral range (A)', center/m, 1., 1000.)
c
c  compute the order separation
       if (center .lt. 3800.) then 
        sep    = 2.*DELY*(MLAMBDA/m - MLAMBDA/(m+1)) 
        iorder = 2
       else 
        sep    = DELY*(MLAMBDA/m - MLAMBDA/(m+1)) 
        iorder = 1
       end if 
c
c  query for cross disperser
        iorder = integer_get('Enter cross disperser order',
     $       iorder, 1, 2)
c
      end if
c
c
c
c
c
      if ((iflag.eq.0).or.(iflag.eq.2)) then
c
       write(6,*) ' --- ENTER SPECTROGRAPH SETTINGS ---'
c
c  slit height query 
       height = real_get('Enter slit height (arcsecs)', height, 
     $                    0., sep)
c
c  slit width query 
       width = real_get('Enter slit width (arcsecs)', 
     $                  width, 0.1, 10.)
c
c  dark current query
       dark = real_get(
     $        'Enter dark count (electrons/unbinned pixel/hour)', 
     $        dark, 0., 1000.)
c
c  read noise query
       read = real_get('Enter readout noise (electrons/pixel)',
     $        read,  0., 100.)
c
c  dispersion binning query
       binc = integer_get('Enter binning factor (dispersion)',
     $        binc, 1, 2048)
c
c  cross-dispersion query
       binr = integer_get('Enter binning factor (cross dispersion)',
     $        max(1, rows/3), 1, 2048)
c
      end if
c
c
c
c
c
c
      if ((iflag.eq.0).or.(iflag.eq.3)) then
c
       write(6,*) ' --- ENTER OBSERVING CONDITIONS ---'
c
c  seeing query
       seeing = real_get('Enter seeing FWHM [arcsecs]', 
     $                    seeing, 0.25, 5.)
c
c  lunar phase query
       phase = real_get('Lunar phase (days before/after new)', 
     $                   phase, -15., 15.)
       phase = abs(phase)
c
c  air mass query
       air = real_get('Enter airmass', air, 1.0, 10.0)
c
c  integration time query
       time = real_get('Enter exposure time (seconds)', 
     $                  time, 0., 24*3600.)
c
      end if 
c
c
c
c
c
c
c
      if ((iflag.eq.0).or.(iflag.eq.4)) then
c
       write(6,*) ' --- ENTER OBJECT BRIGHTNESS ---'
c  get the magnitude of the Object, get the mag system
c
       mstar = real_get('Enter magnitude of Object', mstar, -5., 25.)
       mtype = integer_get('What type of magnitude? (Johnson=1, AB=2)',
     $                     mtype, 1, 2)
      end if
c
c
c
c
c
c
c  it is time... compute the star counts and noise
      write(6,*) ' ------------------------------------------------'
c
      nj0     = flux(wave)
      nab0    = ab(wave)
      if (mtype .eq. 1) then 
       n0 = nj0 
      else 
       n0 = nab0 
      end if 
      slit0   = slit(width/seeing, height/seeing, 0., 0.)
      rows    = int(3*seeing/SCALE_PERP+0.999)
      columns = max(2., min(width,3.*seeing)/SCALE_PARA)
      nsky    = float((int(height/SCALE_PERP+0.999) - rows))/rows     
c
c  communicate slit width projection (resolution elements)
      write(*,1000) columns
 1000 format(' Slit width projects to ', f3.1, ' pixels')
c
c  communicate rows of data
      write(*,1010) rows
 1010 format(' Object rows     ', i3)
c
c  communicate total cross dispersion rows
      write(*,1020) max(0, nint(nsky*rows))
 1020 format(' Object+Sky rows ', i3)
c       
c  warn if number of sky rows is zero
      if (nsky .le. 0.) then 
       nsky = 1.e8
       write(*,*) ' *** Slit is too short for sky subtraction.' 
      end if
c
c  grab through put from within subroutine and communicate
      thru = thruput(center,iorder)
      thru = blaze(wave, center, fsr)*thru
      write(*, 1001) nint(wave), 100.*thru
 1001 format(' Single order efficiency at ', i5, 'A is ', f5.1, '%')
c
      pixel    = binc*(wave/R)
      slarea   = 3*seeing*width 
      extinct  = trans(wave)
      slit1    = slit0*10**(-0.4*extinct*air)
c
      star     = n0*(10**(-0.4*mstar))*area*thru*slit1*pixel*time
      dstar    = sqrt(star)
      noise    = read/sqrt(float(binr*binc))
      magsky   = msky(wave, phase)
      projslit = columns*(wave/R)
c
c  compute the sky counts and noise
      sky   = nj0*(10**(-0.4*magsky))*slarea*area*thru*pixel*time
      dsky  = sqrt(sky)
c
c  compute the full signal to noise
      ndark  = binc*dark*rows*time/3600.
      ddark  = sqrt(ndark)
      tnoise = sqrt(star+(1.+1./nsky)*(noise*noise+sky+ndark))
      sn     = star/tnoise
c
      write (6,1160) star,dstar,sky,dsky,ndark,ddark,noise,
     1       star,tnoise,sn, pixel, sn*sqrt(columns), projslit
c
 1160 format (/,11x,'Object counts =',f8.0,1x,f7.1,
     1        /,11x,'Sky counts    =',f8.0,3x,f5.1,
     2        /,11x,'Dark counts   =',f8.0,3x,f5.1,
     3        /,11x,'Readout       =',11x,f5.1,
     4        /,11x,'Net Object    =',f8.0,1x,f7.1,
     5       //,1x,'Net S/N: ',f6.0,' per ', 3pf4.0,' mA pixel',
     6        /,1x,'         ',0pf6.0,' per ', 3pf4.0,
     $              ' mA resolution element')
c
c
c  query for next action
      write(6,*) ' ------------------------------------------------'
        iflag = integer_get('Change? 0=all 1=wave 2=spec 3=obs 4=mag',
     $                       iflag, 0, 5)
      if (iflag.eq.5) stop ' normal termination'
      go to 130
c
      end
c
c..............................................................................
c..............................................................................
c
c
	function flux(wave)
c
c
c  no I/O
c..............................................................................
c
	dimension xwave(16),xflux(16)
	data xwave /3062,3312,3562,3812,4062,4212,4462,4712,4962,5425,5925,
     1       6425,6925,7750,8350,11050/
	data xflux /569,586,590,1326,1770,1707,1530,1356,1257,1054,886,749,
     1       641,502,435,263/
c
c
	do 100 i=2,20
           if (wave.le.xwave(i)) go to 200
 100    continue
c
 200    xw1=xwave(i-1)
	xw2=xwave(i)
	xf1=xflux(i-1)
	xf2=xflux(i)
	frac=(wave-xw1)/(xw2-xw1)
	flux=xf1+(xf2-xf1)*frac
c
	return
	end
c
c..............................................................................
c..............................................................................
c
c
	function thruput(wave,iorder)
c
c
c  HIRES throughput estimates revised by S.Vogt, 23 July 93
c  I/O - queries for cross disperser number (1 or 2) 
c      - communicates spectrograph throughput efficiency
c..............................................................................
c
        parameter (NWAVE = 11) 
        parameter (NCD = 1) 
        parameter (NORDER = 2) 
c
        parameter (FUDGE=1.00)
c
	dimension xwave(NWAVE),xthru(NCD, NORDER, NWAVE)
        data xwave /3000, 3200, 3500, 3800, 4000, 
     $              4500, 5000, 6000, 7000, 8000, 
     $              9500/
c
c                        CD
c                         1
c
c                       ORDER 
c                     1      2
c
	data xthru /1.e-4, 0.003,     ! 3000 A
     $              2.e-4, 0.008,     ! 3200 A
     $              0.008, 0.020,     ! 3500 A
     $              0.020, 0.019,     ! 3800OA A
     $              0.035, 0.018,     ! 4000 A
     $              0.061, 0.009,     ! 4500 A
     $              0.077, 0.003,     ! 5000 A
     $              0.080, 3.e-4,     ! 6000 A
     $              0.068, 3.e-4,     ! 7000 A
     $              0.051, 3.e-4,     ! 8000 A
     $              0.017, 3.e-4 /    ! 9500 A
c
c
c  choose default cross disperser  
        icd = 1 
c
c
	do 100 i=2,11
	  if (wave.le.xwave(i)) go to 200
 100      continue
c
c  linearly interpolate the efficiency
 200    xw1=xwave(i-1)
	xw2=xwave(i)
	xt1=xthru(icd,iorder,i)
	xt2=xthru(icd,iorder,i-1)
	frac=(wave-xw1)/(xw2-xw1)
	thruput=(xt2-(xt2-xt1)*frac)/FUDGE
c
c  and communicate
        write (*,1010) nint(wave), 100*thruput
 1010   format(' Spectrograph efficiency at ', 
     $       i5, 'A (blaze peak) is ', f4.1, '%')
c
	return
	end
c
c..............................................................................
c..............................................................................
c
c
	function trans(wave)
c
c
c  Extintion at Mauna Kea, taken from CFHT Bulletin, 19, 16 (1988). 
c  I/O - communicate extinction at this wavelenth
c..............................................................................
c
	dimension xwave(25),xthru(25)
	data xwave /3000,3100,3200, 3300, 3400,
     $              3500,3600,3700, 3800, 3900,
     $              4000,4250,4500, 4750, 5000,
     $              5250,5500,5750, 6000, 6500,
     $              7000,8000,9000,10000,12000/
c
	data xthru /4.90, 1.37, 0.82, 0.57, 0.51, 
     $              0.42, 0.37, 0.33, 0.30, 0.27, 
     $              0.25, 0.21, 0.17, 0.14, 0.13, 
     $              0.12, 0.12, 0.12, 0.11, 0.11,
     $              0.10, 0.07, 0.05, 0.04, 0.03/
c
	do 100 i=2,25
	  if (wave.le.xwave(i)) go to 200
  100     continue
c
  200   xw1=xwave(i-1)
	xw2=xwave(i)
	xt1=xthru(i)
	xt2=xthru(i-1)
	frac=(wave-xw1)/(xw2-xw1)
	trans=xt2-(xt2-xt1)*frac
c
        write (*,1000) int(wave), trans
 1000   format(' Extinction at ', i5, 
     $       'A is ', f4.2, ' magnitudes/airmass')
c
	return
	end
c
c..............................................................................
c..............................................................................
c
c
	real function msky(wave, phase)
c
c
c     Sky brightness at Mauna Kea. 
c..............................................................................
c
        parameter (NWAVE=6, NPHASE=5)
        dimension xwave(NWAVE), xphase(NPHASE)
	dimension xsky(NWAVE,NPHASE)
c
	data xwave /3500, 4200, 5500, 6700, 7800, 22000/
        data xphase /0., 3., 7., 10., 14./     
	data xsky / 22.4, 23.0, 21.9, 21.2, 19.9, 12.0,
     $              21.5, 22.4, 21.7, 20.8, 19.9, 12.0, 
     $              19.9, 21.6, 21.4, 20.6, 19.7, 12.0,
     $              18.5, 20.7, 20.7, 20.3, 19.5, 12.0,
     $              17.0, 19.5, 20.0, 19.9, 19.2, 12.0/
c
	do 100 j=NWAVE-1, 1, -1
           if (wave.ge.xwave(j)) go to 200
 100    continue
        j = j+1 
c
 200    continue 
        do 300 k=NPHASE-1, 1, -1
           if(phase.ge.xphase(k)) go to 400
 300    continue
        k = k+1
c
 400    t = (wave - xwave(j))/(xwave(j+1)-xwave(j))
        u = (phase - xphase(k))/(xphase(k+1)-xphase(k))
c       
        msky = (1.-t)*(1.-u)*xsky(j,k) + t*(1-u)*xsky(j+1,k)
     $        + t*u*xsky(j+1,k+1) +(1-t)*u*xsky(j,k+1)
c        
        write (*, 1000) nint(wave), msky
 1000   format(1x,'Sky brightness at ',  i5, 'A is ', f4.1,
     $       ' mag/arcsec**2')
	return
	end
c
c..............................................................................
c..............................................................................
c
c
      real function blaze(lambda, center, fsr) 
c
c
c     I/O - communicate % of blaze function peak at lambda
c
c..............................................................................
c
      real lambda
      real center, fsr
c
      parameter (PI = 3.14159268)

      gamma = PI * (center - lambda)/fsr
c
      if (gamma .eq. 0.) Then 
         blaze = 1.00
      else
         blaze = (sin(gamma)/gamma)**2
      end if 
c
      write(*,1000) 100*blaze
 1000 format(' Blaze function is ', f4.1, '%')
c
      return 
      end 
c
c..............................................................................
c..............................................................................
c
c
      function ab(wave)
c
c
c..............................................................................
      real wave
c     
      real FLAM0, h, c
      parameter (FLAM0=3.54E-9)
      parameter (h = 6.626E-27)
      parameter (c = 3.00E10)
c     
      ab = FLAM0*wave*1E-8/(c*h)
c
      return
      end
c
c..............................................................................
