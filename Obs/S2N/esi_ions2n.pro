;+ 
; NAME:
;       esi_ions2n
;
; PURPOSE:
;       This procedure produces a plot of either the S/N per pixel,
;       the system efficiency, or the noise from various sources as a
;       function of wavelength with the HIRES spectrograph on Keck I.
;       It is a customized version of a more general routine, intended
;       for invocation from an ION (IDL On the NET) script.
;
; CALLING SEQUENCE:
;       esi_ions2n, plottype, binning, exptime, mag, slitwidth, seeing, $
;               airmass, moondays, minwave, maxwave, deltawave, $
;               XSIZE=xsize,  YSIZE=ysize, SLIT_THRU=slit_thru, RESOLU=resolu, $
;               TABLE=table
;
; INPUTS:
;       plottype: (string) indicates which type of plot to generate.
;       Valid options are:
;               - 'S2N' for signal-to-noise vs. wavelength
;               - 'Eff' for system efficiency vs. wavelength
;               - 'Noise' for noise vs. wavelength
;       Any other value will result in no plot (although certain
;       parameters will be computed and returned.
;
;       binning: (string) indicates CCD binning to use.  Options are:
;               - '1x1'
;               - '2x1'
;               - '3x1'
;               - '2x2'
;
;       exptime: (int) exposure time [sec]
;
;       mag: (double) AB magnitude of the source
;
;       slitwidth: (double) width of the spectrograph slit [arcsec]
;
;       seeing: (double) size of the seeing profile [arcsec]
;
;       airmass: (double) airmass of observation
;
;       moondays: (double) age of moon [days].  New moon corresponds
;       to moondays=0, full to moondays=14.
;
;       minwave: (double) starting wavelength for computation and plot [Angstroms]
;
;       maxwave: (double) ending wavelength for computation and plot [Angstroms]
;
;       deltawave: (double) ending wavelength for computation and plot [Angstroms]
;
; INPUT KEYWORD PARAMETERS:
;       XSIZE: Size of plot in screen x-pixels (default=1000)
;
;       YSIZE: Size of plot in screen y-pixels (default=600)
;
; OUTPUT KEYWORD PARAMETERS:
;       SLIT_THRU: named variable to receive the fraction of light
;       lost due to slit size [percent] 
;
;       RESOLU: named variable to receive the resolution of resulting
;       spectrum [R value] 
;
;       TABLE: named variable to receive a structure containing the
;       wavelength and S/N as arrays
;
; EXAMPLES:
;       1) Generate a plot of S/N for a mag=20.0 object over the
;       wavelength range 4000-8000A.  Also return the slit loss as
;       "slitloss", the resolution as "r", and the wavelength and S/N
;       as a structure called "table":
;
;       plottype = 'S2N'
;       binning = '1x1'
;       exptime = 3600
;       mag = 20.0
;       slitwidth = 1.0
;       seeing = 1.0
;       airmass = 1.0
;       moondays = 0.0
;       minwave = 4000.0
;       maxwave = 8000.0
;       dwave = 100.0
;       esi_ions2n, plottype, binning, exptime, mag, slitwidth, seeing, $
;              airmass, moondays, minwave, maxwave, dwave,
;              SLIT_THRU=slitloss, RESOLU=r, TABLE=table
;
; REVISION HISTORY:
;   04-Jun-2008 JXP     v1.0    Original version
;   07-Jun-2008 GDW     v1.1    Slight rearrangement and addition of header doco
;   24-May-2011 GDW     v1.2    Set mtype=2 for AB mags; 
;                               change esi_calcs2n -> keck_calcs2n;
;                               reformat y-axis label; 
;                               use !x.crange for plotting
;-
;------------------------------------------------------------------------------

pro esi_ions2n, plottype, binning, exptime, mag, slitwidth, seeing, $
                  airmass, minwave, maxwave, deltawave, $
                  XSIZE=xsize,  YSIZE=ysize, SLIT_THRU=slit_thru, RESOLU=resolu, $
                  TABLE=table

  if  N_params() LT 11 then begin 
      print,'Syntax: ', $
            'esi_ions2n, plottype, binning, exptime, mag, slitwidth, seeing,', $
            'airmass, moondays, minwave, maxwave, deltawave,', $
            'XSIZE=xsize,  YSIZE=ysize, SLIT_THRU=slit_thru, RESOLU=resolu,', $
            'TABLE=table'
      return
  endif 

  ;; Initialize
  x_initesi, str_instr, INFIL=infil, STR_TEL=str_tel, SLIT=slit
  str_obs = x_obsinit(infil)

  ;; define STATE
  state = {             $
          nwv: 0L, $
          dwv: 100., $
          wave: fltarr(1000), $
          s2n: fltarr(1000), $
          wvmn: 4000., $
          wvmx: 10000., $
          infil: '', $
          slits: ['0.3', '0.5', '0.75', '1.0'], $
          binning: ['1x1', '2x1', '3x1', '2x2'], $
          deckidx: 0L, $
          pixel: 0., $
          flg_plot: 0, $
          str_instr: str_instr, $ ; PLOTTING LINES
          str_tel: str_tel, $
          str_obs: str_obs, $
          pos: [0.1,0.1,0.95,0.95], $ ; Plotting
          psfile: 0, $
          help: strarr(50), $
          svxymnx: fltarr(4), $
          xymnx: fltarr(4), $
          tmpxy: fltarr(4), $
          size: lonarr(2) $
  }

  ;; Setup
  state.str_instr.bins = long(strmid(strtrim(binning,2),0,1))
  state.str_instr.bind = long(strmid(strtrim(binning,2),2,1))
  state.str_obs.exptime = float(exptime)>1
  state.str_obs.mstar = float(mag) 
  state.str_obs.seeing = float(seeing) > 0.4
  ;state.str_obs.mphase = 0 > long(moondays) < 14
  state.str_obs.airmass = (float(airmass)>1.)
  state.str_instr.swidth = float(slitwidth)
  state.str_obs.mtype = 2 ;; require AB magnitude system

; GDW 2011-Oct-12 change to match HIRES method...
;  resolu = state.str_instr.R / ( state.str_instr.swidth / $
;                                 state.str_tel.plate_scale / $
;                                 (state.str_instr.pixel_size/100) )
  resolu = state.str_instr.R / ( state.str_instr.swidth / $
                                 state.str_instr.scale_para)
  state.dwv = deltawave
  state.wvmn = minwave > 2000.
  state.wvmx = 11000. < maxwave > (state.wvmn+100)

  ;; Wave array
  state.nwv = long((state.wvmx-state.wvmn)/state.dwv) + 1
  state.wave[0:state.nwv-1] = state.wvmn + findgen(state.nwv)*state.dwv

  ;; Calc
  keck_calcs2n, state.wave[0:state.nwv-1], 3, $
    STATE=state, /nopr, S2N=s2n,  FSTRCT=s2n_fstrct
  state.s2n[0:state.nwv-1] = s2n

  ;; Update slit throughput
  slit_thru=(1.-s2n_fstrct.slit0)*100
  
  wave= state.wave[0:state.nwv-1]
  sn= state.s2n[0:state.nwv-1]
  mxx = max(wave, min=mnx)

  table = { $
          wave: wave, $
          s2n: s2n $
          }

  case strtrim(plottype,2) of
      'S2N': state.flg_plot = 0
      'Eff': state.flg_plot = 1
      'Noise': state.flg_plot = 2
      else: return
  endcase

  ;; Set plot window
  clr = getcolor(/load)

  case state.flg_plot of
      0: begin ; S2N

          mxy = max(sn, min=mny)
          state.xymnx = [mnx, mny, mxx, mxy]
          
          binc = state.str_instr.bind
          pixel    = binc*(3e5/state.str_instr.R)
          plot, [0], [0], color=clr.black, background=clr.white, charsize=1.7,$
            xmargin=[10,2], ymargin=[5,2], xtitle='!17Wavelength (Ang)', $
            ytitle='S/N (per '+strtrim(string(pixel,'(f10.2)'),2)+' km/s pix)', $
            yrange=[mny*0.9,mxy*1.1], thick=4, $
            xrange=[mnx*0.9,mxx*1.05], ystyle=1, xstyle=1, psym=1, /nodata

          oplot, [wave], [sn],  psym=-1, color=clr.blue
          
          ;; Label
          xpos = mnx + (mxx-mnx)*0.95
          dy = mxy-mny
      end
      1: begin ; Efficiency
          plot, [0], [0], color=clr.black, background=clr.white, charsize=1.7,$
            xmargin=[10,2], ymargin=[5,3], xtitle='!17Wavelength (Ang)', $
            ytitle='Efficiency Curve', $
            yrange=[0., 0.5], thick=4, $
            xrange=[mnx*0.9,mxx*1.05], ystyle=1, xstyle=1, psym=1, /nodata
          oplot, wave, s2n_fstrct.thru, psym=-1, color=clr.black, thick=3
      end
      2: begin ; Noise
          mxy = max([s2n_fstrct.noise^2, s2n_fstrct.star, s2n_fstrct.sky])
                    
          state.xymnx = [mnx, 0., mxx, mxy]

          plot, [0], [0], color=clr.black, background=clr.white, charsize=1.4,$
            xmargin=[14,4], ymargin=[5,3], xtitle='!17Wavelength (Ang)', $
            ytitle='Variance', $
            yrange=[1.,mxy*1.05], thick=4, $
            xrange=[mnx*0.9,mxx*1.05], ystyle=1, xstyle=1, psym=1, /nodata

          ;; Sky + Obj
          oplot, wave, s2n_fstrct.star, psym=-1, color=clr.black, thick=3
          oplot, wave, s2n_fstrct.sky, psym=-1, color=clr.blue, thick=3
          ;; Readno + dark
          oplot, !x.crange, replicate(s2n_fstrct.noise^2,2), linestyle=2, $
            color=clr.green, thick=2
          oplot, !x.crange, replicate(s2n_fstrct.ndark,2), linestyle=2, $
            color=clr.brown, thick=2

          ;; Label
          xyouts, mnx*0.95, mxy*0.9, 'Object', color=clr.black, chars=2.1
          xyouts, mnx*0.95, mxy*0.8, 'Sky', color=clr.blue, chars=2.1
          xyouts, mnx*0.95, mxy*0.7, 'Dark', color=clr.brown, chars=2.1
          xyouts, mnx*0.95, mxy*0.6, 'ReadNo', color=clr.green, chars=2.1

      end
      else: stop
  endcase
  
  return
end
