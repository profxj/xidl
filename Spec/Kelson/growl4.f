
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine growl4( ndata, raw0, disco100, telescope, graname,
     &        nbig, linelist,arc,nimages,zero,wtol,fwhm1,red,lamps,
     &  dichroic,redside,err,
     &        slice, debug, nsname, lamb1, lamb2, d0, wfft0,bth0,scat)

        real xco(4096), xc(2048)
        real raw0(ndata)
        real raw(4096)
        real draw(4096), y2(4096)
        real xp(64), getmax, cprev, subsky(128)
        real disco100(100), disco(50), chisq, lamb1, lamb2
        logical slice, skk, wfftonly
        integer wfft,wfft0
        integer lobt, hibt, zero, err, nsname
        integer minlines
        character*(*) telescope

        real bincen(64), smdev(2048)
        logical skip(64), arc, getmore, red, debug(10), reject(2048)
        integer rq(2048)

        integer nbs(10)
        real bth0
        logical isaqr(2048)

        real ldb(2048), fwhm16
        real lspec(4096), llam(4096)
        character*40 lname(0:2048)
        real lstrong(0:2048), d0, midp(64), com(64), comr(64)
        real gimmes(64)

        real sraw(4096), srcor(4096), cref(32)
        real qsraw(4096)
        complex sfraw(2049), sfcor(2049)
        complex fraw(2049), fcor(2049)
        complex fref(2049)

        integer ia(2), ib(2), lista(50)
        character*80 line, lampline
        character*(*) linelist, graname, lamps

        real covar(50,50)

        real*8 u(2048,24)
        real*8 v(24,24)
        real*8 w(24)
        real*8 shifts(2)
        real hold(2048)
        integer index(4096)
        integer spord(2048)
        integer pind(2048)
        integer rr(2)

        logical hg, ne, ar, he, already, dichroic, redside, zn, cd

        real sraw2(4096), lspec2(4096)

        integer axlen(7), naxis, dtype, ierr, im(10)
        integer row, nslits, ref
        integer imax, dmax(2048), half, attach(2048)
        integer slit1(2), slit2(2), rslit
        integer ix1m(64),ix2m(64), refmax(2048), ihold(2048)
        real s2n(64)
        real wmax(2048), fwhm(2048), fwhm2(2048), dev(2048)
        real tmax(2048), som(64)
        real cen(2048), wt(2048), ymax(2048), temp(2048), ht(2048)
        real fx(2048),fy(2048),fw(2048), pt(24)
        real h(2048)

        real dumx(2048), dumy(2048), dumw(2048), duml(50)
        integer fi(2048)

        logical dum, used(2048), kill(2048), triplet(2048), goodgrab

        real tcen(2048,2)
        integer sbins(2),dbins(2)

        real temp1(2048)
        real temp2(2048)
        real temp3(2048)
        real tri1(2048)
        real tri2(2048)
        real tri3(512,3,3,3)
        real tri4(2048,3)
        real nbt

        real lambda

        external fleg, fpoly, fuser, fleg1

        wfft = wfft0
        ncvm = 50

        do j = 1,50
         disco(j) = disco100(j)
        end do

        do i=1,4096
        xco(i) = float(i)
         if (i.le.ndata) then
           raw(i) = raw0(i)
           write(29,*) i,raw0(i)
          else
           raw(i) = 0.0
         end if
        end do
        do i=1,2048
        used(i) = .false.
        triplet(i) = .false.
        end do

        goodgrab=.false.

        numlfit = 0
        fbin = 4.0
        if ( d0.eq.0.0) then
        if ( graname.eq."cosmic") disp0 = 3.25
        else
        disp0 = d0
        end if

        if ( arc.and.nimages.eq.1) then
        if (telescope.eq."KECK") then
          hg=.false.
          zn=.false.
          cd=.false.
          ne=.false.
          ar=.false.
          if ( lamps(1:1).eq.'1') hg=.true.
          if ( lamps(3:3).eq.'1') ne=.true.
          if ( lamps(5:5).eq.'1') ar=.true.
          if ( lamps(7:7).eq.'1') zn=.true.
          if ( lamps(9:9).eq.'1') cd=.true.
        else if (telescope.eq."COSMIC") then
          hg=.true.
          ne=.false.
          ar=.true.
          he=.false.
          zn=.false.
          cd=.false.
        else
          hg=.true.
          ne=.true.
          ar=.true.
          he=.true.
          zn=.true.
          cd=.true.
        end if
        write(6,1001) hg,ne,ar,zn,cd
 1001        format('Lamps:   Hg - ',l1,2x,'Ne - ',l1,2x,'Ar - ',l1,2x,
     &         'Zn - ',l1,2x,'Cd - ',l1)
        end if

        write(6,*)
        write(6,30) telescope, graname,disp0
  30        format(3x,'Telescope= ',a8,3x,' Grating= ',a8,
     &           ' Approx. Dispersion= ', f4.2,' A/pix  ',$)
        if ( red) then
         write(6,31)
        else
         write(6,*)
        end if
  31        format('in the red')

        write(6,*)

        yp1 = raw(2) - raw(1)
        ypn = raw(ndata) - raw(ndata-1)
        call spline(xco, raw, ndata, yp1, ypn, draw)

        do i=1,80
         if (linelist(i:i).eq.' ') goto 444
        end do  
  444   continue
        linelist = linelist(1:i-1)

c        write(6,*) 'opening ',linelist

cccccccc

        ip = 0
        open(unit=2,file=linelist,status='old')
        i = 0
        do ii=1,9999
        read(2,'(a80)',end=600) line
        i = i + 1
        do j=1,80
        if (line(j:j).ne.' ') goto 598
        end do
  598        jw1 = j
        do j=jw1,80
        if (line(j:j).eq.' ') goto 597
        end do
  597        jw2 = j - 1
        read(line(jw1:jw2),*) ldb(i)
        do j=jw2+1,80
        if (line(j:j).ne.' ') goto 596
        end do
  596        jw = j
        lname(i)='unknown'
        read(line(jw:jw+39),'(a40)') lname(i)
        used(i) = .false.
        kill(i) = .false.
        if ( arc.and.nimages.eq.1) then
         lstrong(i) = 0.0
         if ( lname(i)(1:2).eq.'Hg') lstrong(i) = 300.0
         if ( lname(i)(1:2).eq.'Zn') lstrong(i) = 1000.0
         if ( lname(i)(1:2).eq.'Cd') lstrong(i) = 1000.0
         if ( lname(i)(1:2).eq.'Kr') lstrong(i) = 50.0
         if ( lname(i)(1:2).eq.'Ne') lstrong(i) = 1000.0
         if ( lname(i)(1:2).eq.'He') lstrong(i) = 1000.0
         if ( lname(i)(1:2).eq.'Ar') lstrong(i) = 1000.0
         if (dichroic.and..not.redside) then
          if ( lname(i)(1:2).eq.'Ne') lstrong(i) = 0.0
          if ( lname(i)(1:2).eq.'He') lstrong(i) = 100.0
          if ( lname(i)(1:2).eq.'Ar') lstrong(i) = 100.0
         end if
         if ( lname(i)(1:2).eq.'AI') lstrong(i) = 100.0
         if ( lname(i)(1:3).eq.'AII') lstrong(i) = 100.0
         if ( lname(i)(1:2).eq.'Kr'.and.ldb(i).gt.7e3) lstrong(i)=1e3
         if ( lname(i)(1:2).eq.'Ar'.and.ldb(i).gt.7e3) lstrong(i)=1e3
         if ( lname(i)(1:2).eq.'Ar'.and.ldb(i).gt.9.2e3) lstrong(i)=3e2
         if ( int(ldb(i)).eq.3649) lstrong(i) = 100.0
         if ( int(ldb(i)).eq.3663) lstrong(i) = 10.0
         if ( int(ldb(i)).eq.4046) lstrong(i) = 100.0
         if ( int(ldb(i)).eq.4077) lstrong(i) = 10.0
         if ( int(ldb(i)).eq.4358) lstrong(i) = 300.0
         if ( int(ldb(i)).eq.4916) lstrong(i) = 500.0
         if ( int(ldb(i)).eq.5460) lstrong(i) = 500.0
         if ( int(ldb(i)).eq.6123) lstrong(i) = 5.0
        if (.not.hg.and.lname(i)(1:2).eq.'Hg') kill(i) = .true.
        if (.not.hg.and.lname(i)(1:2).eq.'Kr') kill(i) = .true.
        if (.not.cd.and.lname(i)(1:2).eq.'Cd') kill(i) = .true.
        if (.not.zn.and.lname(i)(1:2).eq.'Zn') kill(i) = .true.
        if (.not.ne.and.lname(i)(1:2).eq.'Ne') kill(i) = .true.
        if (.not.he.and.lname(i)(1:2).eq.'He') kill(i) = .true.
        if (.not.ar.and.lname(i)(1:1).eq.'A') kill(i) = .true.
c         if (lname(i)(1:2).eq.'He') kill(i)=.true.
         if (.not.hg.and.lname(i)(1:2).eq.'Hg') lstrong(i) = -99
         if (.not.hg.and.lname(i)(1:2).eq.'Kr') lstrong(i) = -99
         if (.not.ne.and.lname(i)(1:2).eq.'Ne') lstrong(i) = -99
         if (.not.ar.and.lname(i)(1:1).eq.'A') lstrong(i) = -99
         if (.not.cd.and.lname(i)(1:2).eq.'Cd') lstrong(i) = -99
         if (.not.zn.and.lname(i)(1:2).eq.'Zn') lstrong(i) = -99
c         if (lname(i)(1:2).eq.'He') lstrong(i) = -99
        else
        lstrong(i) = -3.0
        if ( lname(i)(1:5).eq.'[O I]') then
         lstrong(i) = 2.0
         if (nint(ldb(i)).eq.5577) lstrong(i) = 10000.0
         if (nint(ldb(i)).eq.6300) lstrong(i) = 5000.0
         if (nint(ldb(i)).eq.6364) lstrong(i) = 2000.0
        else if ( lname(i)(1:2).eq.'Hg') then
         if (int(ldb(i)).eq.3649) lstrong(i) = 100.0
         if (int(ldb(i)).eq.4046) lstrong(i) = 500.0
         if (int(ldb(i)).eq.4077) lstrong(i) = 50.0
         if (int(ldb(i)).eq.4358) lstrong(i) = 500.0
         if (int(ldb(i)).eq.4916) lstrong(i) = 500.0
         if (int(ldb(i)).eq.5460) lstrong(i) = 500.0
        else if ( lname(i)(1:4).eq.'Na I') then
          if (telescope.eq."COSMIC") then
                lstrong(i) = 5000.0
                if (disp0.lt.3.0) lstrong(i) = 4000.0
          else if (telescope.eq."KECK") then
                lstrong(i) = 500.0
          end if
          if (disp0.gt.1.5) ldb(i) = 5892.935
        else
        if ( lname(i)(5:5).eq.'P') then
         ip1 = 100
         if ( lname(i)(7:7).eq.'(') then
          read(lname(i)(8:9),*) realip1
          ip1 = int(realip1)
         end if
         if ( lname(i)(8:8).eq.'(') then
          read(lname(i)(9:10),*) realip1
          ip1 = int(realip1)
         end if
         if ( lname(i)(6:6).eq.'1') then
           lstrong(i) = 20*(ldb(i)/2500)**3 * exp(-(ip1-3.5)**2)
         end if
         if ( lname(i)(6:6).eq.'2') then
           lstrong(i) = 10*(ldb(i)/2500)**3 * exp(-(ip1-3.5)**2)
         end if
         if (ip1.gt.8) lstrong(i) = 0
         if ( (graname.ne.'300gmm'.and.ip1.lt.8).or.
     &     (graname.eq.'300gmm'.and. lname(i)(6:6).eq.'1'.and.ip1.lt.6))
     &       then
          ip = ip + 1
          pind(ip) = i
         end if
        end if
        isaqr(i) = .false.
        if ( lname(i)(5:6).eq.'Q1') then
           lstrong(i) = 40*(ldb(i)/5500)**3
           isaqr(i) = .true.
        end if
        if ( lname(i)(5:6).eq.'R1')  then
           lstrong(i) = 20*(ldb(i)/5500)**3
           isaqr(i) = .true.
        end if
        if ( lname(i)(5:7).eq.'Q1(')  then
           lstrong(i) = 40*(ldb(i)/5500)**3
           isaqr(i) = .true.
        end if
        if ( lname(i)(5:7).eq.'R1(')  then
           lstrong(i) = 20*(ldb(i)/5500)**3
           isaqr(i) = .true.
        end if
        if ( lname(i)(5:7).eq.'Q2(')  then
           lstrong(i) = 40*(ldb(i)/5500)**3
           isaqr(i) = .true.
        end if
        if ( lname(i)(5:7).eq.'R2(')  then
           lstrong(i) = 20*(ldb(i)/5500)**3
           isaqr(i) = .true.
        end if
c        if (isaqr(i).and.d0.gt.4.5) then
c              isaqr(i) = .false.
c              lstrong(i) = lstrong(i) * 10
c        end if
        if ( lname(i)(1:1).eq.'O')  then
           lstrong(i) = 80.0
           isaqr(i) = .true.
        end if
        if ( ldb(i).gt.7550.0.and.ldb(i).lt.7650) lstrong(i)=0.0

        end if
        end if

c        if (telescope.eq."COSMIC") then
c           lstrong(i)=lstrong(i)*exp(-0.5*(ldb(i)-6000.0)**2/2000**2)
c        end if

        if (dichroic) then
         if (redside.and.ldb(i).lt.6500.0) then
                lstrong(i)=lstrong(i)*(ldb(i)/7000.0)**(10)
         end if
         if (.not.redside.and.ldb(i).gt.6500.0) then
                lstrong(i)=lstrong(i)*(ldb(i)/7000.0)**(-10)
         end if
        end if



c        if ( lname(i)(1:5).eq.'[O I]') then
c         write(6,*) ldb(i),lstrong(i),dichroic,redside
c        end if

        end do
  600        continue
        close(2)

        nlines = i
        lstrong(0) = -9
        lname(0) = 'xxxxxxxxxxxxxxxxxxxx'

        write(6,*) 'There were ',nlines, ndata, ip

ccccccccc
        if ( (arc.and.nimages.eq.1).or.(.not.arc.and.nimages.eq.1.and.
     &        zero.lt.nbig).or.
     &        	(arc.and.nimages.gt.1.and.zero.gt.nbig)) then

        write(6,601)
  601        format('Starting Line Correlation')

c        set up the reference array

        do ik=1,nlines
        xc(ik) = ldb(ik)
        end do
        call quick( xc, nlines, index)

        do j=1,50
        disco(j) = 0.0
        end do
        wfftonly=.false.
        if (wfft .gt. 5) then
          wfftonly=.true.
          wfft = wfft - 5
        end if
        ncents = min(wfft,5)
        nbs(1) = 1
        nbs(2) = 2
        nbs(3) = 4
        nbs(4) = 8
        nbs(5) = 16
        nbs(ncents+1) = nbs(ncents)
c        ncents = ncents + 1

        do icent=1,ncents

        write(6,*) 'FFT level=',icent
        write(4,*) 'FFT level=',icent

        if (icent.eq.1) then
         mfit = 2
         nuse = ndata
         interval = nuse
         jn0 = 1
         jn1 = 1
        else
         nuse = int(float(ndata)/float(nbs(icent)))
         interval = nuse/2
         jn0 = 1
         jn1 = ndata-interval
        endif

        nbin = nuse
        nhalf = nbin / 2
        nhalf1 = nhalf + 1

c        do nb=1,nbs(icent)
        nb = 0
        do inb=jn0,jn1,interval
        nb = nb + 1
        n0 = (nb-1) * interval + 1
        n1 = min(ndata,n0+nuse-1)
        write(6,*) icent,nb,n0,n1

        fmax = 0.0
        if ( icent.eq.1) then
         rwav0 = xc(1)
         rwav1 = xc(nlines)
          if ( slice.and.nsname.ne.0.0.and.(.not.debug(2).and.
     &          .not.debug(3).and..not.debug(4).and..not.debug(1))) then
           rwav0 = lamb1
           rwav1 = lamb2
          end if
          if ( lamb2.gt.lamb1) then
           rwav0 = lamb1
           rwav1 = lamb2
          end if
         rwavmax = rwav0
         rdelta = 50.0
        else
c          rwav0 = rwav-disp0*ndata/2.0 + n0*disp0
c          rwav1 = rwav-disp0*ndata/2.0 + n1*disp0
c          rdelta = 2*nuse*disp0
          bmid = ndata/2.0 + 0.5
          bmid0= ndata/2.0
          xn0 = (n0-bmid)/bmid0
          rwav0=0.0
          call fleg1( xn0, duml, mfit)
          do mf=1,mfit
           rwav0 = rwav0+duml(mf)*disco(mf)
          end do
          xn1 = (n1-bmid)/bmid0
          rwav1=0.0
          call fleg1( xn1, duml, mfit)
          do mf=1,mfit
           rwav1 = rwav1+duml(mf)*disco(mf)
          end do
          rdelta = 2*(rwav1-rwav0)
          write(6,*) rwav0,rwav1,rdelta
        end if

        do i=n0,n1
        sraw(i-n0+1) = raw(i)
        sraw2(i-n0+1) = raw(i)
        end do

        do i=1,2048
        tcen(i,1) = 0.0
        tcen(i,2) = 0.0
        end do

c        call quick(sraw,nuse,index)
        do i=1,n1-n0+1
         nij=5*int(fwhm1)
         do j=1,nij
          ij = min(max(1,i+j-nij/2),n1-n0)
          qsraw(j)=sraw2(ij)
c          if (icent.eq.3.and.nb.eq.1) then
c           write(6,*) i,ij,sraw2(ij)
c          end if
         end do
        call quick(qsraw,nij,index)
        qs=qsraw(int(0.05*nij)+1)
        sraw(i) = sraw2(i) - qs
c        sraw(i) = sraw2(i) - sraw2(index(0.05*nuse))
c        if (icent.eq.3.and.nb.eq.1) then
c         write(6,*) i,nij,sraw(i),qs
c        end if
        end do

        icen = 0
        tlast = 0.0
        do rwavu = rwav0, rwav1, rdelta

        t2 = log10(float(nuse)-1)/log10(2.0)
        nt2 = int(t2)+2
        nbin = 2**(nt2)

        if ( icent.gt.1) then
        end if
        do i=nuse+1,nbin
        sraw(i) = 0.0
        end do

        nhalf = nbin / 2
        nhalf1 = nbin / 2 + 1

        if (icent.eq.1) then
         z0 = rwavu - disp0*nuse/2.0
         z1 = z0 + (nuse-1) * disp0
         do i=1,nbin
         llam(i) = z0 + (i-1) * disp0
         end do
        else
         z0 = rwav0
         z1 = rwav1
         do i=1,nbin
          xxx = (i+n0-bmid)/bmid0
          llam(i)=0.0
          call fleg1( xxx, duml, mfit)
          do mf=1,mfit
           llam(i)= llam(i)+duml(mf)*disco(mf)
          end do
         end do
        end if

        do i=1,nbin
        lspec(i) = 0.0
        srcor(i) = 0.0
        end do

        totflux = 0.0
        do i=1,nlines
         if ( ldb(i).ge.z0.and.ldb(i).le.z1.and..not.kill(i)) then
          if ( ldb(i).ge.rwav0.and.ldb(i).le.rwav1) then
           r = (ldb(i) - z0)/disp0 + 1
           dispu = disp0*(1.0+0.0*((r-nuse/2.0)/(nuse/2.0)))
           r = (ldb(i) - z0)/dispu + 1
           i0 = r - 10
           i1 = r + 10
           sss = lstrong(i)
            do j=i0,i1
             lspec(j)=lspec(j)+max(profile(fwhm1,j,r)*(3.0*sss),0.)
            end do
          end if
         end if
        end do

        bmid = nbin/2.0 + 0.5
        bmid0= nbin/2.0
        do jj=nuse,nbin
         lspec(jj) = 0.0
         sraw(jj) = 0.0
        end do

        call fft( nbin, sraw, (nbin/2+1), sfraw)
        call invfft( nhalf1, sfraw, nbin, sraw)
        call fft( nbin, lspec, (nbin/2+1), fref)
        call correlate( nhalf1, fref, sfraw, sfcor)
        call invfft( nhalf1, sfcor, nbin, srcor)

        call shiftme( 1, nbin, srcor)

        yp1 = srcor(1) - srcor(nbin)
        ypn = -yp1
        call spline(xco, srcor, nbin, yp1, ypn, y2)

        imax=1
        if ( icent.lt.2) then
         ii1 = 10
         ii2 = nbin-20
        else
         ii1 = nbin/2-nbin/8
         ii2 = nbin/2+nbin/8
        end if
        do j=ii1,ii2
         if ( srcor(j).gt.srcor(imax)) imax = j
        end do

        gimme0 = getmax(nbin,xco,imax,srcor,y2,f)
        gimme = gimme0 - nbin/2.
        call thresh(sraw,nuse,1,nuse,s2n(nb),skk, com(nb))
        if (icent.eq.1) then
         rwavc=int((rwavu+gimme*disp0)/
     &        	(1.5*fwhm1)+0.5)*(1.5*fwhm1)
c         write(6,*) rwavc, gimme0, f
        else
         ddr = llam(nbin/2)- llam(nbin/2-1)
         rwavc = (rwav0+com(nb)*ddr) + ddr*gimme
         gimmes(nb) = gimme
         com(nb) = (n0+n1)/2.0
         comr(nb) = llam(int(nuse/2+gimme))
         midp(nb) = (n0+n1)/2.0
         write(4,*) n0,n1,gimme0,gimme,rwavc,com(nb),comr(nb),s2n(nb)
         write(6,*) n0,n1,gimme0,gimme,rwavc,com(nb),comr(nb),s2n(nb)
        end if

        ich=0
        do jc=1,icen
        if ( rwavc.eq.tcen(jc,1)) then
        tcen(jc,2) = tcen(jc,2) + 1
        ich = 1
        end if
        end do

        if ( ich.eq.0) then
        icen = icen + 1
        tcen(icen,1) = rwavc
        tcen(icen,2) = 1
        end if

        if ( f.gt.fmax) then
         fmax = f
         rwavmax = rwavc
        end if

        if (icent.eq.1) then
         disp1 = disp0
        end if

        end do

        testwav0 = rwav0
        testwav1 = rwav1
        if ( slice.and.nsname.ne.0.0.and.(.not.debug(2).and.
     &        .not.debug(3).and..not.debug(4).and..not.debug(1))) then
         difl = 0.05*(lamb2 - lamb1)
         testwav0 = lamb1+difl
         testwav1 = lamb2-difl
        end if
        if ( lamb2.gt.lamb1) then
         difl = 0.05*(lamb2 - lamb1)
         testwav0 = lamb1+difl
         testwav1 = lamb2-difl
        end if
        tnwav = 0
        do i=1,icen
        write(6,*) icent,nb,(tcen(i,j),j=1,2),tnwav,tcen(iwav,1)
        write(4,*) icent,nb,(tcen(i,j),j=1,2),tnwav,tcen(iwav,1)
        if ( tcen(i,2).gt.tnwav.and.
     &        	(tcen(i,1).gt.testwav0.and.
     &        	tcen(i,1).lt.testwav1)) then
         iwav = i
         tnwav = tcen(i,2)
        end if
        end do

        nbin = nuse
        bmid = nbin/2.0+0.5
        bmid0= nbin/2.0

        write(6,*) rwavmax,fmax,tnwav
        write(4,*) rwavmax,fmax,tnwav

c        rwav = rwavmax

        if (icent.eq.1) then
         rwav = tcen(iwav,1)
         write(6,*) rwav
         write(4,*) rwav
         write(4,*) rwav,disp1
        end if
        write(6,*)

        end do

        if (icent.eq.1) then
         disco(1) = rwav
         disco(2) = disp1*ndata/2
         nterms=2
        else
         bmid = ndata/2.0+0.5
         bmid0= ndata/2.0

c         if (nbs(icent).gt.2) then
         if (nb.gt.2) then
          fbord = nbin/2
         else
          fbord = nbin
         end if

         s2nlim = 1.5*nb
         if (d0.gt.4.5) s2nlim=100.0
         if (arc.and.nimages.eq.1) s2nlim=1e3
c         s2nlim = 1.5*nbs(icent)
c         s2nlim = 4.0

         kn = 0
         write(6,*) bmid, bmid0, nbin, fbord
c         do jn = 1,nbs(icent)
         do jn = 1,nb
          if (s2n(jn).gt.s2nlim.and.comr(jn).gt.0.0.and.
     &          abs(gimmes(jn)).lt.fbord) then
            kn = kn + 1
c            dumx(kn) = (midp(jn) - bmid)/bmid0
            dumx(kn) = (com(jn) - bmid)/bmid0
            dumy(kn) = comr(jn)
            dumw(kn) = 1.0
c            dumw(kn) = 1.0/s2n(jn)
          end if
         write(6,*) jn, midp(jn), com(jn), comr(jn), gimmes(jn),
     &               s2n(jn)
         write(4,*) jn, midp(jn), com(jn), comr(jn), gimmes(jn),
     &               s2n(jn)
         end do

         if (icent.eq.2) mfit=min(min(3,nbig),kn)
         if (icent.eq.3) mfit=min(min(4,nbig),kn)
         if (icent.eq.4) mfit=min(min(5,nbig),kn)
         if (icent.eq.5) mfit=min(min(5,nbig),kn)
c         if (icent.eq.ncents) mfit=min(min(4,nbig),kn)
         nterms = mfit
         do ij=1,mfit
          lista(ij) = ij
         end do
         do ij=mfit+1,nterms
          lista(ij) = 0
         end do

         write(6,*) (dumx(j),j=1,kn)
         write(4,*) (dumx(j),j=1,kn)
         write(6,*) (dumy(j),j=1,kn)
         write(4,*) (dumy(j),j=1,kn)
         call lfit( dumx, dumy, dumw, kn, disco, nterms, lista, mfit,
     &                        covar, ncvm, chisq, fleg1)
         write(6,*) (disco(j),j=1,mfit)
         write(4,*) (disco(j),j=1,mfit)
         do jn=1,kn
          yy=0.0
          call fleg1( dumx(jn), duml, mfit)
          do mf=1,mfit
           yy = yy+duml(mf)*disco(mf)
          end do
          if (icent.eq.ncents) then
          write(24,*) jn,dumx(jn),dumy(jn),yy,yy-dumy(jn)
          end if
          write(4,*) jn,dumx(jn),dumy(jn),yy,yy-dumy(jn)
          write(6,*) jn,dumx(jn),dumy(jn),yy,yy-dumy(jn)
         end do
        end if

        end do

        write(4,*) (disco(j),j=1,nterms)
c        if (.not.arc.and.wfft.gt.3) goto 1290
c        stop

c        if ( red) rwav=rwavmax

c        disco(1) = rwav
c        disco(2) = 0.96*disp1*bmid0
c        disco(3) = 0.037897*disco(2)
c        disco(4) = -0.0116*disco(3)
c        disco(5) = -0.22*disco(4)

        else if ( nimages.gt.1) then

         nterms = nbig
         rwav = disco(1)
         disp1 = disco(1)/(ndata/2.)

        end if
ccccccccc

        if (wfftonly) goto 1290

        nbin = ndata
        bmid = nbin/2.0+0.5
        bmid0= nbin/2.0

        do i=1,ndata
        sraw(i) = raw(i)
        end do

        write(13,*) 'starting frame',nimages
        call quick( sraw, ndata, index)
        perc = 0.02

        nbot = 0
        bottom = 0.0
        do i=ndata*perc/2,ndata*perc
        nbot = nbot + 1
        bottom = bottom + sraw(i)
        end do
        bottom = bottom/(nbot)

        b2 = 0.0
        do i=1,ndata*perc
        b2 = b2 + (sraw(i) - bottom)**2
        end do

        b2 = sqrt(b2 / (nbot - 1.))
        b2 = sqrt(bottom)

        if (bth0.eq.0) then
         bth = 5.0
c         if ( rwav.gt.7500) bth = bth * 1.0
c         if ( rwav.lt.7300) bth = bth * 2. / 3.
c         bth=bth/2d0
c         write(6,*) bth*b2+bottom
c         if ( arc.and.nimages.gt.1) bth = 2.
         if ( arc.and.nimages.eq.1) bth = 2.
c         if ( .not.arc.and.rwav.lt.7300) bth = 7.5
         if (d0.gt.3) bth = bth*2
         if (telescope.eq."COSMIC") bth = 2.0
         if (telescope.eq."KECK") bth = 4.0
        else
         bth = bth0
c         if (nimages.gt.1) bth = 2
        end if
        write(6,*) 'BTH=',bth
        write(4,*) 'BTH=',bth

        nsearch = 0
        fwhm16= max(fwhm1,2.0)
  602        continue
        nsearch = nsearch + 1

c        ndata1 = ndata-1

c        do thr=hi,bth,-1.0

        nfw = 2
        do ifw=1,nfw
        nm = 0

        do id=ndata,1,-1

        thr = bth
        if ( (.not.arc.or.(arc.and.nimages.gt.1)).and.
     &        	index(id).lt.pix5577) thr = 2

        nib = 10 * fwhm1

        ihloc = index(id)
        ib1 = max((ihloc-nib/2),1)
        ib2 = min((ib1+nib),ndata)
        ib1 = ib2 - nib

        ib1 = max((ihloc-nib/2),1)
        nib = ib2-ib1+1
        mib = nib/20+1

        do j=1,nib
        subsky(j) = raw(ib1+j-1)
        end do

        call quick( subsky, nib, index)
        bottom = subsky(mib)
        hh = sqrt(max((raw(index(id)) - bottom),0.0))
c        write(6,*) index(id),raw(index(id)), bottom, thr, hh
        if (hh.gt.thr.and.index(id).gt.5.and.index(id).lt.ndata-5) then

        nm = nm + 1

        dmax(nm) = index(id)

        tmax(nm) = getlinemax(ndata,xco,dmax(nm),raw,draw,ymax(nm),
     &        	fwhm0,bottom)
        h(nm) = sqrt(max((raw(index(id)) - bottom),0.0))
c        write(6,*) nm,index(id),raw(index(id)), bottom, thr, hh,fwhm0

        fwhm(nm) = fwhm0

        if ( h(nm).lt.thr) then
         nm = nm - 1
         goto 210
        end if

        s = abs( dmax(nm) - tmax(nm) )
        if ( s.gt.0.6*fwhm16) then
         nm = nm - 1
         goto 210
        end if

        do km=1,nm-1
        s = abs( tmax(km) - tmax(nm) )
        if ( s.le.0.6*fwhm16) then
         nm = nm -1
         goto 210
        end if
        end do

        if ( fwhm(nm).gt.fwhm16.or.fwhm(nm).lt.1.5) then
         nm = nm - 1
         goto 210
        end if

        write(6,'(i4,4f9.3)') nm,tmax(nm),h(nm),fwhm(nm),bottom
        write(13,'(i4,4f9.3)') nm,tmax(nm),h(nm),fwhm(nm),bottom

c        if ( nm.eq.nlines) goto 211

  210        continue
        end if

        end do
c        ndata1 = index(id)

        nmold = nm + 1

  211        continue
        hi=h(1)

        write(6,*) 'Peaks ',nm
        write(4,*) 'Peaks ',nm
        if ( nm.eq.0) then
         if (fwhm16.lt.25) then
           fwhm16 = fwhm16 * 1.5
           write(6,*) 'upping the fwhm by x1.5 to',fwhm16
           goto 602
         else
           write(6,*) 'no lines'
           err=2
           return
         end if
        end if

        if ( ifw.le.nfw-1) then
        do jfw=1,nm
        xc(jfw) = fwhm(jfw)
        reject(jfw) = .false.
        end do
        call quick(xc,nm,index)
        if (nm.gt.1) fwhm16 = xc(nm/2) * 2.0
        write(4,*) fwhm1, fwhm16
        write(6,*) fwhm1, fwhm16
        end if

        minlines = min(20,int(20*disp0+0.5))
        if (.not.arc.and.rwav.lt.7300) minlines = int(nbig)
        if (arc.and.nimages.eq.1) minlines = int(nbig)
        if (arc.and.nimages.gt.1) minlines = zero
c        if ( nm.lt.20*disp0) then
        write(6,*) 'nm,minlines',nm,minlines
c        if (.not.arc) then
        if ( nm.lt.minlines) then
         fwhm16 = fwhm16 * 1.5
         write(6,*) 'upping the fwhm by x1.5 to',fwhm16
         goto 602
        end if
c        end if

        end do

        write(6,*) (disco(j),j=1,nterms)
        write(4,*) (disco(j),j=1,nterms)

        do i=1,nm
c        wmax(i) = rwav + (tmax(i)-bmid)*disp1
        xx = (tmax(i) - bmid)/bmid0
        call fleg1(xx,pt,nterms)
        wmax(i) = 0.0
        do j=1,nterms
        wmax(i) = wmax(i) + pt(j)*disco(j)
        end do
c        write(6,*) i, tmax(i), wmax(i)
        end do

c        nlinmin = nbig + 1
c        if ( nimages.gt.1.and.zero.lt.nbig) nlinmin=zero+1

        if ( nimages.eq.1.and.nsearch.lt.4.and.nm.lt.nlinmin) then
         bth = bth / 2.
         goto 602
        end if

        if ( nimages.eq.1.and.nm.le.1) then
         write(6,*) 'Stopping. No lines. Try fixing fwhm',fwhm16
         err = 2
         stop
        end if
        if ( nimages.gt.1.and.nm.lt.1) then
         write(6,*) 'Stopping. No lines'
         err = 3
         return
        end if

c        if ( arc.and.nimages.eq.1) then
c         if ( nm.gt.max(100,nlines)) nm=nlines
c        end if

        do i=1,nm
        temp1(i) = tmax(i)
        ht(i) = h(i)
        fwhm2(i) = fwhm(i)
        end do

        tlobt = 0.0
        if ( arc.and.nimages.eq.1) tlobt = 0.0
        hibt = 1
        nmax = 4
        if ( nimages.eq.1) then
         tol = 10*disp0
         if ( arc) tol=40*disp0
         if ( .not.arc.and.dichroic.and..not.redside) tol=40*disp0
        else
         tol = 5*disp0
        end if
  305        continue

  306        continue

c        write(6,*) tol

        do i=1,nm
        used(refmax(i)) = .false.
        refmax(i) = 0
        end do

        if (arc.and.nimages.eq.1) then
         tktop = 10e4
        else
         tktop = 1e3
        end if

        do 240 i=1,nm

        refmax(i) = 0
        smin = 9999.9
        do 230 tk=10e4,tlobt,-100

        do 220 j=1,nlines
        if ( .not.isaqr(j)) then
         if ( lstrong(j).ge.tk.and..not.used(j).and.ht(i).gt.bth) then
          s = abs(wmax(i) - ldb(j))
          if ( s.lt.smin.and.s.lt.tol) then
           if ( lstrong(j).ge.lstrong(refmax(i))) then
            used(refmax(i)) = .false.
            smin = s
            refmax(i) = j
            used(j) = .true.
           end if
          end if
         end if
        end if
  220        continue

  230        continue

  240        continue

        do i=1,nm
        write(4,212) i,lstrong(refmax(i)),
     &        temp1(i), wmax(i), ldb(refmax(i)), ht(i), lname(refmax(i))
        diff = abs(wmax(i) - ldb(refmax(i)))
        if ( temp1(i).gt.ndata-fwhm16) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        if ( temp1(i).lt.fwhm16) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        if ( diff.gt.tol) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        if ( fwhm(i).gt.1*fwhm16.and.arc) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        end do
        write(4,*)

        call quick( temp1, nm, spord)
        call rectfy( ht, nm, spord, hold)
        call irectfy( refmax, nm, spord, ihold)
        call rectfy( wmax, nm, spord, hold)
        call rectfy( fwhm2, nm, spord, hold)

        do i=1,nlines
        used(i) = .false.
        end do
        do i=1,nm
        used(refmax(i))=.true.
        write(13,*) i,temp1(i),ht(i),wmax(i)
        end do

        nr = 0
c        do nbt=10,-5,-0.5
        nbt = 0
        do i=1,nm
        if ( refmax(i).gt.0.and.lstrong(refmax(i)).gt.nbt) then
        nr = nr + 1
        fx(nr) = (temp1(i)-bmid)/bmid0
        fy(nr) = ldb(refmax(i))
        fw(nr) = 1.0
        fi(nr) = i
        write(6,212) nr, lstrong(refmax(i)), temp1(i), wmax(i),
     &      ldb(refmax(i)), ht(i), lname(refmax(i)), fwhm2(i)
        write(4,212) nr, lstrong(refmax(i)), temp1(i), wmax(i),
     &       ldb(refmax(i)), ht(i), lname(refmax(i)), fwhm2(i)
        end if
c        end do

  212        format(i5,f7.1,4f11.4,4x,a10,f7.2)

        end do
        write(4,*)

c     &        (red.and.arc.and.nimages.gt.1.and.zero.eq.nbig)).and.
        if (graname.eq.'xcosmic'.and.((.not.arc.and.nimages.eq.1).or.
     &        (arc.and.nimages.gt.1.and.zero.eq.nbig)).and.
     &        (nr.lt.999).and.fwhm1.le.12) then
        write(6,*) 'Checking for P1 spacings'

        irr = 0
        do rk=10.0,0.0,-0.5
        do i=1,nr
        if ( lstrong(refmax(fi(i))).ge.rk.and.ht(fi(i)).gt.bth) then
        irr = irr + 1
        rr(irr) = i
        if (irr.eq.2) then
         write(6,*) ldb(refmax(fi(rr(1)))), ldb(refmax(fi(rr(2))))
         sdif = abs(temp1(fi(rr(2)))-temp1(fi(rr(1))))
         if ( sdif.lt.250) then
         irr = 1
         else
         goto 1008
         end if
        end if
        end if
        end do
        end do

 1008        continue

c        wave0 = ldb(refmax(fi(rr(1))))
c        wave1 = ldb(refmax(fi(rr(2))))

        maxy0 = fi(rr(1))
        maxy1 = fi(rr(2))
        wave0 = wmax(maxy0)
        wave1 = wmax(maxy1)

        refm0 = refmax(maxy0)
        refm1 = refmax(maxy1)

        if ( lstrong(refmax(fi(rr(1)))) .gt.4) goodgrab=.true.
        write(6,*) lstrong(refmax(fi(rr(1)))),ldb(refmax(fi(rr(1))))
        write(6,*) lstrong(refmax(fi(rr(2)))),ldb(refmax(fi(rr(2))))
        write(6,*) 'goodgrab: ',goodgrab

        ww2 = 1.0
        if ( wave1-wave0.lt.500.0) ww2=2.0

        space = temp1(maxy1) - temp1(maxy0)
        write(6,*) temp1(maxy0), temp1(maxy1)
        write(13,*)

        do i=1,nm
c        if ( refmax(i).ne.refm0.and.refmax(i).ne.refm1) then
c        if ( refmax(i).ne.refm0) then
        used(refmax(i)) = .false.
        refmax(i) = 0
c        end if
        end do

        avs = 0.0
        navs = 0
        ssm = 0.0
        do ntri=1,3
        do mtri=1,3

        do i=1,nm-2
        tri1(i) = (temp1(i) - temp1(fi(rr(1)))) / space
        tri3(i,1,ntri,mtri) = (temp1(i)-temp1(fi(rr(1)))) / space
        tri3(i,2,ntri,mtri) = (temp1(i+ntri)-temp1(i)) / space
        tri3(i,3,ntri,mtri) = (temp1(i+mtri+ntri)-temp1(i+ntri)) / space
        write(13,*) temp1(i),(tri3(i,j,ntri,mtri),j=1,3)
        end do
        write(13,*)
        end do
        end do

        do i=1,ip-2
        tri2(i) = ldb(pind(i))
        tri4(i,1) = (ldb(pind(i)) - wave0)/(wave1-wave0)
        tri4(i,2) = (ldb(pind(i+1)) - ldb(pind(i)))/(wave1-wave0)
        tri4(i,3) = (ldb(pind(i+2)) - ldb(pind(i+1)))/(wave1-wave0)
c        tri4(i,1) = (ldb(pind(i)) - wave0)/(wave1-wave0)
c        tri4(i,2) = (ldb(pind(i+2)) - ldb(pind(i)))/(wave1-wave0)
c        tri4(i,3) = (ldb(pind(i+4)) - ldb(pind(i+2)))/(wave1-wave0)
        write(13,*) ldb(pind(i)),(tri4(i,j),j=1,3)
        end do
        write(13,*)

        p1tol = 0.02
        p3tol = 0.02

        p1tol = 0.010
        p3tol = 0.0005

        iold = 1
        jold = 1
 1007        continue
c        iold = 1
        do i=iold,ip-2
        if (.not.used(pind(i)).and.
     &        .not.used(pind(i+1)).and.
     &        .not.used(pind(i+2))) then
        s1min = 9.9
        s2min = 9.9
        s3min = 9.9
        if (ldb(pind(i)).gt.5000.0) then
        do j=jold,nm-2

        do ntri=1,3
        do mtri=1,3

        s1 = abs((tri3(j,1,ntri,mtri) - tri4(i,1))/tri4(i,1))
        s2 = abs((tri3(j,2,ntri,mtri) - tri4(i,2)))
        s3 = abs((tri3(j,3,ntri,mtri) - tri4(i,3)))
        if ( (s1.lt.s1min).and.(s1.lt.p1tol
     &        .and.s1.lt.s1min.and.s2.lt.s2min.and.s3.lt.s3min
     &        	.and.s2.lt.p3tol.and.s3.lt.p3tol)) then
           s1min = s1
           s2min = s2
           s3min = s3
           j1min = j
           j2min = j+ntri
           j3min = j+mtri+ntri
           used(pind(i)) = .false.
           used(pind(i+1)) = .false.
           used(pind(i+2)) = .false.
           triplet(pind(i)) = .false.
           triplet(pind(i+1)) = .false.
           triplet(pind(i+2)) = .false.
           refmax(j) = pind(i)
           refmax(j+ntri) = pind(i+1)
           refmax(j+mtri+ntri) = pind(i+2)
           attach(pind(i)) = j
           attach(pind(i+1)) = j+ntri
           attach(pind(i+2)) = j+mtri+ntri
           used(pind(i)) = .true.
           used(pind(i+1)) = .true.
           used(pind(i+2)) = .true.
           triplet(pind(i)) = .true.
           triplet(pind(i+1)) = .true.
           triplet(pind(i+2)) = .true.
           navs = navs + 3
c           jold = j + 3
c           iold = i + 3
c           goto 1007
        end if
        end do
        end do

        end do
        end if

        if ( used(pind(i)).and.used(pind(i+1)).and.used(pind(i+2))) then
        write(6,'(6f9.2,3f8.5)') ldb(refmax(j1min)),ldb(refmax(j2min)),
     &        	ldb(refmax(j3min)),
     &        temp1(j1min), temp1(j2min), temp1(j3min),
     &        s1min, s2min, s3min
        write(4,'(6f9.2,3f8.5)') ldb(refmax(j1min)),ldb(refmax(j2min)),
     &        	ldb(refmax(j3min)),
     &        temp1(j1min), temp1(j2min), temp1(j3min),
     &        s1min, s2min, s3min
        end if

        end if
        end do

        navs = 0
        do j=1,nm
        if ( attach(refmax(j)).ne.j) refmax(j) = 0
        if ( refmax(j).gt.0) navs = navs + 1
        end do

        if ( navs.eq.0) write(4,*) ' nothing'
        if ( navs.eq.0) write(6,*) ' nothing'
        if ( navs.gt.0) write(6,*) 'navs',navs

c        if ( navs.eq.0) then
        if ( int(ldb(refm0)).eq.5577.or.
     &        	int(ldb(refm0)).eq.6300.or.navs.lt.9) then
        if ( .not.triplet(refm0)) then
        used(refm0) = .true.
        refmax(maxy0) = refm0
        write(6,*) ldb(refmax(maxy0))
        navs = navs + 1
        end if
        end if
        if ( int(ldb(refm1)).eq.5577.or.
     &        	int(ldb(refm1)).eq.6300.or.navs.lt.9) then
        if ( .not.triplet(refm1)) then
        used(refm1) = .true.
        refmax(maxy1) = refm1
        navs = navs + 1
        write(6,*) ldb(refmax(maxy1))
        end if
        end if
c        end if

        nterms=min(nterms,nbig)
c        if (.not.arc.and.navs.lt.5) nterms=min(3,nbig)
c        if (.not.arc.and.navs.gt.5) nterms=min(3,nbig)
c        if (.not.arc.and.navs.gt.15) nterms=min(4,nbig)
        else if ( arc.and.nimages.eq.1) then
         nterms = nterms
c         nterms = 3
        else if ( arc.and.nimages.gt.1) then
         nterms = nbig
        else if ( .not.arc.and.nimages.gt.1) then
         nterms = nbig
        end if

        tol = 20.0

        getmore = .true.
        nr = 0
        lobt = 2
        idiot = 0
        nbig = max(3,nbig)
        scat = 99.0
        www = 20.0
        if ( arc.and.nimages.eq.1) www=1.5*www
        if ( .not.arc.and.nimages.gt.1.and.zero.lt.nbig) then
        	www= 2
        end if
c        bth = max(bth -10.,10.0)
c        if ( .not.arc.and.rwav.lt.7200) www=www*3
        if ( arc.and.zero.eq.nbig.and.
     &        	nimages.gt.1.and.rwav.lt.7200) www=3.0
c        if ( red.and.arc.and.nimages.eq.1) tol=5.0
c        if ( red.and.arc.and.nimages.eq.1) www=5.0
c        if ( .not.red.and.arc.and.nimages.eq.1) tol=10.0
c        if ( .not.red.and.arc.and.nimages.eq.1) www=5.0
 1010        continue
        already = .false.
        idiot = idiot + 1
        istupid = 0
 1011        continue

c        write(6,*) disco
        minpix = 4096
        maxpix = 0
        nr = 0
        s = 0.0
        do i=1,nm
c        write(6,*) i,temp1(i),refmax(i),ldb(refmax(i))
        if ( refmax(i).gt.0) then
         s = abs(wmax(i) - ldb(refmax(i)))
         if (s.lt.www*disp0.and..not.reject(i)
     &                      .and..not.isaqr(refmax(i))) then
          nr = nr + 1
        if ( temp1(i).lt.minpix) minpix=temp1(i)
        if ( temp1(i).gt.maxpix) maxpix=temp1(i)
          fx(nr) = (temp1(i)-bmid)/bmid0
          fy(nr) = ldb(refmax(i))
          fi(nr) = i
          fw(nr) = fwhm2(i)/ht(i)
         else
          used(refmax(i)) = .false.
          refmax(i) = 0
         end if
        end if
        end do

        if ( arc.and.nimages.eq.1.and.nr.lt.4) nterms = 2
        if ( .not.arc) then
         if ( nr.lt.nterms) nterms = nr-1
         if ( nr.lt.4) nterms = 2
         if ( nimages.eq.1.or.zero.eq.nbig) then
         if ( idiot.gt.3) then
          if ( nr.ge.4.and.scat.gt.1.0*wtol.and.fwhm1.lt.12) nterms = 3
          if ( nr.gt.2*nterms.and.scat.lt.1.0*wtol)
     &        		nterms = min(nterms+1,nbig)
c         else
c          if ( nr.gt.6) nterms = min(4,nbig-1)
         end if
        else if ( nimages.gt.1.and.zero.lt.nbig) then
         nterms = nbig
        end if
        end if
        if ( arc.and.nimages.gt.1.and.zero.eq.nbig) then
         nterms = min(nterms,nr-1)
         if ( nr.lt.4) nterms = 2
        end if
        if ( istupid.gt.4.and.nr.gt.7) nterms = min(nterms+1,nbig)

        mfit = nterms
        if (arc.and.nimages.gt.1.and.zero.lt.nbig) mfit = min(zero,nr-1)
        if ( .not.arc.and.nimages.gt.1
     &          .and.zero.lt.nbig.and.zero.gt.0) then
                mfit = min(zero,nr-1)
                nterms = nbig
        end if
        if ( mfit.lt.1) mfit = 1

        do ij=1,mfit
        lista(ij) = ij
        end do
        do ij=mfit+1,nterms
        lista(ij) = 0
        end do

        write(6,*) ' before lfit', nr

        if (nr.eq.0) then
         write(6,*) ' fuck: went horribly wrong!'
         return
        end if

        nrold = nr

        call lfit( fx, fy, fw, nr, disco, nterms, lista, mfit,
     &                        covar, ncvm, chisq, fleg1)
        nterms_old = nterms
        already = .true.
        numlfit = numlfit + 1

        if ( istupid.eq.0.0) then
        oldscat2 = oldscat
        oldscat = scat
        end if

        write(4,*) idiot,www,chisq,lobt
        write(4,*) (disco(i),i=1,nterms)
        write(4,*)
        write(6,*) idiot,www,chisq,lobt
        write(6,*) (disco(i),i=1,nterms)

        scat = 0.0
        scat22 = 0.0
        sumwt = 0.0
        sumwt22 = 0.0
        nsm=0
        do i=1,nm
         xx = (temp1(i) - bmid)/bmid0
         call fleg1(xx,pt,nterms)
         wmax(i) = 0.0
         do j=1,nterms
         wmax(i) = wmax(i) + pt(j)*disco(j)
         end do

         if ( refmax(i).gt.0) then
          nsm= nsm+1
          scat = scat + (wmax(i)-ldb(refmax(i)))**2*ht(i)**2/fwhm2(i)
          scat22 = scat22 + (wmax(i)-ldb(refmax(i)))**2
          sumwt = sumwt + ht(i)**2/fwhm2(i)
          sumwt22 = sumwt22 + 1
          dev(i) = abs(wmax(i) - ldb(refmax(i)))
c          smdev(nsm) = 1 / ht(i)
          smdev(nsm) = dev(i)
          rq(nsm) = i
          write(6,291) temp1(i),ht(i),wmax(i),
     &         	ldb(refmax(i)),fwhm2(i),lname(refmax(i))(1:12)
          write(4,291) temp1(i),ht(i),wmax(i),
     &         	ldb(refmax(i)),fwhm2(i),lname(refmax(i))(1:12)
  291         format(4f10.3,f7.2,3x,a12)
         end if
        end do

        scat = sqrt(scat/sumwt)
        scat22 = sqrt(scat22/sumwt22)
        write(4,289) scat, scat22
        write(6,289) scat, scat22
  289        format('RMS= ',2f15.10,1x,'(weighted,unweighted)')

        if ( idiot.gt.100.and.fwhm1.le.12.and.
c     &        .not.arc.and.
     &        	(scat.gt.0.5*fwhm16*disp0/11.75
     &          .or.scat22.gt.0.5*fwhm16*disp0/11.75)) then
        nreject = 0
        call quick(smdev,nsm,index)
        call irectfy( rq, nsm, index, ihold)
        do i=nsm,1,-1
        if ( abs(smdev(i)).gt.(disp0*fwhm16/2.35/5.0)
     &                 .and.abs(smdev(i)).gt.scat) then
c        if ( abs(dev(i)).gt.(disp0*fwhm16/2.35/5.0)
c     &          .and.abs(dev(i)).gt.scat) then
          used(refmax(rq(i))) = .false.
          refmax(rq(i)) = 0
          reject(rq(i)) = .true.
          nreject = nreject + 1
          goto 1013
         end if
        end do
 1013        continue
        if ( nreject.ge.1) then
         write(6,*) 'nreject=',nreject
         write(4,*) 'nreject=',nreject
         istupid = istupid + 1
         if (istupid.gt.7) goto 1012
         goto 1011
        end if
 1012        continue
        end if

        do i=1,nm
        if (refmax(i).gt.0) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        end do

c        write(6,*) scat, oldscat, oldscat2
c        write(6,*) 'here i am'

        if ( idiot.gt.3.and.nr.gt.nterms+1.and.nimages.ge.1
     &        .and.(scat.lt.wtol.or.scat.eq.oldscat
     &        .or.(scat.eq.oldscat2.and.repeat.ge.4))) then
        	nterms = min(nterms + 1,nbig)
        	www = max(www-1,5.0)
        	if (nterms.ne.nterms_old) already = .false.
        end if

        if (idiot.gt.1.and.scat.lt.oldscat) tol = max(wtol*2,(1.+scat))
        if (idiot.gt.1.and.scat.gt.oldscat.and.scat.gt.wtol.and..not.
     &        	(nimages.gt.1.and.zero.lt.nbig)) then
        	tol = min(tol+10.0,20.0)
        end if

        if (idiot.gt.1.and.scat.lt.wtol) then
         lobt=max(lobt-1,-1)
         tol = max((1.+scat),wtol*2)
        end if

        if (idiot.eq.1.and.scat.lt.wtol) then
         tol = max((1.+scat),tol/2.0)
         if ( nr.gt.200*nterms) then
          nterms = min(nterms + 1,nbig)
        	if (nterms.ne.nterms_old) already = .false.
         end if
        end if

        if ( scat.eq.oldscat) then
        repeat = repeat + 1
        else
        repeat = 0
        end if

c        write(6,*) 'here i am'

        if ( repeat.gt.8) then
         if ( nr.gt.2*nterms) then
          nterms = min(nterms + 1,nbig)
          if (nterms.ne.nterms_old) already = .false.
         end if
        end if

        if ( scat.gt.wtol.and.nterms.eq.nbig) lobt = min(lobt + 1,2)
c        if ( (idiot.le.2.or.scat.gt.1.0).and.arc) lobt = 4

        if ( idiot.gt.0) then

        write(6,*) 'nlines ', nm
        getmore = .false.
        if ( getmore) then

        if ( nr.lt.5) bth = bth / 2.
        if ( nr.gt.2*nterms) bth = max(bth/2.,5.0)
        if ( nr.lt.5) thr = bth 
        if ( nr.gt.2*nterms) thr=bth
        thr = max(thr,5.0)
c        write(6,*) thr

c        do thr=bthold,bth,-1

c        write(6,*) 'here i am'

        do id=ndata,1,-1
c        write(6,*) id

        if ( (.not.arc.or.(arc.and.nimages.gt.1)).and.
     &        	index(id).lt.pix5577) thr = bth
        if ( arc.and.nimages.eq.1) thr=bth

        nib = 10 * fwhm1

        ihloc = index(id)
        ib1 = max((ihloc-nib/2),1)
        ib2 = min((ib1+nib),ndata)
        ib1 = ib2 - nib

        nib = ib2-ib1+1
        mib = nib/20+1
c        write(6,*) nib,mib

        do j=1,nib
        subsky(j) = raw(ib1+j-1)
c        if ( id.eq.637) write(6,*) j,ib1+j-1,subsky(j)
        end do
        call quick( subsky, nib, index)
        bottom = subsky(mib)
c        write(6,*) bottom

        hh = sqrt(max((raw(index(id)) - bottom),0.0))
        if (hh.gt.thr.and.index(id).gt.5.and.index(id).lt.ndata-5) then

        do km=1,nm
        s = abs( tmax(km) - xco(index(id)) )
        if ( s.le.1*fwhm16) goto 1210
        end do

        nm = nm + 1

        dmax(nm) = index(id)

c        write(6,*) nm,dmax(nm)
        tmax(nm)=getlinemax(ndata,xco,dmax(nm),raw,draw,ymax(nm),fwhm0,
     &        	bottom)
c        write(6,*) nm,tmax(nm),bottom,fwhm0
        h(nm) = sqrt(max((raw(index(id)) - bottom),0.0))
        fwhm(nm) = fwhm0

c        write(6,*) nm,tmax(nm),bottom,fwhm0

        if ( h(nm).lt.thr) then
         nm = nm - 1
         goto 1210
        end if

        if ( fwhm(nm).gt.fwhm16.or.fwhm(nm).lt.2) then
         nm = nm - 1
         goto 1210
        end if

        s = abs( dmax(nm) - tmax(nm) )
        if ( s.gt.fwhm16) then
         nm = nm - 1
         goto 1210
        end if

        do km=1,nm-1
        s = abs( tmax(km) - tmax(nm) )
        if ( s.le.fwhm16) then
         nm = nm -1
         goto 1210
        end if
        end do

c        write(6,*) nm,tmax(nm),h(nm),fwhm(nm)

 1210        continue
        end if
c        end if

c        if ( red.and.nm.gt.(5*nterms)) goto 1250

        end do

c        end do

c        if ( nm.gt.nrold+10) nm=nrold+10
 1250        continue

c        write(6,*) 'here i am'
        getmore=.true.
c        if ( idiot.gt.2) getmore=.false.
c        if ( nr.gt.30) getmore = .false.
c        if ( bth.le.20.0.and.nr.gt.25) getmore = .false.
c        if ( bth.le.8.0) getmore = .false.

c        if ( red.and.bth.lt.8.) getmore = .false.
c        if ( red.and.nr.gt.(6*nterms)) getmore = .false.

        bthold = bth


        do i=1,nm
        temp1(i) = tmax(i)
        ht(i) = h(i)
        fwhm2(i) = fwhm(i)
        end do

c        call quick( temp1, nm, spord)
c        call rectfy( ht, nm, spord, hold)

        end if

        do i=1,nm
        xx = (temp1(i) - bmid)/bmid0
        call fleg1(xx,pt,nterms)
        wmax(i) = 0.0
        do j=1,nterms
        wmax(i) = wmax(i) + pt(j)*disco(j)
        end do
        write(4,*) i,temp1(i),wmax(i),ht(i)
        end do

        if ( idiot.eq.1) then
        do i=1,nlines
        if ( arc.and.nimages.eq.1) then
         lstrong(i) = 3.0
         if ( lname(i)(1:2).eq.'Hg'.and.ldb(i).lt.5500.0) lstrong(i)=3
         if ( lname(i)(1:2).eq.'Kr'.and.ldb(i).lt.5500.0) lstrong(i)=2
         if ( lname(i)(1:1).eq.'A'.and.ldb(i).lt.7500.0) lstrong(i)=3
         if ( ldb(i).gt.8000.0.and.lname(i).eq.'Kr') lstrong(i) = 4
         if ( int(ldb(i)).eq.8110.0) lstrong(i) = -5
         if ( int(ldb(i)).eq.8115.0) lstrong(i) = -5
         if (.not.hg.and.lname(i)(1:2).eq.'Hg') lstrong(i) = -99
         if (.not.hg.and.lname(i)(1:2).eq.'Kr') lstrong(i) = -99
         if (.not.cd.and.lname(i)(1:2).eq.'Cd') lstrong(i) = -99
         if (.not.zn.and.lname(i)(1:2).eq.'Zn') lstrong(i) = -99
         if (.not.ne.and.lname(i)(1:2).eq.'Ne') lstrong(i) = -99
         if (.not.ar.and.lname(i)(1:1).eq.'A') lstrong(i) = -99
         if (lname(i)(1:2).eq.'He') lstrong(i) = -99
        else
        lstrong(i) = -3.0
        used(i) = .false.

        if ( lname(i)(1:5).eq.'[O I]') then
         lstrong(i) = 2.0
         if (nint(ldb(i)).eq.5577) lstrong(i) = 50000.0
         if (nint(ldb(i)).eq.6300) lstrong(i) = 1000.0
         if (nint(ldb(i)).eq.6364) lstrong(i) = 500.0
        else if ( lname(i)(1:2).eq.'Hg') then
         if (int(ldb(i)).eq.3649) lstrong(i) = 100.0
         if (int(ldb(i)).eq.4046) lstrong(i) = 500.0
         if (int(ldb(i)).eq.4077) lstrong(i) = 500.0
         if (int(ldb(i)).eq.4358) lstrong(i) = 3000.0
         if (int(ldb(i)).eq.4916) lstrong(i) = 3000.0
         if (int(ldb(i)).eq.5460) lstrong(i) = 3000.0
        else if ( lname(i)(1:4).eq.'Na I') then
          if (telescope.eq."COSMIC") then
                lstrong(i) = 10000.0
                if (disp0.lt.3.0) lstrong(i) = 4000.0
          else if (telescope.eq."KECK") then
                lstrong(i) = 500.0
          end if
          if (disp0.gt.1.5) ldb(i) = 5892.935


c        if ( lname(i)(1:5).eq.'[O I]') then
c         lstrong(i) = 1.5
c         if (nint(ldb(i)).eq.5577) lstrong(i) = 4.0
c         if (nint(ldb(i)).eq.6300) lstrong(i) = 3.0
c         if (nint(ldb(i)).eq.6364) lstrong(i) = 3.0
c        end if
c        if ( lname(i)(1:2).eq.'Hg') then
c         lstrong(i) = 6
c        end if
c        if ( lname(i)(1:4).eq.'Na I') then
c         lstrong(i) = 5.0
c         ldb(i) = 5892.935
c        end if

        else if ( lname(i)(5:5).eq.'P') then
         if ( lname(i)(6:6).eq.'1') lstrong(i) = 2.0
         if ( lname(i)(6:6).eq.'2') lstrong(i) = 1.0
         ip1 = 0
         if ( lname(i)(7:7).eq.'(') then
         read(lname(i)(8:9),*) realip1
         ip1 = int(realip1)
         end if
         if ( lname(i)(8:8).eq.'(') then
         read(lname(i)(9:10),*) realip1
         ip1 = int(realip1)
         end if
         if (ip1.gt.8) lstrong(i) = -5
                if (graname.ne.'300gmm') then
         	lstrong(i) = lstrong(i)-ip1/3
          end if
        end if
        if ( (fwhm1.le.12).and.(d0.lt.4.5)) then
        if ( lname(i)(5:6).eq.'Q1') lstrong(i) = -2.0
        if ( lname(i)(5:6).eq.'Q2') lstrong(i) = -2.0
        if ( lname(i)(5:7).eq.'Q1(') lstrong(i) = -2.0
        if ( lname(i)(5:7).eq.'Q2(') lstrong(i) = -2.0
        if ( lname(i)(5:5).eq.'R') lstrong(i) = -2.0
        else
c        isaqr(i) = .false.
        if ( lname(i)(5:6).eq.'Q1') lstrong(i) = 2.0
        if ( lname(i)(5:6).eq.'Q2') lstrong(i) = 2.0
        if ( lname(i)(5:7).eq.'Q1(') lstrong(i) = 2.0
        if ( lname(i)(5:7).eq.'Q2(') lstrong(i) = 2.0
        if ( lname(i)(5:5).eq.'R') lstrong(i) = 2.0
        end if
        if ( lname(i)(1:1).eq.'O') lstrong(i) = -5.0
c        if ( red.and.lname(i)(1:1).eq.'O') lstrong(i) = 0.5
        end if
        end do

        nlines = i

        end if

        end if

        write(4,*) 'tol ',tol, www, nm
        write(6,*) 'tol ',tol, bth
        nmatch = 0
        lobt = -1
        if ( nterms.gt.3) then
        throbt = 0
        else
        throbt = bth
        end if
        do 1240 i=1,nm
        reject(i) = .false.

        if ( ht(i).ge.throbt) then

        refmax(i) = 0
        smin = 9999.9
        do 1230 k=10,lobt,-1

        do 1220 j=1,nlines
         if ( lstrong(j).ge.k.and..not.used(j)) then
         s = abs(wmax(i) - ldb(j))
         if ( s.lt.smin.and.(s.lt.fwhm16*disp0)) then
c     &        	(nr.le.9.and.s.lt.tol).or.
c     &        	(idiot.gt.3.and.nr.gt.9..and.s.lt.tol).or.
c     &        (idiot.gt.3.and.nr.gt.9.and.
c     &        	scat.lt.tol.and.s.lt.10*scat))) then
          if ( lstrong(j).ge.lstrong(refmax(i))) then
           used(refmax(i)) = .false.
           smin = s
           refmax(i) = j
           used(j) = .true.
        nmatch = nmatch + 1
          end if
         end if
        end if
 1220        continue

 1230        continue
        end if

 1240        continue
        write(6,*) 'tol ',tol, nmatch
        if ( numlfit.gt.20) goto 1020

        do i=1,nm
        write(4,212) i,lstrong(refmax(i)),
     &        temp1(i), wmax(i), ldb(refmax(i)), ht(i), lname(refmax(i))
        diff = abs(wmax(i) - ldb(refmax(i)))
        if ( temp1(i).gt.ndata-fwhm16) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        if ( temp1(i).lt.fwhm16) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        if ( diff.gt.tol) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        if ( fwhm(i).gt.1.*fwhm16.and.arc.and.nimages.eq.1) then
         used(refmax(i)) = .false.
         refmax(i) = 0
        end if
        end do
        write(4,*)

        nr = 0

        if (idiot.gt.2.and.scat.lt.wtol) then
         bth = max(bth-10.,5.)
        end if
        if (scat.lt.wtol.and.lobt.eq.-1.and.scat.eq.oldscat.and.
     &        bth.le.5.0.and.(nterms.eq.nbig.and.already)) goto 1020
c        if (idiot.gt.10.and.scat.le.wtol) goto 1020
        if (idiot.gt.20.and.(scat.gt.wtol.or.nterms.lt.nbig)) then
         write(6,*) 'We certainly do have a problem!', wtol
         err=2
         return
        end if

        if (nimages.eq.1.and.scat.gt.wtol
     &        .and.nterms.eq.nbig.and.nbig.gt.5) then
              nterms=max(nterms-1,2)
          already = .false.
        end if

        if ((nterms.le.nbig.and.already).or.scat.lt.wtol)
     &        goto 1010
        if (nterms.lt.nbig.and.idiot.lt.20) goto 1010
 1020        continue

        xx = (1.0 - bmid)/bmid0
        call fleg1(xx,pt,nterms)
        wavbegin = 0.0
        do j=1,nterms
        wavbegin = wavbegin + pt(j)*disco(j)
        end do

        xx = (real(ndata) - bmid)/bmid0
        call fleg1(xx,pt,nterms)
        wavend = 0.0
        do j=1,nterms
        wavend = wavend + pt(j)*disco(j)
        end do

        err = 0
        write(4,1289) scat, scat22
        write(6,1289) scat, scat22
 1289        format('RMS (final)= ',2f15.10,1x,'(weighted,unweighted)')

 1290   continue

        do j = 1,100
        if (j.le.50) then
         disco100(j) = disco(j)
        else
         disco100(j) = 0.0
        end if
        end do

        bmid = ndata/2.0 + 0.5
        bmid0= ndata/2.0

        do ix =1,ndata
         xx = (ix - bmid)/bmid0
         call fleg1(xx,pt,nterms)
         lambda= 0.0
         do j=1,nterms
         lambda= lambda+ pt(j)*disco(j)
         end do
         write(40,*) ix,lambda,raw0(ix)
        end do

        return
        end

