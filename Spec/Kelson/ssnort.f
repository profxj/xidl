        subroutine ssnort( nrows, ncols, oraw, pref, ibord, fwhma,
     &        xord, yord, verb, chi, ierr, coef,
     &        raw, sraw, srcor, y2, nh1, sfraw, sfcor)

        include "nmax.f"

        parameter( nxo=10, nyo=10)
        integer ron(nxo,nyo)

        integer nrows, ncols, ibord, ierr, xord, yord
        real pref, chi, fwhma
        logical verb

        real oraw(nrows*ncols)
        real*8 coef(100)
        real raw(ncols,nrows)

        integer rx1, rx2, orx1, orx2
        real xco(nmax), refer(nmax), sxco(nmax), tchi
        real thr, xcen2(16,rmax), y2(ncols,nrows), cen(16,rmax)
        integer imax(256)
        real zerol(nmax)
        real masterpt, cprev, refsend, xsend, look
        real xf(16,rmax), yf(16,rmax), cf(16,rmax), wt(16,rmax)
        integer ls(16,rmax)
        real bincen(16)
        logical skip(16), bad(16)

        real amask(2,100)
        integer nmask

        real xmask(2,100), xweight(nmax)
        integer nxmask
        real ymask(2,100), yweight(rmax)
        integer nymask

        real smask(2,100)
        integer nsmask

        real sraw(ncols,nrows), srcor(ncols,nrows), cref
        complex sfraw(nh1,nrows), sfcor(nh1,nrows)
        complex fref(nmax,16,2)

        integer ia(2), ib(2)
        integer row, nslits, row1, row2
        integer refrow, cenrow

        integer ix1m(16), ix2m(16), index(nmax)
        real com(16), s2n(16), wcol(nmax)

        real rimcen(16), srimcen(16)
        integer sbins(2),dbins(2)



        if (verb) write(6,*) "entered ssnort",nrows,ncols

        nxmask = 0
        nymask = 0

        chi = 0.0
        ierr = 0
        do i=1,nmax
        xco(i) = float(i)
        zerol(i) = 0.0
        end do

        if (verb) write(6,*) "here"
        do i=1,ncols
         do j=1,nrows
          raw(i,j) = oraw((j-1)*ncols+i)
         end do
        end do

        do i=1,ncols
        xweight(i) = 1.0
        end do

        do i=1,nrows
        yweight(i) = 1.0
        end do

        do j=1,nxmask
         do i=xmask(1,j),xmask(2,j)
          xweight(i) = 0.0
         end do
        end do

        do j=1,nymask
         do i=ymask(1,j),ymask(2,j)
          yweight(i) = 0.0
         end do
        end do

        rx1 = 0.0
        rx2 = 0.0

        islit = 1 
        nslits = 1

        if ( rx1.ge.rx2) then
        rx1 = 1
        rx2 = ncols
        end if

        ndata = (rx2-rx1)+1
        t2 = log10(real(ndata)-1)/log10(2.0)
        nt2 = int(t2)+1
        ndata = 2**(nt2)
        midrow = nrows/2

        nsex = 8
        s2nt = 2.0
        intl = 64
        intr = 16

        noskip = 0

        sbins(1) = 1
        sbins(2) = 1
        if ( nsex.gt.1) sbins(2) = nsex
        dbins(1) = 1
        dbins(2) = 1

        nimages = 1
        oldnslits = nslits

        row1 = 1
        row2 = nrows
        refrow = int(pref*nrows)
        cenrow = int(0.5*(row1+row2))

        maxisb = 2
        do isb=1,maxisb
        if (verb) write(6,*) "Working level",isb

        nsb = sbins(isb)
        subd = dbins(isb)
        irefch=0

  604        continue

        nnobs = 0
        do nob=1,int(subd)*nsb
        nsp = 1
        nbin = (rx2-rx1)+1
        if ( isb.gt.1) then
         nsp=2
         nbin = ((rx2-rx1)+1)/nsb
        end if

        ix1 = rx1 + (nob-1)*nbin/subd
        ix2 = ix1 + nbin/subd - 1

        nbin = ix2-ix1+1
        n2n = 2**(int(log10(float(nbin-1))/log10(2.0))+1)
        nhalf = n2n/2
        nhalf1 = nhalf + 1
        nbottom = max(refrow-8,row1+ibord)
        ntop = min(refrow+8,row2-ibord)

        do k=ix1,ix2
c        if (verb) write(6,*) "generating reference for",ix1,ix2

        refer(k-ix1+1) = 0.0
        nnn = 0

        do j=nbottom,ntop
        nnn = nnn + 1
        wcol(nnn) = raw(k,j)
        end do

        call quick( wcol, nnn, index)
        superw = 1.0

        if (nymask.eq.0) then
         superw = xweight(k)
        else
         if ( yweight(j).eq.0.and.xweight(k).eq.0) superw = 0.0
        end if

        refer(k-ix1+1) = wcol(nnn/2) * superw
        end do

        do k=nbin+1,n2n
        refer(k) = 0.0
        end do

        if ( nsp.gt.1) then

        thr = fwhma
        if (nrows.le.15) thr=20.0
        call cleanup( refer, n2n, nclean, thr,1)
        call thresh( refer, n2n, 1, n2n, s2n(nob), skip(nob), com(nob))

        ix1m(nob) = ix1
        ix2m(nob) = ix2
        skip(nob) = .false.
        com(nob) = com(nob) + ix1 - 1
        if ( s2n(nob).lt.s2nt) skip(nob)=.true.
        if ( com(nob).lt.2*fwhma.or.
     &        		com(nob).gt.ix2-2*fwhma) skip(nob)=.true.
        if (.not.skip(nob)) noskip = noskip + 1
        if (verb) write(6,'(i2,f6.1,f10.2,2i5,i3)')
     &                        nob,s2n(nob),com(nob),ix1,ix2,skip(nob)
        end if

        if ( .not.skip(nob)) nnobs = nnobs + 1

        call fft( n2n, refer, nhalf1, fref(1,nob,isb))

        end do

        if ( isb.eq.2.and.nnobs.lt.minobs) then
         s2nt = s2nt * 0.95
         goto 604
        end if

        end do

        chiold = 999.0
  990        continue
        look = 0.0
        diffold = 1e24

        istart = 1

        do 400 isb=1,maxisb
        if (verb) write(6,*) "Working level",isb

        nsb = sbins(isb)

        nsp = 1
        nbin = ncols

        if ( isb.gt.1) then
         nsp=2
         nbin = ((rx2-rx1)+1)/nsb
        end if
        n2n = 2**(int(log10(float(nbin-1))/log10(2.0))+1)
        nhalf = n2n / 2.
        nhalf1 = nhalf + 1.

        nob = 0

        ix1 = 1
        neff = 0

        subd = dbins(isb)

        xsendold = 0.0
        do 300 nob=1,int(subd)*nsb

        do j=row1,row2
        xcen2(nob,j) = 0.0
        end do

        if ( nsp.eq.1.or..not.skip(nob)) then

        ix1 = rx1 + (nob-1)*nbin
        ix2 = ix1 + nbin- 1
        if (verb) write(6,*) "Working interval",ix1,ix2,n2n,nhalf1

        neff = neff + 1

        if ( isb.eq.1) then
         bincen(nob) = 0.5*(ix1+ix2) - look
        else
         ix1 = ix1m(nob) - look
         ix2 = ix2m(nob) - look
         bincen(nob) = com(nob) - look
         bad(nob) = .false.
         binwid = (ix2-ix1)/8.
         if (bincen(nob).lt.8) bad(nob) = .true.
         if (bincen(nob).gt.ndata-8) bad(nob) = .true.
         if ( skip(nob)) bad(nob) = .true.
        end if

        if ( nsp.eq.1.or..not.bad(nob)) then


        if (verb) write(6,*) "Cleaning and Getting FFTs"
        do i=row1+ibord, row2-ibord
        l = i - row1 + 1
        do k=ix1,ix2
          sraw(k-ix1+1,l) = 0.0
          do kk=max(row1-i,-2),min(row2-i,2)
           if ( k.gt.0.and.k.le.ndata) then
            superw = 1.0
            if (nymask.eq.0) then
             superw = xweight(k)
            else
             if ( yweight(i+kk).eq.0.and.xweight(k).eq.0) superw = 0.0
            end if
            value = raw(k,i+kk) * superw
            sraw(k-ix1+1,l) = sraw(k-ix1+1,l) + value
           else if (k.lt.1) then
            sraw(k-ix1+1,l) = sraw(k-ix1+1,l) + 0.0
           else if (k.gt.ndata) then
            sraw(k-ix1+1,l) = sraw(k-ix1+1,l) + 0.0
           end if
          end do
          sxco(k-ix1+1) = float(k)
        end do

        do k=ix1,ix1+8
        sraw(k-ix1+1,l) = sraw(k-ix1+1,l) * real(k-ix1)/8.0
        end do
        do k=ix2-8,ix2
        sraw(k-ix1+1,l) = sraw(k-ix1+1,l) * real(ix2-k)/8.0
        end do
        do k=nbin+1,n2n
        sraw(k,l) = 0.0
        end do

        if ( nsp.gt.1) call cleanup(sraw(1,l), n2n, nclean, thr,1)
        call fft( n2n, sraw(1,l), nhalf1, sfraw(1,l))

        end do

        ia(1) = refrow
        ib(1) = row2-ibord
        ia(2) = row1+ibord
        ib(2) = refrow-1

        if (verb) write(6,*) "Getting cross-correlations and splining"
        do l=1,2
        do i=ia(l),ib(l)
        k = i - row1 + 1
        call correlate(nhalf1, fref(1,nob,isb), sfraw(1,k), sfcor(1,k))
        call invfft( nhalf1, sfcor(1,k), n2n, srcor(1,k))
        call shiftme( nsb, n2n, srcor(1,k))
        yp1 = srcor(1,k) - srcor(n2n,k)
        ypn = -yp1
        call spline(xco, srcor(1,k),n2n,yp1,ypn,y2(1,k))
        end do
        end do

        k = refrow - row1 + 1

        if ( isb.eq.1) then
        imax(k)=n2n/2
        do j=n2n/8,7*n2n/8
        if ( srcor(j,k).gt.srcor(imax(k),k))
     &        	imax(k) = j
        end do
        xsend = imax(k)
        else if ( isb.eq.2) then
          xsend= n2n/2
        do j=-intl,intl
        if ( srcor(j+n2n/2,k).gt.srcor(xsend,k))
     &        	xsend = j+n2n/2
        end do
         if ( nimages.ge.1) then
         refsend = xsend
         end if
        xsendold=xsend-n2n/2
        end if

        cen(nob,k)=getmax(n2n,xco,nint(xsend), srcor(1,k),y2(1,k),f)
        cref = cen(nob,k)
        cen(nob,k) = n2n/2. - cen(nob,k) - bincen(nob)
        cprev = cref
        masterpt = cref
        end if

        if ( isb.eq.1) look = cref - (1.0+n2n)/2.0

        ia(1) = refrow
        ib(1) = row2-ibord
        ia(2) = refrow-1
        ib(2) = row1+ibord

        if (verb) write(6,*) "Getting phase-shifts"
        do l=1,2
        cprev = cref
        iv = abs((ib(l)-ia(l)))/ (ib(l)-ia(l))

        do i=ia(l),ib(l),iv
        j = i - row1 + 1

        cen(nob,j)= getmax(n2n,xco,nint(cprev),srcor(1,j),y2(1,j),f)

        cprev = cen(nob,j)

        cen(nob,j) = n2n/2. - cen(nob,j)
        xcen2(nob,i) = bincen(nob) + cen(nob,j)

        if (i.eq.cenrow) then
           rimcen(nob) = xcen2(nob,i)
           srimcen(nob) = xcen2(nob,i)
        end if

        if ( isb.eq.2.and.i.eq.cenrow) then
         diff = look + (xcen2(nob,i) - rimcen(nob))
         if ( abs(diff).gt.intr) bad(nob)=.true.
        end if

        end do
        end do

        end if
  300        continue

  400        continue

        nbad = 0
c        write(20,*) look
        if ( maxisb.gt.1) then
        do nob=1,nsex
        diff = xcen2(nob,cenrow) - rimcen(nob)
        if (verb) write(6,302) nob,skip(nob),
     &                xcen2(nob,cenrow),rimcen(nob),diff,
     &          s2n(nob), bad(nob)
        if (skip(nob).or.bad(nob)) nbad = nbad + 1
  302        format(i4,i3,f10.3,f10.3,f10.3,f10.3,i3)	
        end do
        end if

        ngood = nsex - nbad

c        on(1,1) = 1
c        on(2,1) = 1
c        on(1,2) = 1
c        on(2,2) = 1
c        on(1,3) = 1
c        on(2,3) = 1
c        if ( ngood.gt.2) on(3,1) = 1
c        if ( ngood.gt.2) on(3,2) = 1
c        on(1,4) = 1
c        if ( ngood.gt.2) on(2,4) = 1
c        on(1,5) = 1
c        if ( ngood.gt.2) on(2,5) = 1
c        on(1,6) = 1

        ipf = 0
        do j=row1+ibord,row2-ibord

        inm = 1
        do nm=1,nsmask
        if ( j.ge.smask(1,nm).and.
     &        		j.le.smask(2,nm)) inm=0
        end do

        do nm=1,nmask
        if ( j.ge.amask(1,nm).and.
     &        		j.le.amask(2,nm)) inm=0
        end do

        rcenter = (rx2+rx1)/2.0
        rwidth = (rx2-rx1+1.0)/2.0
        if ( inm.eq.1) then
        ipf = ipf + 1
        do i=1,nsex
        xf(i,ipf) = (rimcen(i)-rcenter)/rwidth
        yf(i,ipf) = 2.0*(real(j)-0.5*(row2+row1))/(row2-row1)
        cf(i,ipf) = (xcen2(i,j)-rimcen(i))/ncols
        ls(i,ipf) = 1
        wt(i,ipf) = s2n(i)/real(nsex)
        end do
        end if

        end do

  800        continue

        if (verb) write(6,802) neff,ipf
  802        format('Gonna try to fit to ',i2,' and ',i5)

        do i=1,nxo
        do j=1,nyo
        ron(i,j) = 0
        if (i.le.xord+1.and.j.le.yord+1) ron(i,j) = 1
        end do
        end do

        if (verb) write(6,*) "going into reduce"
        call reduce( skip, bad, nsex, ipf, ron, islit, nslits,
     &        xf, yf, cf, wt, coef, chi, ierr, nrows, ncols)


        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine cleanup( array, nlen, nclean, th, bg)

        include "nmax.f"
        real array(nmax), draw(nmax), xs(nmax), sar(nmax)
        real narray(nmax), ndraw(nmax), nar(nmax)
        real fwhm(nmax), dum, dum2, dum3, th
        real larray(nmax)
        real tarray(nmax)
        integer nlen,nclean, index(nmax), bg

c        write(6,*) 'cleaning',th,nlen

        if (bg.eq.1) then
        do i=1,nlen
         do j=1,25
          k=min(max(i+j-13,1),nlen)
          larray(j) = array(k)
         end do
         call quick(larray,25,index)
         tarray(i) = larray(1)
        end do

        do i=1,nlen
         array(i) = array(i) - tarray(i)
        end do
        end if

        do i=1,nlen
        xs(i) = float(i)
        sar(i) = array(i)
        narray(i) = -array(i)
        nar(i) = narray(i)
        end do
c        write(6,*) (array(i),i=1,nlen)
c        write(6,*) (sar(i),i=1,nlen)

        yp1=0
        ypn=0
        call spline( xs, array, nlen, yp1, ypn, draw)
        yp1=0
        ypn=0
        call spline( xs, narray, nlen, yp1, ypn, ndraw)

        call quick(sar,nlen,index)
c        write(6,*) 'maxing'

        do i=nlen,nlen-3,-1
         dum3 = getlinemax( nlen, xs, index(i), array, draw,
     &        	dum2, fwhm(i), bottom)
c         if ( fwhm(i).lt.th/3.0) then
         if ( fwhm(i).lt.2.0) then
          dy = array(3*th+index(i)) - array(-3*th+index(i))
          dx = 6*th
          do j=-1,1
           k = j*th
           array(index(i)) = th*(j+1)*dy/dx + array(-3*th+index(i))
          end do
         else
          goto 10
         end if
        end do
  10        continue

        sumn = 0.0
        do i=1,nlen
        sumn = sumn + sar(i)
        if ( sumn.gt.0.0) goto 15
        end do
  15        continue
        ilo = i

        do i=ilo,ilo+2
        if ( index(i).gt.10.and.index(i).lt.nlen-10) then
         dum3 = getlinemax( nlen, xs, index(i), narray, ndraw,
     &        	dum2, fwhm(i), bottom)
         if ( fwhm(i).lt.th) then
          dy = array(3*th+index(i)) - array(-3*th+index(i))
          dx = 6*th
          do j=-1,1
           k = j*th
           array(index(i)) = th*(j+1)*dy/dx + array(-3*th+index(i))
          end do
         else
          goto 20
         end if
         end if
        end do
  20        continue

        return
        end

