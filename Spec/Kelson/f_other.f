cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real function profile( fwhm, j, x0)

	s = fwhm / 2.35
	x = real(j)
	pi = 3.14159265

	profile = exp( -0.5*((x-x0)/s)**2)/sqrt(2.0*pi*s**2)
	
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C#######################################################################
C
      SUBROUTINE  QUICK (DATUM, N, INDEX)
C
C=======================================================================
C
C A quick-sorting algorithm suggested by the discussion on pages 114-119
C of THE ART OF COMPUTER PROGRAMMING, Vol. 3, SORTING AND SEARCHING, by
C D.E. Knuth, which was referenced in Don Wells' subroutine QUIK.  This
C is my own attempt at encoding a quicksort-- PBS.
C
C Arguments
C
C DATUM (INPUT/OUTPUT) is a vector of dimension N containing randomly 
C        ordered real data upon input.  Upon output the elements of 
C        DATUM will be in order of increasing value.
C
C 
C INDEX (OUTPUT) is an integer vector of dimension N.  Upon return to
C       the calling program the i-th element of INDEX will tell where
C       the i-th element of the sorted vector DATUM had been BEFORE 
C       DATUM was sorted.
C
C=======================================================================
C
      IMPLICIT NONE
      INTEGER MAXSTK, N
      PARAMETER (MAXSTK=28)
C
C Parameter
C
C MAXSTK is the maximum number of entries the stack can contain.
C         A limiting stack length of 14 restricts this quicksort 
C         subroutine to vectors of maximum length of order 32,768 
C         (= 2**15).

      REAL DATUM(N)
      INTEGER INDEX(N), STKLO(MAXSTK), STKHI(MAXSTK)
C
      REAL DKEY
      INTEGER I, HI, LO, NSTAK, LIMLO, LIMHI, IKEY
C
C Initialize INDEX.
C
      DO I=1,N
         INDEX(I)=I
      END DO
C
C Initialize the pointers.
C
      NSTAK=0
      LIMLO=1
      LIMHI=N
C
  100 DKEY=DATUM(LIMLO)
      IKEY=INDEX(LIMLO)
C     TYPE *, 'LO =', LIMLO, '   HI =', LIMHI
C
C Compare all elements in the sub-vector between LIMLO and LIMHI with
C the current key datum.
C
      LO=LIMLO
      HI=LIMHI
  101 CONTINUE
C
      IF (LO .EQ. HI)GO TO 200
C
      IF (DATUM(HI) .LE. DKEY) GO TO 109
      HI=HI-1
C
C The pointer HI is to be left pointing at a datum SMALLER than the
C key, which is intended to be overwritten.
C
      GO TO 101
C
  109 DATUM(LO)=DATUM(HI)
      INDEX(LO)=INDEX(HI)
      LO=LO+1
  110 CONTINUE
C
      IF (LO .EQ. HI) GO TO 200
C
      IF (DATUM(LO) .GE. DKEY) GO TO 119
C
      LO=LO+1
      GO TO 110
C
  119 DATUM(HI)=DATUM(LO)
      INDEX(HI)=INDEX(LO)
      HI=HI-1
C
C The pointer LO is to be left pointing at a datum LARGER than the
C key, which is intended to be overwritten.
C
      GO TO 101
C
  200 CONTINUE
C
C LO and HI are equal, and point at a value which is intended to
C be overwritten.  Since all values below this point are less than
C the key and all values above this point are greater than the key,
C this is where we stick the key back into the vector.
C
      DATUM(LO)=DKEY
      INDEX(LO)=IKEY
C     DO 1666 I=LIMLO,LO-1
C1666 TYPE *, DATUM(I)
C     TYPE *, DATUM(LO), ' KEY'
C     DO 2666 I=LO+1,LIMHI
C2666 TYPE *, DATUM(I)
C
C At this point in the subroutine, all data between LIMLO and LO-1, 
C inclusive, are less than DATUM(LO), and all data between LO+1 and 
C LIMHI are larger than DATUM(LO).
C
C If both subarrays contain no more than one element, then take the most
C recent interval from the stack (if the stack is empty, we're done).
C If the larger of the two subarrays contains more than one element, and
C if the shorter subarray contains one or no elements, then forget the 
C shorter one and reduce the other subarray.  If the shorter subarray
C contains two or more elements, then place the larger subarray on the
C stack and process the subarray.
C
      IF (LIMHI-LO .GT. LO-LIMLO) GO TO 300
C
C Case 1:  the lower subarray is longer.  If it contains one or no 
C elements then take the most recent interval from the stack and go 
C back and operate on it.
C
      IF (LO-LIMLO .LE. 1) GO TO 400
C
C If the upper (shorter) subinterval contains one or no elements, then
C process the lower (longer) one, but if the upper subinterval contains
C more than one element, then place the lower (longer) subinterval on
C the stack and process the upper one.
C
      IF (LIMHI-LO .GE. 2) GO TO 250
C
C Case 1a:  the upper (shorter) subinterval contains no or one elements,
C so we go back and operate on the lower (longer) subinterval.
C
      LIMHI=LO-1
      GO TO 100
C
  250 CONTINUE
C
C Case 1b:  the upper (shorter) subinterval contains at least two 
C elements, so we place the lower (longer) subinterval on the stack and
C then go back and operate on the upper subinterval.
C 
      NSTAK=NSTAK+1
      IF (NSTAK .GT. MAXSTK) THEN
c        CALL STUPID ('Stack overflow in QUICK.  Increase MAXSTK.')
c        CALL OOPS
	stop
      END IF
      STKLO(NSTAK)=LIMLO
      STKHI(NSTAK)=LO-1
      LIMLO=LO+1
C     DO 3666 I=1,NSTAK
C3666 TYPE *, 'STACK: ', I, STKLO(I), STKHI(I)
      GO TO 100
C
  300 CONTINUE
C
C Case 2:  the upper subarray is longer.  If it contains one or no 
C elements then take the most recent interval from the stack and 
C operate on it.
C
      IF (LIMHI-LO .LE. 1) GO TO 400
C
C If the lower (shorter) subinterval contains one or no elements, then
C process the upper (longer) one, but if the lower subinterval contains
C more than one element, then place the upper (longer) subinterval on
C the stack and process the lower one.
C
      IF (LO-LIMLO .GE. 2) GO TO 350
C
C Case 2a:  the lower (shorter) subinterval contains no or one elements,
C so we go back and operate on the upper (longer) subinterval.
C
      LIMLO=LO+1
      GO TO 100
C
  350 CONTINUE
C
C Case 2b:  the lower (shorter) subinterval contains at least two 
C elements, so we place the upper (longer) subinterval on the stack and
C then go back and operate on the lower subinterval.
C 
      NSTAK=NSTAK+1
      IF (NSTAK .GT. MAXSTK) THEN
c        CALL STUPID ('Stack overflow in QUICK.  Increase MAXSTK.')
c        CALL OOPS
	stop
      END IF
      STKLO(NSTAK)=LO+1
      STKHI(NSTAK)=LIMHI
      LIMHI=LO-1
C     DO 4666 I=1,NSTAK
C4666 TYPE *, 'STACK: ', I, STKLO(I), STKHI(I)
      GO TO 100
C
  400 CONTINUE
C
C Take the most recent interval from the stack.  If the stack happens 
C to be empty, we are done.
C
      IF (NSTAK .LE. 0) THEN
         RETURN                           ! Normal return
      END IF
C     TYPE *, 'POP: ', NSTAK, STKLO(NSTAK), STKHI(NSTAK)
      LIMLO=STKLO(NSTAK)
      LIMHI=STKHI(NSTAK)
      NSTAK=NSTAK-1
      GO TO 100
C
      END!
C

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END

C#######################################################################
C
      SUBROUTINE RECTFY (X, NSTAR, INDEX, HOLD)
      IMPLICIT NONE
C
      REAL X(*), HOLD(*)
      INTEGER INDEX(*)
C
      INTEGER I, NSTAR
C
      DO I=1,NSTAR
         HOLD(I)=X(I)
      END DO
      DO I=1,NSTAR
         X(I)=HOLD(INDEX(I))
      END DO
      RETURN
      END!
C

C#######################################################################
C
      SUBROUTINE iRECTFY (X, NSTAR, INDEX, HOLD)
      IMPLICIT NONE
C
      integer X(*), HOLD(*)
      INTEGER INDEX(*)
C
      INTEGER I, NSTAR
C
      DO I=1,NSTAR
         HOLD(I)=X(I)
      END DO
      DO I=1,NSTAR
         X(I)=HOLD(INDEX(I))
      END DO
      RETURN
      END!
C

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real function getlinemax( npts, xco, imax, rcor, y2,
     &		f,fwhm, bot)

	include "nmax.f"
	real xco(npts), rcor(npts), y2(npts), temp(nmax)
	real d1, d2, d(20), f, fwhm, bot
	integer imax, npts, index(nmax)
	real del

	i0 = imax-30
	i1 = imax+30
	if ( i0.lt.1) i0 = 1
	if ( i1.gt.npts) i1 = npts
	do i=i0,i1
	temp(i-i0+1) = rcor(i)
	end do
	ni = i1-i0+1
	call quick( temp, ni, index)
	bot = temp(3)

	d(2) = xco(imax)

	del = 1.

	kmax = imax

c	do j=imax-int(del)/2,imax+int(del)/2
c	write(8,*) j,imax, rcor(j)
c	if ( rcor(j).gt.rcor(kmax)) kmax = j
c	end do
c
c	del = 1.
	d(2) = kmax

	do i=1,20

	d(1) = d(2) - del/2.
	d(3) = d(2) + del/2.
	call splint(xco, rcor, y2, npts, d(1), f1)
	call splint(xco, rcor, y2, npts, d(2), f2)
	call splint(xco, rcor, y2, npts, d(3), f3)

c	write(6,*) d(1), d(2), d(3)
c	write(6,*) f1,f2,f3
c	write(6,*) 

	if ( f1.gt.f2) d(2)=d(1)
	if ( f3.gt.f2) d(2)=d(3)

	del = del / 2.

	if ( del.lt.0.001) goto 20

	end do

   20	continue
c	write(6,*) 'center',d2
c	write(6,*) 'center',bot

	f = f2
	d2 = d(2)
c	write(6,*) 'center',d2, bot

	do x2=d2,d2-16.0,-0.005
	call splint(xco, rcor, y2, npts, x2, f2)
c	write(6,*) d2,f,x2,f2
	if ( (f2-bot).le.(f-bot)/2.0) goto 30
	end do
   30	continue
c	write(6,*) x2

	do x3=d2,d2+16.0,0.005
	call splint(xco, rcor, y2, npts, x3, f3)
c	write(6,*) d2,f,x3,f3
	if ( (f3-bot).le.(f-bot)/2.0) goto 40
	end do
   40	continue
c	write(6,*) x3 
	fwhm = x3-x2
c	write(6,*) fwhm, x2, x3

	do x2=d2,d2-16.0,-0.002
	call splint(xco, rcor, y2, npts, x2, f2)
c	write(6,*) d2,f,x2,f2
	if ( (f2-bot).le.0.6*(f-bot)) goto 130
	end do
  130	continue
c	write(6,*) x2

	do x3=d2,d2+16.0,0.002
	call splint(xco, rcor, y2, npts, x3, f3)
c	write(6,*) d2,f,x3,f3
	if ( (f3-bot).le.0.6*(f-bot)) goto 140
	end do
  140	continue
c	write(6,*) x3

	getlinemax = 0.5*(x2+x3)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine correlate(nhalf1, fref, fraw,fcor)

	include "nmax.f"
	integer nhalf1
	complex fref(nhalf1), fraw(nhalf1), fcor(nhalf1)

	do i=1,nhalf1
c	write(6,*) i,fref(i),fraw(i)
	 fcor(i) = fref(i) * conjg(fraw(i)) / (nhalf1-1.0)
	end do

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
c	include "nmax.f"
      PARAMETER (NMAX=4096)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real function getmax( npts, xco, imax, rcor, y2, f)

	include "nmax.f"
	real xco(npts), rcor(npts), y2(npts)
	real d1, d2, d(20)
	integer imax

	real del


	d(2) = xco(imax)

	del = 4.

	kmax = imax

c	do j=imax-int(del)/2,imax+int(del)/2
c	write(8,*) j,imax, rcor(j)
c	if ( rcor(j).gt.rcor(kmax)) kmax = j
c	end do

	del = 1.
c	d(2) = kmax

	do i=1,20

	d(1) = d(2) - del/2.
	d(3) = d(2) + del/2.

	call splint(xco, rcor, y2, npts, d(1), f1)
	call splint(xco, rcor, y2, npts, d(2), f2)
	call splint(xco, rcor, y2, npts, d(3), f3)

c	write(6,*) d(1), d(2), d(3)
c	write(6,*) f1,f2,f3
c	write(6,*) 

	if ( f1.gt.f2) d(2)=d(1)
	if ( f3.gt.f2) d(2)=d(3)

	del = del / 2.

	if ( del.lt.0.0005) goto 20

	end do

   20	continue

	f = f2

	getmax = d(2)


	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dfpoly(x,p,np)
	implicit real*8 (a-h,o-z)
      dimension p(np)
      p(1)=1.
       do 11 j=2,np
	 p(j) = p(j-1)*x
11      continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fpoly(x,p,np)
      dimension p(np)
      p(1)=1.
       do 11 j=2,np
	 p(j) = p(j-1)*x
11      continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fleg1(x,pl,nl)
      real pl(nl)
      pl(1)=1.
      pl(2)=x
      if(nl.gt.2) then
        twox=2.*x
        f2=x
        d=1.
        do 11 j=3,nl
          f1=d
          f2=f2+twox
          d=d+1.
          pl(j)=(f2*pl(j-1)-f1*pl(j-2))/d
11      continue
      endif
      return
      end

      subroutine fleg(x,pl,nl)
      real*8 pl(nl),twox,d,f2
      pl(1)=1.
      pl(2)=x
      if(nl.gt.2) then
        twox=2.*x
        f2=x
        d=1.
        do 11 j=3,nl
          f1=d
          f2=f2+twox
          d=d+1.
          pl(j)=(f2*pl(j-1)-f1*pl(j-2))/d
11      continue
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine shiftme( m, n, array)

	include "nmax.f"
	real array(n), temp(nmax)

	if ( m.gt.0) then

	 do i=1,n/2
	 j = i + n/2
	 temp(i) = array(j)
	 array(j) = array(i)
	 array(i) = temp(i)
	 end do

	else

	do j=1,m
	 do i=1,n
	 array(i+n) = array(i)
	 end do
	end do

	end if

c	do i=1,n
c	write(8,*) i,array(i)
c	end do

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fuser(x,p,np)

	real*8 p
      dimension p(np)

	call fleg(x,p,np)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	subroutine reduce( skip, bad, nxs, nyf, on, nref,
     &	   nslits, xf, yf, cf, wt, coef, chi, ierr, nrows, ncols)

	include "nmax.f"

	parameter( nxo=10, nyo=10)
	integer on(nxo,nyo), nrows, ncols

	integer xord
	integer yord, nslits

	logical skip(16), bad(16)

	real xf(16,rmax), yf(16,rmax), cf(16,rmax), wt(16,rmax)
	real*8 w(100,16*rmax), ww(16*rmax)
	real*8 coef(100), y(16*rmax), deter
	real*8 px(10), py(10)

c	write(6,*) 'in reduce!'

  400	continue

	in = 0
	nn = 0

	do j=1,nyf

	nwt = 1

	do i=1,nxs

c	if (abs(cf(i,j)-1024.5).gt.1008.0) goto 90

	call fleg( yf(i,j), py, nyo)
	call fleg( xf(i,j), px, nxo)

	if (.not.skip(i).and..not.bad(i)) then
	in = in + 1

	y(in) = cf(i,j)
	ww(in) = wt(i,j)

	is = 0
	is = is + 1
	w(is,in) = 1.0

	do iy=1,nyo
	do ix=1,nxo
	if ((ix.gt.1.or.iy.gt.1).and.on(ix,iy).gt.0) then
	is = is + 1
	w(is,in) = px(ix) * py(iy)
	end if
	end do
	end do

	nn = nn + 1

	end if

  90	continue

	end do

  95	continue
	end do

	nco = is

	nob = in

	if ( nco.gt.100) then
	write(6,*) 'Too many degrees of freedom',nco
	stop
	else if ( nob.eq.0) then
	write(6,*) 'Uh-oh',nob
	end if
	call kk( nco, nob, w, y, ww, coef, deter)

	chi = 0.0
	nsum = 0
	wsum = 0.0
	do j=1,nyf
	do i=1,nxs
	if (.not.skip(i).and..not.bad(i)) then
	 nsum = nsum + 1
	 yt = yf(i,j)
	 x = xf(i,j)
	 ctem = cf(i,j)*ncols
	 xtem = delta(on, nref, nref, nslits,
     &		 coef, x, yt)*ncols
	 chi = chi + (ctem - xtem)**2*wt(i,j)**2
	 wsum = wsum + wt(i,j)**2
	end if
  91	continue
	end do
	end do

	chi = sqrt(chi / wsum)

c	chi = sqrt(chi/real(nsum))
c	write(7,20) chi
  20	format('RMS = ',f15.9)

	ierr= 0

c	if (deter.eq.0.0) then
c	 ierr = 1
c	end if

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	SUBROUTINE KK(NDIMEN,NDATA,X,Y,W,A,deter)

C		SUBR. SOLVES  A FOR Y=SUMMATION(A*X)
C		WITH LEAST SSQARES METHOD

	implicit real*8 (a-h,o-z)
        parameter(nmax=100)
	INTEGER*4 NDATA,NDIMEN
	REAL*8  X(nmax,NDATA),Y(NDATA),A(nmax), W(NDATA)
	REAL*8  XXSUM(nmax,nmax),XYSUM(nmax),DET
        real*8 rmatinv


	LOGICAL	DEBUG/.TRUE./

	debug=.false.	

	DO I=1,nmax
	DO J=1,nmax

	XYSUM(J)=0.
	XXSUM(I,J)=0.

	END DO
	END DO

	SW1 = 0.0

	DO N=1,NDATA
	DO I=1,NDIMEN

c		XYSUM(I)=XYSUM(I)+X(I,N)*Y(N)
		XYSUM(I)=XYSUM(I)+X(I,N)*Y(N)*W(N)**2
		SW1 = SW1 + W(N)**2

		DO J=1,NDIMEN

c		  XXSUM(I,J)=XXSUM(I,J)+X(I,N)*X(J,N)
		  XXSUM(I,J)=XXSUM(I,J)+X(I,N)*X(J,N)*W(N)**2

		END DO

	END DO
	END DO

C		augment the matrix

	DO J=1,NDIMEN

c		XXSUM(NDIMEN+1,J)=XYSUM(J)/FLOAT(NDATA)
c		XXSUM(J,NDIMEN+1)=XYSUM(J)/FLOAT(NDATA)
		XXSUM(NDIMEN+1,J)=XYSUM(J)/FLOAT(NDATA)/SW1
		XXSUM(J,NDIMEN+1)=XYSUM(J)/FLOAT(NDATA)/SW1

		DO I=1,NDIMEN

c		  XXSUM(I,J)=XXSUM(I,J)/FLOAT(NDATA)
		  XXSUM(I,J)=XXSUM(I,J)/FLOAT(NDATA)/SW1

		END DO

	END DO


	IF (DEBUG) THEN
		DO J=1,NDIMEN
		WRITE(*,'(1X,11G)') (XXSUM(I,J),I=1,NDIMEN+1)
		END DO
	END IF

	EPS=1.d-40
	INDIC=0
	NARRAY=nmax
	
	DETER=RMATINV(NDIMEN,XXSUM,A,EPS,INDIC,NARRAY)

	IF (DEBUG) WRITE(*,'('' Determinant ='',G)') DETER

	IF (DEBUG) WRITE(*,'(/1X,11G//)') (A(I),I=1,NDIMEN)

	IF (DEBUG) THEN
		DO J=1,NDIMEN
		WRITE(*,'(1X,11G)') (XXSUM(I,J),I=1,NDIMEN)
		END DO
	END IF

	RETURN
	END


      real*8 FUNCTION RMATINV(N,A,X,EPS,INDIC,NRC)                               
C                                                                       
C             FUNCTION RMATINV
C                                                                       
C        THIS FUNCTION RETURNS THE VALUE OF THE DETERMINANT OF A        
C        MATRIX.  IN ADDITION, THE INVERSE MATRIX MAY BE CALCULATED     
C        IN PLACE, AND THE SOLUTION VECTOR OF THE CORRESPONDING LINEAR  
C        SYSTEM COMPUTED.  GAUSS-JORDAN ELIMINATION WITH  MAXIMUM       
C        PIVOT STRATEGY IS EMPLOYED, USING DOUBLE PRECISION ARITH-      
C        METIC.  IF THE MATRIX EXCEEDS THE MAXIMUM SIZE (50 BY 50),     
C        OR IF IT IS SINGULAR, A TRUE ZERO IS RETURNED.                 
C                                                                       
C        CALLING SEQUENCE:
C        RMATINV(N,A,X,EPS,INDIC,NRC)
C             N IS THE SIZE OF THE MATRIX (N BY N)
C             A IS THE MATRIX                                           
C             X IS THE SOLUTION VECTOR                                  
C             EPS IS A SMALL NUMBER TO BE USED AS A TEST FOR SINGULARITY
C             INDIC IS THE CONTROL PARAMETER:                           
C                  IF IT IS NEGATIVE, THE INVERSE IS COMPUTED IN PLACE  
C                  IF IT IS ZERO,THE MATRIX IS ASSUMED TO BE AUGMENTED, 
C                  AND THE SOLUTION AND INVERSE ARE COMPUTED            
C                  IF IT IS POSITIVE, ONLY THE SOLUTION IS COMPUTED     
C             NRC IS THE DIMENSION OF A IN THE CALLING PROGRAM          
C                                                                       
C        SUBPROGRAMS REQUIRED: NONE                                     
C                                                                       
	implicit real*8 (a-h,o-z)
      parameter(nmax=100)
      DIMENSION IROW(nmax),JCOL(nmax),JORD(nmax),Y(nmax)
      dimension A(NRC,NRC),X(N)
      MAX=N                                                             
      IF(INDIC.GE.0) MAX=N+1  
      IF(N.LE.nmax) GO TO 5                                               
      RMATINV=0.0
      RETURN                                                            
    5 DETER=1.0
      DO 18 K=1,N                                                       
      KM1=K-1
      PIVOT=0.0
      DO 11 I=1,N                                                       
      DO 11 J=1,N                                                       
      IF(K.EQ.1) GO TO 9                                                
      DO 8 ISCAN=1,KM1                                                  
      DO 8 JSCAN=1,KM1                                                  
      IF(I.EQ.IROW(ISCAN)) GO TO 11                                     
      IF(J.EQ.JCOL(JSCAN)) GO TO 11                                     
    8 CONTINUE                                                          
    9 IF( ABS(A(I,J)).LE. ABS(PIVOT)) GO TO 11                                  
      PIVOT=A(I,J)
      IROW(K)=I
      JCOL(K)=J
   11 CONTINUE                                                          
      IF( ABS(PIVOT).GT.EPS) GO TO 13                                           
      RMATINV=0.0
      RETURN                                                            
   13 IROWK=IROW(K)
      JCOLK=JCOL(K)
      DETER=DETER*PIVOT
      DO 14 J=1,MAX                                                     
   14 A(IROWK,J)=A(IROWK,J)/PIVOT
      A(IROWK,JCOLK)=1.0/PIVOT
      DO 18 I=1,N                                                       
      AIJCK=A(I,JCOLK)
      IF(I.EQ.IROWK) GO TO 18                                           
      A(I,JCOLK)=-AIJCK/PIVOT
      DO 17 J=1,MAX                                                     
   17 IF(J.NE.JCOLK) A(I,J)=A(I,J)-AIJCK*A(IROWK,J)
   18 CONTINUE                                                          
      DO 20 I=1,N                                                       
      IROWI=IROW(I)
      JCOLI=JCOL(I)
      JORD(IROWI)=JCOLI
   20 IF(INDIC.GE.0) X(JCOLI)=A(IROWI,MAX)
      INTCH=0
      NM1=N-1
      DO 22 I=1,NM1                                                     
      IP1=I+1
      DO 22 J=IP1,N                                                     
      IF(JORD(J).GE.JORD(I)) GO TO 22                                   
      JTEMP=JORD(J)
      JORD(J)=JORD(I)
      JORD(I)=JTEMP
      INTCH=INTCH+1
   22 CONTINUE                                                          
      IF(INTCH/2*2.LE.INTCH) DETER=-DETER
   24 IF(INDIC.LE.0) GO TO 26                                           
      RMATINV=DETER
      RETURN                                                            
   26 DO 28 J=1,N                                                       
      DO 27 I=1,N                                                       
      IROWI=IROW(I)
      JCOLI=JCOL(I)
   27 Y(JCOLI)=A(IROWI,J)
      DO 28 I=1,N                                                       
   28 A(I,J)=Y(I)
      DO 30 I=1,N                                                       
      DO 29 J=1,N                                                       
      IROWJ=IROW(J)
      JCOLJ=JCOL(J)
   29 Y(IROWJ)=A(I,JCOLJ)
      DO 30 J=1,N                                                       
   30 A(I,J)=Y(J)
      RMATINV=DETER
      RETURN                                                            
      END                                                               

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real function delta( on, is, rs, ns, coef, x, y)

	parameter( nxo=10, nyo=10)
	integer on(nxo,nyo), rs

	real*8 coef(100),del
	real*8 px(10), py(10)

c	write(6,*) (coef(k),k=1,10)
c	del = coef(is)
	del = coef(1)

	c = 1

	x2 = x
	y2 = y

	call fleg( x2, px, nxo)
	call fleg( y2, py, nyo)

	do iy=1,nyo
	do ix=1,nxo
	if ((ix.gt.1.or.iy.gt.1).and.on(ix,iy).gt.0) then
	c = c + 1
	del = del + coef(c) * px(ix) * py(iy)
	end if
	end do
	end do

	delta = del

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine thresh( array, n, i1, i2, h, skip, com)

	include "nmax.f"

	real array(n), com, array2(nmax), h
	integer index(nmax)
	logical skip

	sumx = 0.0

	do i=i1,i2
	sumx = sumx + array(i)
	array2(i) = array(i)
	end do

	call quick( array2, n, index)

	iii = 1
  41	continue
c	bot = array2(iii)
	bot = 0.5 * (array2(iii) + array2(iii+1))
	if ( bot.le.0.0) then
	 iii = iii + 1
	 goto 41
	end if

	xn = i2 - i1 + 1
	xmean = sumx / xn

	bot = 0.5 * (array2(iii+0.02*xn) + array2(iii+1+0.02*xn))

c	if ( xmean.lt.0.0) then
c	do i=i1,i2
c	write(14,*) i,array(i)
c	end do
c	end if
c	write(6,*) i1,i2,bot

	v2 = 0.0
	var = 0.0
	com = 0.0
	sumx2 = 0.0

	vn = 0.0
	vf = 0.0
	do i=i1,i2
	if ( array(i).ge.bot/2.0) then
	var = var + (array(i) - xmean)**2/(xn-1.)
	if (i.gt.1.and.i.lt.n) then
	com = com + (array(i+1)-array(i-1))**2*float(i-i1+1)
	v2 = v2 + (array(i+1)-array(i-1))**2
	vf = vf + (array(i+1)-array(i-1))**2*max((array(i)-bot),0.0)
	vn = vn + max((array(i)-bot),0.0)
	end if
	end if
	end do

c	write(6,*) v2,vn,var
	if ( v2.le.0.0.or.vn.le.0.0) then
	 write(6,*) 'v2=0'
	 stop
	end if
	com = com/v2 + i1 - 1

	v2 = sqrt(v2)/2./sqrt(xn-2.)

	sig = sqrt(var)

	skip = .false.

c	write(6,*) vf,vn,bot
c	h = v2/sqrt(bot)
	h = sqrt(vf/vn/bot)


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real*8 function dgauss(j,sig)

	implicit real*8 (a-h,o-z)

	pi = 3.14159265

	x = dble(j)

	dgauss = exp(-0.5*(x/sig)**2)/dsqrt(2.*pi*sig**2)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dpoly(x,p,np)
	implicit real*8 (a-h,o-z)
      dimension p(np)
      p(1)=1.
       do 11 j=2,np
	 p(j) = p(j-1)*x
11      continue
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LFIT(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,CHISQ,
     *FUNCS)
      PARAMETER (MMAX=50)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MA),
     *    COVAR(NCVM,NCVM),BETA(MMAX),AFUNC(MMAX)
      KK=MFIT+1
      DO 12 J=1,MA
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
          PAUSE 'Improper set in LISTA'
        ENDIF
12    CONTINUE
      IF (KK.NE.(MA+1)) PAUSE 'Improper set in LISTA'
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE
      DO 18 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-A(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE
      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
        A(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
      DO 24 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24    CONTINUE
      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      RETURN
      END


      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END

      subroutine gaussj(a,n,np,b,m,mp)
      parameter (nmax=50)
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix.'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      end


C
C#######################################################################
C
      SUBROUTINE  dQUICK (DATUM, N, INDEX)
C
C=======================================================================
C
C A quick-sorting algorithm suggested by the discussion on pages 114-119
C of THE ART OF COMPUTER PROGRAMMING, Vol. 3, SORTING AND SEARCHING, by
C D.E. Knuth, which was referenced in Don Wells' subroutine QUIK.  This
C is my own attempt at encoding a quicksort-- PBS.
C
C Arguments
C
C DATUM (INPUT/OUTPUT) is a vector of dimension N containing randomly 
C        ordered real data upon input.  Upon output the elements of 
C        DATUM will be in order of increasing value.
C
C 
C INDEX (OUTPUT) is an integer vector of dimension N.  Upon return to
C       the calling program the i-th element of INDEX will tell where
C       the i-th element of the sorted vector DATUM had been BEFORE 
C       DATUM was sorted.
C
C=======================================================================
C
      IMPLICIT NONE
      INTEGER MAXSTK, N
      PARAMETER (MAXSTK=28)
C
C Parameter
C
C MAXSTK is the maximum number of entries the stack can contain.
C         A limiting stack length of 14 restricts this quicksort 
C         subroutine to vectors of maximum length of order 32,768 
C         (= 2**15).

      REAL*8 DATUM(N)
      INTEGER INDEX(N), STKLO(MAXSTK), STKHI(MAXSTK)
C
      REAL*8 DKEY
      INTEGER I, HI, LO, NSTAK, LIMLO, LIMHI, IKEY
C
C Initialize INDEX.
C
      DO I=1,N
         INDEX(I)=I
      END DO
C
C Initialize the pointers.
C
      NSTAK=0
      LIMLO=1
      LIMHI=N
C
  100 DKEY=DATUM(LIMLO)
      IKEY=INDEX(LIMLO)
C     TYPE *, 'LO =', LIMLO, '   HI =', LIMHI
C
C Compare all elements in the sub-vector between LIMLO and LIMHI with
C the current key datum.
C
      LO=LIMLO
      HI=LIMHI
  101 CONTINUE
C
      IF (LO .EQ. HI)GO TO 200
C
      IF (DATUM(HI) .LE. DKEY) GO TO 109
      HI=HI-1
C
C The pointer HI is to be left pointing at a datum SMALLER than the
C key, which is intended to be overwritten.
C
      GO TO 101
C
  109 DATUM(LO)=DATUM(HI)
      INDEX(LO)=INDEX(HI)
      LO=LO+1
  110 CONTINUE
C
      IF (LO .EQ. HI) GO TO 200
C
      IF (DATUM(LO) .GE. DKEY) GO TO 119
C
      LO=LO+1
      GO TO 110
C
  119 DATUM(HI)=DATUM(LO)
      INDEX(HI)=INDEX(LO)
      HI=HI-1
C
C The pointer LO is to be left pointing at a datum LARGER than the
C key, which is intended to be overwritten.
C
      GO TO 101
C
  200 CONTINUE
C
C LO and HI are equal, and point at a value which is intended to
C be overwritten.  Since all values below this point are less than
C the key and all values above this point are greater than the key,
C this is where we stick the key back into the vector.
C
      DATUM(LO)=DKEY
      INDEX(LO)=IKEY
C     DO 1666 I=LIMLO,LO-1
C1666 TYPE *, DATUM(I)
C     TYPE *, DATUM(LO), ' KEY'
C     DO 2666 I=LO+1,LIMHI
C2666 TYPE *, DATUM(I)
C
C At this point in the subroutine, all data between LIMLO and LO-1, 
C inclusive, are less than DATUM(LO), and all data between LO+1 and 
C LIMHI are larger than DATUM(LO).
C
C If both subarrays contain no more than one element, then take the most
C recent interval from the stack (if the stack is empty, we're done).
C If the larger of the two subarrays contains more than one element, and
C if the shorter subarray contains one or no elements, then forget the 
C shorter one and reduce the other subarray.  If the shorter subarray
C contains two or more elements, then place the larger subarray on the
C stack and process the subarray.
C
      IF (LIMHI-LO .GT. LO-LIMLO) GO TO 300
C
C Case 1:  the lower subarray is longer.  If it contains one or no 
C elements then take the most recent interval from the stack and go 
C back and operate on it.
C
      IF (LO-LIMLO .LE. 1) GO TO 400
C
C If the upper (shorter) subinterval contains one or no elements, then
C process the lower (longer) one, but if the upper subinterval contains
C more than one element, then place the lower (longer) subinterval on
C the stack and process the upper one.
C
      IF (LIMHI-LO .GE. 2) GO TO 250
C
C Case 1a:  the upper (shorter) subinterval contains no or one elements,
C so we go back and operate on the lower (longer) subinterval.
C
      LIMHI=LO-1
      GO TO 100
C
  250 CONTINUE
C
C Case 1b:  the upper (shorter) subinterval contains at least two 
C elements, so we place the lower (longer) subinterval on the stack and
C then go back and operate on the upper subinterval.
C 
      NSTAK=NSTAK+1
      IF (NSTAK .GT. MAXSTK) THEN
c        CALL STUPID ('Stack overflow in QUICK.  Increase MAXSTK.')
c        CALL OOPS
	stop
      END IF
      STKLO(NSTAK)=LIMLO
      STKHI(NSTAK)=LO-1
      LIMLO=LO+1
C     DO 3666 I=1,NSTAK
C3666 TYPE *, 'STACK: ', I, STKLO(I), STKHI(I)
      GO TO 100
C
  300 CONTINUE
C
C Case 2:  the upper subarray is longer.  If it contains one or no 
C elements then take the most recent interval from the stack and 
C operate on it.
C
      IF (LIMHI-LO .LE. 1) GO TO 400
C
C If the lower (shorter) subinterval contains one or no elements, then
C process the upper (longer) one, but if the lower subinterval contains
C more than one element, then place the upper (longer) subinterval on
C the stack and process the lower one.
C
      IF (LO-LIMLO .GE. 2) GO TO 350
C
C Case 2a:  the lower (shorter) subinterval contains no or one elements,
C so we go back and operate on the upper (longer) subinterval.
C
      LIMLO=LO+1
      GO TO 100
C
  350 CONTINUE
C
C Case 2b:  the lower (shorter) subinterval contains at least two 
C elements, so we place the upper (longer) subinterval on the stack and
C then go back and operate on the lower subinterval.
C 
      NSTAK=NSTAK+1
      IF (NSTAK .GT. MAXSTK) THEN
c        CALL STUPID ('Stack overflow in QUICK.  Increase MAXSTK.')
c        CALL OOPS
	stop
      END IF
      STKLO(NSTAK)=LO+1
      STKHI(NSTAK)=LIMHI
      LIMHI=LO-1
C     DO 4666 I=1,NSTAK
C4666 TYPE *, 'STACK: ', I, STKLO(I), STKHI(I)
      GO TO 100
C
  400 CONTINUE
C
C Take the most recent interval from the stack.  If the stack happens 
C to be empty, we are done.
C
      IF (NSTAK .LE. 0) THEN
         RETURN                           ! Normal return
      END IF
C     TYPE *, 'POP: ', NSTAK, STKLO(NSTAK), STKHI(NSTAK)
      LIMLO=STKLO(NSTAK)
      LIMHI=STKHI(NSTAK)
      NSTAK=NSTAK-1
      GO TO 100
C
      END!
C
      SUBROUTINE REALFT(DATA,N,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
        DATA(2*N+1)=DATA(1)
        DATA(2*N+2)=DATA(2)
      ELSE
        C2=0.5
        THETA=-THETA
        DATA(2*N+1)=DATA(2)
        DATA(2*N+2)=0.0
        DATA(2)=0.0
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      N2P3=2*N+3
      DO 11 I=1,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        DATA(2)=DATA(2*N+1)
      ELSE
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END

      subroutine fft(nin,in,nout,out)
 
      real*4 in(nin)
      real*4 h1(32768)
      real*8 dble
	integer nin, nout
      complex out(nout),dcmplx,c
 
c		subroutine is interface between program
c		and NRECIPES fast fourier transform routine
 
      rnin=sqrt(float(nin))
      nhalf=nin/2
      do 10 i=1,nin
      h1(i)=(in(i))/rnin
10    continue
 
      nin2=nin/2

      isign=1
      call realft(h1,nin2,isign)

 
      do i=1,nout
	if (i.eq.1) then
	      out(i)=cmplx( (h1(i)), 0.0)
	end if
	if (i.gt.1.and.i.le.nhalf) then
	      out(i)=cmplx( (h1(2*i-1)), (h1(2*i)) )
	end if
	if (i.eq.nhalf+1) then
	      out(i)=cmplx( (h1(2)), 0.0)
	end if
	if (i.gt.nhalf+1) then
	      out(i)=0.0
	end if
      end do
 
      return
      end 

	subroutine invfft(nin,in,nout,out)


	complex in(nin)
	real*4  out(nout)
	real*4  h1(32768)
	real*8  dimag,dreal

	if (nout.ne.2*(nin-1)) then
		write(32,'(''INVFFT: Dimensions wrong '',2i)') nin,nout
		stop
	end if

	fac=sqrt(float(nout))

	do i=1,nin-1
c	h1(2*i-1)=dreal(in(i))/fac
	h1(2*i-1)=real(in(i))/fac
	end do

c	h1(2)=dreal(in(nin))/fac
	h1(2)=real(in(nin))/fac
	do i=2,nin-1
c	h1(2*i)=dimag(in(i))/fac
	h1(2*i)=imag(in(i))/fac
	end do

	isign=-1

	nout2=nout/2

	call realft(h1,nout2,isign)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C A factor 2.0 was added in the next expression. This factor seems to
C have been forgotten in an earlier version of this routine. The error
C is probably due to the difference in definitions used by Marijn and
C Numerical Recipes. The old version had been correct if the
C last line of the comment on the routine REALFT (Num. Rec. page 400)
C had read :
C   " Result in this case must be divided by 1/(2N)"
C However, it reads
C   " Result in this case must be divided by 1/N"
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	do i=1,nout
	  out(i) = 2.0 * h1(i)
	end do

	return
	end 

