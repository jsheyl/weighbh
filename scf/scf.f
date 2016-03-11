

C***********************************************************************
C
C
                             PROGRAM scf
C
C
C***********************************************************************
C
C
C    A code to evolve self-gravitating systems using a self-consistent
C    field approach.  This version has been optimized for supercomputers
C    and is fully vectorized.  The code is written in standard FORTRAN,
C    although some CRAY-specific vector intrinsic routines have been
C    used.
C
C    The computational system of units is determined by the input data.
C    No explicit assumtions about the value of the gravitational 
C    constant have been made; it is read in as a parameter.
C    Particles are not required to have identical masses.
C
C
C                       Version 1: January 1, 1991
C
C
C                    Lars Hernquist, U.C. Santa Cruz
C
c
C
C=======================================================================
C
C
C     This is the top-level evolution program scf.  Its tasks are:
C
C          1) to initialize file structures and global variables;
C          2) to input parameters and the initial system state;
C          3) to advance the state of the system for a given number
C             of timesteps;
C          4) to perform a diagnostic analysis of the system at
C             each time step (energy, angular momentum, etc.);
C          5) to periodically record the state of the system;
C          6) and to terminate the simulation and close data files.
C
C
C=======================================================================
C
C
C     Basic global variables/parameters:
C
C          ax,ay,az    : accelerations of bodies.
C          clm, dlm,   : radial functions used to evaluate expansions.
C          elm, flm
C          costh       : cosine of polar angular coordinate of bodies.
C          cputime     : cpu time (secs) used during the simulation.
C          cputime0    : cumulative cpu time at start of run.
C          cputime1    : cumulative cpu time at end of run.
C          dplm        : values of derivatives of Legendre functions 
C                        for bodies.
C          dtime       : the timestep.
C          fixacc      : option to force conservation of linear
C                        momentum by setting acceleration of c.o.m.=0.
C          G           : the gravitational constant.
C          headline    : identification string for the run.
C          inptcoef    : option to read in expansion coefficients.
C          irsort      : vector of sorting indices.
C          lmax        : number of angular eigenfunctions.
C          mass        : masses of bodies.
C          nbodies     : total number of bodies.
C          nbodsmax    : maximum number of bodies.
C          nmax        : number of radial eigenfunctions.
C          noutbod     : frequency of system record outputs.
C          noutlog     : frequency of outputs to log file.
C          nsteps      : number of time-steps to integrate the system.
C          one         : the constant 1.
C          onesixth    : the constant 1/6.
C          outpcoef    : option to write out expansion coefficients.
C          phi         : azimuthal angular coordinate of bodies.
C          plm         : values of Legendre functions for bodies.
C          pi          : the constant pi.
C          pot         : potentials of bodies (self-gravity).
C          potext      : potentials of bodies (external field).
C          r           : radial coordinate of each body.
C          selfgrav    : option to turn off (.FALSE.) system self-
C                        gravity.
C          temp1,temp2,: temporary arrays used to improve the
C          temp3,temp4,  efficiency of the calculation of 
C          temp5,temp6   self-gravity.
C          tnow        : current system time.
C          tpos        : current position time.
C          tvel        : current velocity time.
C          twoopi      : the constant 2./pi.
C          ultrass,    : vectors containing values of ultraspherical
C          ultrasp1      polynomials for the bodies.
C          vx,vy,vz    : velocities of bodies.
C          x,y,z       : positions of bodies.
C          xi          : radially transformed coordinate of bodies.
C          zeroeven    : option to zero out all even terms in the
C                        basis function expansions.
C          zeroodd     : option to zero out all odd terms in the
C                        basis function expansions.
C
c
cXXX
c
c	indexm 		:array to order masses
c	iistar		:array carries stellar type flags to identify evolved stars
c	nmgrpmax	:maximum number of mass groups
c	nmgrp		:variable input number of mass groups < nmgrpmax
c	tlim		:array of stellar evolutionary timescales
c	masslim		:boundaries of mass groups
c	vbar		:array of mean velocities for delta v
c	vsig		:array of standard dev for delta v
c	ximf		:salpeter index
c	sigma		:cluster dispersion is SI
c	r0		:cluster core radius in pc
c	t0		:time scale for evolution
c	tscale		:r0/sigma
c
cXXX
c
C
C-----------------------------------------------------------------------
C
C   Definitions specific to input/output.
C
C          uterm, upars, ulog, ubodsin,   : logical i/o unit numbers.
C            ubodsout,utermfil,uoutcoef,
C            uincoef,ubodsel
C          parsfile, logfile, ibodfile,   : character names of files.
C            obodfile,termfile,outcfile,
C            incfile,elfile
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

C   Initialize state of the system.
C   -------------------------------
        CALL initsys
C            -------

C   Advance system state for a given number of steps.
C   -------------------------------------------------

        DO 100 n=1,nsteps
c

           CALL stepsys(n)
C               -------

 100    CONTINUE

C   Terminate the simulation.
C   -------------------------
        CALL endrun
C            ------

        STOP
        END
C***********************************************************************
C
C
        SUBROUTINE accp_LH
C
C
C***********************************************************************
C
C
C     Subroutine to compute accelerations, potential, and density.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        INTEGER k,l,m,n,lmin,lskip
        LOGICAL firstc
        REAL*8 anltilde,knl,sinth,smphi,cmphi,phinltil,deltam0,factrl,
     &         gammln,arggam,argult,argult1,sinsum,cossum,coeflm,
     &         dblfact,ttemp5,ar,ath,aphi,axt,ayt,azt

        DIMENSION argult(0:lmax),argult1(0:lmax),
     &            anltilde(0:nmax,0:lmax),knl(0:nmax,0:lmax),
     &            sinsum(0:nmax,0:lmax,0:lmax),dblfact(lmax+1),
     &            cossum(0:nmax,0:lmax,0:lmax),coeflm(0:lmax,0:lmax),
     &            ar(nbodsmax),ath(nbodsmax),aphi(nbodsmax)

        EQUIVALENCE (ar(1),ax(1)),(ath(1),ay(1)),(aphi(1),az(1))

        DATA firstc/.TRUE./

        SAVE firstc,dblfact,argult,argult1,anltilde,knl,coeflm,lmin,
     &       lskip

C=======================================================================

        IF(firstc) THEN

           firstc=.FALSE.

           dblfact(1)=1.

           DO 5 l=2,lmax
              dblfact(l)=dblfact(l-1)*(2.*l-1.)
 5         CONTINUE

           DO 20 n=0,nmax
              DO 10 l=0,lmax
                 knl(n,l)=0.5*n*(n+4.*l+3.)+(l+1.)*(2.*l+1.)
                 anltilde(n,l)=-2.**(8.*l+6.)*FACTRL(n)*(n+2.*l+1.5)
                 arggam=2.*l+1.5
                 anltilde(n,l)=anltilde(n,l)*(EXP(GAMMLN(arggam)))**2
                 anltilde(n,l)=anltilde(n,l)/(4.*pi*knl(n,l)*FACTRL(n+
     &                         4*l+2))
 10           CONTINUE
 20        CONTINUE

           DO 25 l=0,lmax

              argult(l)=2.*l+1.5
              argult1(l)=2.*l+2.5

              DO 23 m=0,l
                 deltam0=2.
                 IF(m.EQ.0) deltam0=1.
                 coeflm(l,m)=(2.*l+1.)*deltam0*FACTRL(l-m)/FACTRL(l+m)
 23           CONTINUE
 25        CONTINUE

           lskip=1
           IF(zeroodd.OR.zeroeven) lskip=2

           lmin=0
           IF(zeroeven) lmin=1

        ENDIF

!$OMP PARALLEL DO
        DO 30 k=1,nbodies
           r(k)=SQRT(x(k)**2+y(k)**2+z(k)**2)
           phi(k)=ATAN2(y(k),x(k))
           xi(k)=(r(k)-1.)/(r(k)+1.)
           costh(k)=z(k)/r(k)
           pot(k)=0.0
           ar(k)=0.0
           ath(k)=0.0
           aphi(k)=0.0
 30     CONTINUE
!$OMP END PARALLEL DO
        
        DO 60 l=0,lmax
           DO 50 m=0,l
              DO 40 n=0,nmax
                 sinsum(n,l,m)=0.0
                 cossum(n,l,m)=0.0
 40           CONTINUE
 50        CONTINUE
 60     CONTINUE

        DO 118 l=lmin,lmax,lskip
        
!$OMP PARALLEL DO
           DO 110 k=1,nbodies
              temp5(k)=r(k)**l/((1.+r(k))**(2*l+1))*mass(k)
 110       CONTINUE
!$OMP END PARALLEL DO
           DO 116 m=0,l

              CALL plgndrv(nbodies,l,m,dblfact,costh,temp1,temp2,plm)
C                  -------

!$OMP PARALLEL DO
              DO 111 k=1,nbodies
                 ttemp5=temp5(k)*plm(k)*coeflm(l,m)
                 temp3(k)=ttemp5*SIN(m*phi(k))
                 temp4(k)=ttemp5*COS(m*phi(k))
 111          CONTINUE
!$OMP END PARALLEL DO

              DO 114 n=0,nmax

                 CALL ultrasv(nbodies,n,argult(l),xi,temp1,temp2,
C                     -------
     &                        ultrasp)

!$OMP PARALLEL DO
                 DO 112 k=1,nbodies
                    sinsum(n,l,m)=sinsum(n,l,m)+temp3(k)*ultrasp(k)
                    cossum(n,l,m)=cossum(n,l,m)+temp4(k)*ultrasp(k)
 112             CONTINUE
!$OMP END PARALLEL DO

                 sinsum(n,l,m)=sinsum(n,l,m)*anltilde(n,l)
                 cossum(n,l,m)=cossum(n,l,m)*anltilde(n,l)

 114          CONTINUE
 116       CONTINUE
 118    CONTINUE

        IF(inptcoef.OR.outpcoef) CALL iocoef(sinsum,cossum)
C                                     ------

        DO 190 l=lmin,lmax,lskip

!$OMP PARALLEL DO
           DO 126 k=1,nbodies
              temp3(k)=0.0
              temp4(k)=0.0
              temp5(k)=0.0
              temp6(k)=0.0
 126       CONTINUE
!$OMP END PARALLEL DO

           DO 180 m=0,l

              CALL dplgndrv(nbodies,l,m,dblfact,costh,temp1,temp2,plm,
C                  --------
     &                      dplm)

!$OMP PARALLEL DO
              DO 130 k=1,nbodies
                 clm(k)=0.0
                 dlm(k)=0.0
                 elm(k)=0.0
                 flm(k)=0.0
 130          CONTINUE
!$OMP END PARALLEL DO

              DO 150 n=0,nmax

                 CALL ultrasv1(nbodies,n,argult(l),xi,temp1,temp2,
C                     --------
     &                         ultrasp,ultrasp1)

!$OMP PARALLEL DO
                 DO 140 k=1,nbodies
                    clm(k)=clm(k)+ultrasp(k)*cossum(n,l,m)
                    dlm(k)=dlm(k)+ultrasp(k)*sinsum(n,l,m)
                    elm(k)=elm(k)+ultrasp1(k)*cossum(n,l,m)
                    flm(k)=flm(k)+ultrasp1(k)*sinsum(n,l,m)
 140             CONTINUE
!$OMP END PARALLEL DO
 150          CONTINUE

!$OMP PARALLEL DO
              DO 170 k=1,nbodies
                 cmphi=COS(m*phi(k))
                 smphi=SIN(m*phi(k))
                 temp3(k)=temp3(k)+plm(k)*(clm(k)*cmphi+dlm(k)*smphi)
                 temp4(k)=temp4(k)-plm(k)*(elm(k)*cmphi+flm(k)*smphi)
                 temp5(k)=temp5(k)-dplm(k)*(clm(k)*cmphi+dlm(k)*smphi)
                 temp6(k)=temp6(k)-m*plm(k)*(dlm(k)*cmphi-clm(k)*smphi)
 170          CONTINUE
!$OMP END PARALLEL DO
 180       CONTINUE

!$OMP PARALLEL DO
           DO 183 k=1,nbodies
              phinltil=r(k)**l/((1.+r(k))**(2*l+1))
              pot(k)=pot(k)+temp3(k)*phinltil
              ar(k)=ar(k)+phinltil*(-temp3(k)*(l/r(k)-(2.*l+1.)/
     &              (1.+r(k)))+temp4(k)*4.*(2.*l+1.5)/(1.+r(k))**2)
              ath(k)=ath(k)+temp5(k)*phinltil
              aphi(k)=aphi(k)+temp6(k)*phinltil
 183       CONTINUE
!$OMP END PARALLEL DO

 190    CONTINUE

!$OMP PARALLEL DO
        DO 200 k=1,nbodies
           sinth=SQRT(1.-costh(k)**2)
           pot(k)=pot(k)*G
           ar(k)=G*ar(k)
           ath(k)= -G*sinth*ath(k)/r(k)
           aphi(k)=G*aphi(k)/(r(k)*sinth)
           cmphi=COS(phi(k))
           smphi=SIN(phi(k))
           axt=sinth*cmphi*ar(k)+costh(k)*cmphi*ath(k)-smphi*aphi(k)
           ayt=sinth*smphi*ar(k)+costh(k)*smphi*ath(k)+cmphi*aphi(k)
           azt=costh(k)*ar(k)-sinth*ath(k)
           ax(k)=axt
           ay(k)=ayt
           az(k)=azt
 200    CONTINUE
!$OMP END PARALLEL DO

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE accp_LHa
C
C
C***********************************************************************
C
C
C     Subroutine to compute accelerations, potential, and density.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        INTEGER k

C=======================================================================

!$OMP PARALLEL DO
        DO 10 k=1,nbodies
           r(k)=SQRT(x(k)**2+y(k)**2+z(k)**2)
           ax(k)=ax(k)-G*x(k)/(r(k)*(1.+r(k))**2)
           ay(k)=ay(k)-G*y(k)/(r(k)*(1.+r(k))**2)
           az(k)=az(k)-G*z(k)/(r(k)*(1.+r(k))**2)
           potext(k)=potext(k)-G/(1.+r(k))
 10     CONTINUE
!$OMP END PARALLEL DO

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE accpot
C
C
C***********************************************************************
C
C
C     Subroutine to compute accelerations, potential, and density.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        INTEGER i

C=======================================================================

        IF(selfgrav) THEN
           CALL accp_LH
C               -------
        ENDIF

!$OMP PARALLEL DO
        DO 5 i=1,nbodies
           potext(i)=0.0
 5      CONTINUE
!$OMP END PARALLEL DO

        IF(.NOT.selfgrav) THEN
!$OMP PARALLEL DO
           DO 10 i=1,nbodies
              ax(i)=0.0
              ay(i)=0.0
              az(i)=0.0
              pot(i)=0.0
 10        CONTINUE
!$OMP END PARALLEL DO

           CALL accp_LHa
C               --------
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE checkinp
C
C
C***********************************************************************
C
C
C     Subroutine to check consistency of input parameters and data,
C     output warnings to the terminal and/or log file, and terminate
C     the simulation if necessary.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C=======================================================================

        IF(nsteps.LT.0.OR.nsteps.GT.10000000)
     &     CALL terror(' input error for parameter nsteps ')
C               ------

        IF(noutbod.LT.0)
     &     CALL terror(' input error for parameter noutbod ')
C               ------

        IF(noutlog.LT.0)
     &     CALL terror(' input error for parameter noutlog ')
C               ------

        IF(dtime.LE.-1.e20.OR.dtime.GT.1.e20)
     &     CALL terror(' input error for parameter dtime ')
C               ------

        IF(G.LE.0.0.OR.G.GT.1.e20)
     &     CALL terror(' input error for parameter G ')
C               ------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE corracc
C
C
C***********************************************************************
C
C
C     Subroutine to correct accelerations so that the center of
C     mass remains fixed at the origin.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i
        REAL*8 axcm,aycm,azcm,mtot

C=======================================================================

        axcm=0.0
        aycm=0.0
        azcm=0.0
        mtot=0.0

!$OMP PARALLEL DO
        DO 10 i=1,nbodies
           mtot=mtot+mass(i)
           axcm=axcm+mass(i)*ax(i)
           aycm=aycm+mass(i)*ay(i)
           azcm=azcm+mass(i)*az(i)
 10     CONTINUE
!$OMP END PARALLEL DO

        axcm=axcm/mtot
        aycm=aycm/mtot
        azcm=azcm/mtot

!$OMP PARALLEL DO
        DO 20 i=1,nbodies
           ax(i)=ax(i)-axcm
           ay(i)=ay(i)-aycm
           az(i)=az(i)-azcm
 20     CONTINUE
!$OMP END PARALLEL DO

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE corrvel(rc)
C
C
C***********************************************************************
C
C
C     Subroutine to synchronize particle coordinates when outputing
C     particle data to body data file or when computing diagnostics.
C  
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 rc
        INTEGER p
        REAL*8 rcsign

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        IF(rc.EQ.'correct') THEN
           rcsign=-1.
        ELSE
           rcsign=1.
        ENDIF

!$OMP PARALLEL DO
        DO 10 p=1,nbodies
           vx(p)=vx(p)+rcsign*ax(p)*0.5*dtime
           vy(p)=vy(p)+rcsign*ay(p)*0.5*dtime
           vz(p)=vz(p)+rcsign*az(p)*0.5*dtime
 10     CONTINUE
!$OMP END PARALLEL DO

C   Update velocity time, system time.
C   ----------------------------------
        tvel=tvel+rcsign*0.5*dtime
        tnow=tvel

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE dplgndrv(nbodies,l,m,dblfact,x,plm1m,temp2,plm,dplm)
C
C
C***********************************************************************
C
C
C     A routine to compute derivatives of Legendre functions.
C
C
C=======================================================================

        INTEGER l,m,k,nbodies
        REAL*8 x(1),dplm(1),plm(1),plm1m(1),temp2(1),dblfact(1)

        CALL plgndrv(nbodies,l,m,dblfact,x,plm1m,temp2,plm)
C            -------

        IF(l.EQ.0) THEN
!$OMP PARALLEL DO
           DO 10 k=1,nbodies
              dplm(k)=0.0
 10        CONTINUE
!$OMP END PARALLEL DO
        ELSE
           IF(l.NE.m) THEN
              CALL plgndrv(nbodies,l-1,m,dblfact,x,dplm,temp2,plm1m)
C                  -------
           ELSE
!     $omp parallel do
              DO 100 k=1,nbodies
                 plm1m(k)=0.0
 100          CONTINUE
!     $omp end parallel do
           ENDIF

!     $omp parallel do
           DO 110 k=1,nbodies
              IF(x(k).NE.1.0.and.x(k).NE.-1.0) THEN
                 dplm(k)=(l*x(k)*plm(k)-(l+m)*plm1m(k))/
     &                   (x(k)*x(k)-1.)
              ELSE
                 IF(m.GT.2) THEN
                    dplm(k)=0.0
                 ENDIF
                 
                 IF(m.EQ.0) THEN
                    dplm(k)=0.5*l*(l+1.)
                    IF(MOD(l,2).EQ.0.AND.x(k).EQ.-1.0) dplm(k)=-dplm(k)
                 ENDIF

                 IF(m.EQ.1) THEN
                    dplm(k)=0.5*l*(l+1.)
                    IF(MOD(l,2).EQ.1.AND.x(k).EQ.-1.0) dplm(k)=-dplm(k)
                    dplm(k)=dplm(k)*1.e10
                 ENDIF

                 IF(m.EQ.2) THEN
                    dplm(k)=-0.25*(l-1.)*l*(l+1.)*(l+2.)
                    IF(MOD(l,2).EQ.0.AND.x(k).EQ.-1.0) dplm(k)=-dplm(k)
                 ENDIF
              ENDIF
 110       CONTINUE
!     $omp end parallel do

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE endrun
C
C
C***********************************************************************
C
C
C     Subroutine to end the simulation.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        REAL second

C=======================================================================

        cputime1=SECOND()

        CALL stopout
C            -------

        RETURN
        END
C***********************************************************************
C
C
        FUNCTION FACTRL(N)
C
C
C***********************************************************************
C
C
C     A function to compute factorials.  (From numerical recipes.)
C
C
C=======================================================================

        INTEGER n,ntop,j
        REAL*8 factrl,a,gammln,arggam

        DIMENSION A(33)

        DATA NTOP,A(1)/0,1./

        IF (N.LT.0) THEN
          PAUSE 'negative factorial'
        ELSE IF (N.LE.NTOP) THEN
          FACTRL=A(N+1)
        ELSE IF (N.LE.32) THEN
          DO 11 J=NTOP+1,N
            A(J+1)=J*A(J)
11        CONTINUE
          NTOP=N
          FACTRL=A(N+1)
        ELSE
          arggam=n+1.
          FACTRL=EXP(GAMMLN(arggam))
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        FUNCTION GAMMLN(XX)
C
C
C***********************************************************************
C
C
C     A routine to compute the natural logarithm of the gamma
C     function.  (Taken from numerical recipes.)
C
C
C=======================================================================

        INTEGER j

        REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER,gammln,xx

        DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     &      -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
        DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/

        X=XX-ONE
        TMP=X+FPF
        TMP=(X+HALF)*LOG(TMP)-TMP
        SER=ONE

        DO 11 J=1,6
          X=X+ONE
          SER=SER+COF(J)/X
11      CONTINUE

        GAMMLN=TMP+LOG(STP*SER)

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE inbods
C
C
C***********************************************************************
C
C
C     Subroutine to read phase coordinates.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        INTEGER i

        OPEN(ubodsin,FILE=ibodfile,STATUS='OLD')

        READ(ubodsin,*) nbodies,tnow

        DO 10 i=1,nbodies
           READ(ubodsin,*) mass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i)
 10     CONTINUE

        CLOSE(ubodsin)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpars
C
C
C***********************************************************************
C
C
C     Subroutine to initialize system parameters that depend on
C     either the input data or defined PARAMETERS.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C=======================================================================

C   Initialize misc. useful numbers.
C   --------------------------------
        one=1.0
        two=2.0
        pi=4.0*ATAN(one)
        twoopi=2./pi
        onesixth=1./6.
        tiny=1.e-30
        zero=0.0

        tpos=tnow
        tvel=tnow

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initsys
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the state of the system.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        REAL second

C=======================================================================

C   Begin timing.
C   -------------
        cputime0=SECOND()

        CALL startout
C            --------
        CALL inparams
C            --------
        CALL inbods
C            ------
        CALL checkinp
C            --------
        CALL initpars
C            --------

c

        CALL accpot
C            ------
        IF(fixacc) CALL corracc
C                       -------

        CALL outstate(0)
C            --------

        CALL initvel
C            -------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initvel
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the velocities of the bodies for the
C     initial timestep.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
!     $omp parallel do
        DO 10 p=1,nbodies
           vx(p)=vx(p)+0.5*dtime*ax(p)
           vy(p)=vy(p)+0.5*dtime*ay(p)
           vz(p)=vz(p)+0.5*dtime*az(p)
 10     CONTINUE
!     $omp end parallel do

C   Update velocity time, system time.
C   ----------------------------------
        tvel=tvel+0.5*dtime
        tnow=tvel

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE inparams
C
C
C***********************************************************************
C
C
C     Subroutine to read in parameters.
C
C     Input parameters:
C
C        headline  : identification string for the run.
C        nsteps    : number of timesteps.
C        noutbod   : output system state once every nsteps/noutbod 
C                    steps.
C        noutlog   : output logfile data once every nsteps/noutlog
C                    steps.
C        dtime     : the timestep.
C        G         : value of gravitational constant, in appropriate
C                    units.
C        selfgrav  : option to turn off (.FALSE.) system self-gravity.
C        inptcoef  : option to read-in expansion coefficients.
C        outpcoef  : option to write-out expansion coefficients.
C        zeroodd   : option to zero all odd terms in the expansion.
C        zeroeven  : option to zero all even terms in the expansion.
C        fixacc    : option to force conservation of linear
C                    momentum by subtracting acceleration of c.o.m.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        CHARACTER *1 pcomment

C=======================================================================

        OPEN(UNIT=upars,FILE=parsfile,STATUS='OLD')

C   Read parameters, close the file.
C   --------------------------------

        READ(upars,'(a)') pcomment

        READ(upars,'(a)') headline
        READ(upars,*) nsteps
        READ(upars,*) noutbod
        READ(upars,*) noutlog
        READ(upars,*) dtime
        READ(upars,*) G
        READ(upars,*) selfgrav
        READ(upars,*) inptcoef
        READ(upars,*) outpcoef
        READ(upars,*) zeroodd
        READ(upars,*) zeroeven
        READ(upars,*) fixacc

        CLOSE(UNIT=upars)
 
        RETURN
        END
C***********************************************************************
C
C
                     SUBROUTINE iocoef(sinsum,cossum)
C
C
C***********************************************************************
C
C
C     Subroutine to input and output expansion coefficients.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        INTEGER n,l,m
        LOGICAL firstc
        REAL*8 sinsum,cossum,tt

        DIMENSION sinsum(0:nmax,0:lmax,0:lmax),
     &            cossum(0:nmax,0:lmax,0:lmax)

        DATA firstc/.TRUE./

        SAVE firstc

C=======================================================================

        IF(firstc) THEN

           firstc=.FALSE.

           IF(outpcoef) OPEN(uoutcoef,FILE=outcfile,STATUS='NEW')
           IF(inptcoef) OPEN(uincoef,FILE=incfile,STATUS='OLD')

        ENDIF

        IF(outpcoef) THEN

           WRITE(uoutcoef,100) tnow

           DO 30 n=0,nmax
              DO 20 l=0,lmax
                 DO 10 m=0,l
                    WRITE(uoutcoef,100) sinsum(n,l,m),cossum(n,l,m)
 10              CONTINUE
 20           CONTINUE
 30        CONTINUE

 100       FORMAT(1x,10(1pe22.13))

        ENDIF

        IF(inptcoef) THEN

           READ(uincoef,*) tt

           IF(tt.NE.tnow) CALL terror(' input error in iocoef ')
C                              ------

           DO 130 n=0,nmax
              DO 120 l=0,lmax
                 DO 110 m=0,l
                    READ(uoutcoef,*) sinsum(n,l,m),cossum(n,l,m)
 110             CONTINUE
 120          CONTINUE
 130       CONTINUE

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE outbods
C
C
C***********************************************************************
C
C
C     Subroutine to output phase coordinates.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        CHARACTER*3 sstring
        CHARACTER*7 filename
        CHARACTER*8 filenam1
        CHARACTER*8 filepar
        CHARACTER*10 nstring

        INTEGER i,istring,nsnap
        REAL*8 lx,ly,lz,energy,ek,ep

        SAVE nsnap,nstring

        DATA nsnap/0/,nstring/'0123456789'/

C=======================================================================

        nsnap=nsnap+1

        sstring(1:1)=nstring(1+nsnap/100:1+nsnap/100)
        istring=1+MOD(nsnap,100)/10
        sstring(2:2)=nstring(istring:istring)
        istring=1+MOD(nsnap,10)
        sstring(3:3)=nstring(istring:istring)
        filepar=obodfile
        filename=filepar(1:4)//sstring(1:3)

        OPEN(UNIT=ubodsout,FILE=filename,STATUS='NEW')

        filepar=elfile
        filenam1=filepar(1:5)//sstring(1:3)

        OPEN(UNIT=ubodsel,FILE=filenam1,STATUS='NEW')

        WRITE(ubodsout,20) nbodies,tnow
 20     FORMAT(1x,1i6,1pe14.6)

        DO 30 i=1,nbodies
           WRITE(ubodsout,111) mass(i),x(i),y(i),z(i),vx(i),vy(i),vz(i),
     &                         pot(i)+potext(i)
           lx=mass(i)*(y(i)*vz(i)-z(i)*vy(i))
           ly=mass(i)*(z(i)*vx(i)-x(i)*vz(i))
           lz=mass(i)*(x(i)*vy(i)-y(i)*vx(i))
           ek=0.5*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
           ep=mass(i)*pot(i)+mass(i)*potext(i)
           energy=ek+ep
           WRITE(ubodsel,111) tnow,lx,ly,lz,energy
 30     CONTINUE

 111    FORMAT(1x,10(1pe14.6))

        CLOSE(ubodsout)
        CLOSE(ubodsel)

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE outlog
C
C
C***********************************************************************
C
C
C     Subroutine to output phase coordinates.
C
C
C=======================================================================

        INCLUDE 'scf.h'

        INTEGER i
        LOGICAL firstc
        REAL second
        REAL*8 lxtot,lytot,lztot,mtot,vxcm,vycm,vzcm,etot,ektot,eptot,
     &         m2tw,t1,clausius,m2claus,cpux,xcm,ycm,zcm,epselfg

        DATA firstc/.TRUE./

        SAVE firstc

C=======================================================================

        IF(firstc) THEN
           firstc=.FALSE.
           cputime=cputime0
        ENDIF

        OPEN(UNIT=ulog,FILE=logfile,STATUS='UNKNOWN')

        DO 10 i=1,2*nsteps
           READ(ulog,120,end=20) t1
 10     CONTINUE
 
 20     CONTINUE

        mtot=0.0
        xcm=0.0
        ycm=0.0
        zcm=0.0
        vxcm=0.0
        vycm=0.0
        vzcm=0.0
        lxtot=0.0
        lytot=0.0
        lztot=0.0
        etot=0.0
        ektot=0.0
        eptot=0.0
        epselfg=0.0
        clausius=0.0

        cpux=SECOND()

!     $omp parallel do
        DO 30 i=1,nbodies
           mtot=mtot+mass(i)
           xcm=xcm+mass(i)*x(i)
           ycm=ycm+mass(i)*y(i)
           zcm=zcm+mass(i)*z(i)
           vxcm=vxcm+mass(i)*vx(i)
           vycm=vycm+mass(i)*vy(i)
           vzcm=vzcm+mass(i)*vz(i)
           lxtot=lxtot+mass(i)*(y(i)*vz(i)-z(i)*vy(i))
           lytot=lytot+mass(i)*(z(i)*vx(i)-x(i)*vz(i))
           lztot=lztot+mass(i)*(x(i)*vy(i)-y(i)*vx(i))
           eptot=eptot+0.5*mass(i)*pot(i)+mass(i)*potext(i)
           epselfg=epselfg+0.5*mass(i)*pot(i)
           ektot=ektot+0.5*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
           clausius=clausius+mass(i)*(x(i)*ax(i)+y(i)*ay(i)+z(i)*az(i))
 30     CONTINUE
!     $omp end parallel do

        xcm=xcm/mtot
        ycm=ycm/mtot
        zcm=zcm/mtot
        vxcm=vxcm/mtot
        vycm=vycm/mtot
        vzcm=vzcm/mtot

        etot=ektot+eptot
        m2tw= -2.*ektot/eptot
        m2claus= -2.*ektot/clausius

        BACKSPACE(ulog)
        WRITE(ulog,120) tnow
        WRITE(ulog,120) mtot
        WRITE(ulog,120) xcm,ycm,zcm
        WRITE(ulog,120) vxcm,vycm,vzcm
        WRITE(ulog,120) lxtot,lytot,lztot
        WRITE(ulog,120) ektot,eptot,epselfg,etot
        WRITE(ulog,120) m2tw,clausius,m2claus
        WRITE(ulog,120) cpux-cputime
        WRITE(ulog,130)

 120    FORMAT(5(1pe18.10))
 130    FORMAT(/)

        CLOSE(ulog)

        cputime=cpux

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outstate(n)
C
C
C***********************************************************************
C
C
C     Subroutine to output information about the system state to
C     the log and body data files.
C
C
C=======================================================================
 
        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

        CALL outterm(' step completed: ',n)
C            -------

        IF(n.EQ.0) THEN

           CALL outbods
C               -------
           CALL outlog
C               ------

        ELSE

           IF(MOD(n,noutlog).EQ.0.OR.(MOD(n,noutbod).EQ.0)) THEN

              CALL corrvel('correct')
C                  -------

              IF(MOD(n,noutlog).EQ.0) CALL outlog
C                                          ------
              IF((MOD(n,noutbod).EQ.0)) CALL outbods
C                                            -------
              CALL corrvel('reset  ')
C                  -------
           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE outterm(message,n)
C
C
C***********************************************************************
C
C
C     Subroutine to output a message to the terminal and to the
C     terminal emulation file.
C
C
C=======================================================================
 
        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER n

C=======================================================================
 
        OPEN(UNIT=utermfil,FILE=termfile,STATUS='OLD')

C   Write the message.
C   ------------------

        IF(n.GE.0) THEN
           WRITE(uterm,*) message,n
           WRITE(utermfil,*) message,n
        ELSE
           WRITE(uterm,40)
           WRITE(uterm,50) message 
           WRITE(uterm,40)
           WRITE(utermfil,40)
           WRITE(utermfil,50) message 
           WRITE(utermfil,40)
        ENDIF

 40     FORMAT(/,1x,72('*'))
 50     FORMAT(/,a)

        CLOSE(UNIT=utermfil)

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE plgndrv(nbodies,l,m,dblfact,x,plm1m,plm2m,plm)
C
C
C***********************************************************************
C
C
C     A function to compute Legendre functions.
C
C
C=======================================================================

        INTEGER l,m,i,k,nbodies
        REAL*8 x,plm,plm1m,plm2m,dblfact

        DIMENSION x(1),plm(1),dblfact(1),plm1m(1),plm2m(1)

        IF(l.EQ.m.OR.l.EQ.m+1) THEN

           IF(m.EQ.0) THEN
!     $omp parallel do
              DO 20 k=1,nbodies
                 plm(k)=1.0
 20           CONTINUE
!     $omp end parallel do
           ELSE
!     $omp parallel do
              DO 30 k=1,nbodies
                 plm(k)=(-1.)**m*dblfact(m)*SQRT(1.-x(k)*x(k))**m
 30           CONTINUE
!     $omp end parallel do
           ENDIF

           IF(l.EQ.m+1) THEN
!     $omp parallel do
              DO 40 k=1,nbodies
                 plm(k)=plm(k)*x(k)*(2.*m+1.)
 40           CONTINUE
!     $omp end parallel do
           ENDIF

        ELSE

           IF(m.EQ.0) THEN
!     $omp parallel do
              DO 50 k=1,nbodies
                 plm2m(k)=1.0
 50           CONTINUE
!     $omp end parallel do
           ELSE
!     $omp parallel do
              DO 60 k=1,nbodies
                 plm2m(k)=(-1.)**m*dblfact(m)*SQRT(1.-x(k)*x(k))**m
 60           CONTINUE
!     $omp end parallel do
           ENDIF

!     $omp parallel do
           DO 70 k=1,nbodies
              plm1m(k)=plm2m(k)*x(k)*(2.*m+1.)
 70        CONTINUE
!     $omp end parallel do

           DO 90 i=m+2,l
!     $omp parallel do
              DO 80 k=1,nbodies
                 plm(k)=(x(k)*(2.*i-1.)*plm1m(k)-(i+m-1.)*plm2m(k))/
     &                  (i-m)
                 plm2m(k)=plm1m(k)
                 plm1m(k)=plm(k)
 80           CONTINUE
!     $omp end parallel do
 90        CONTINUE

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                            FUNCTION second()
C
C
C***********************************************************************
C
C
C     Subroutine to return elapsed cpu time.
C
C
C=======================================================================

        REAL etime,utime,stime,x,second

c        x=etime(utime,stime)

        second=utime+stime     

        RETURN 
        END
C***********************************************************************
C
C
                          SUBROUTINE startout
C
C
C***********************************************************************
C
C
C     Subroutine to open disk files for subsequent input/output.
C
C
C=======================================================================
 
        INCLUDE 'scf.h'

C=======================================================================
 
C   Create terminal emulation file.
C   -------------------------------
        OPEN(UNIT=utermfil,FILE=termfile,STATUS='UNKNOWN')
        WRITE(utermfil,*) ' Start of output '
        CLOSE(UNIT=utermfil)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE steppos
C
C
C***********************************************************************
C
C
C     Subroutine to advance the positions of the bodies for a timestep.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------

!     $omp parallel do
        DO 10 p=1,nbodies
           x(p)=x(p)+vx(p)*dtime
           y(p)=y(p)+vy(p)*dtime
           z(p)=z(p)+vz(p)*dtime
 10     CONTINUE
!     $omp end parallel do

C   Update position time, system time.
C   ----------------------------------
        tpos=tpos+dtime
        tnow=tpos

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE stepsys(n)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the state of the system by one timestep.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n, i
        REAL*8 tmp

C=======================================================================

        CALL steppos
C            -------
        CALL accpot
C     add the black hole
!     $omp parallel do
        DO 5 i=1,nbodies
           tmp=1D-4*n/r(i)**3
           ax(i)=ax(i)-dum*x(i)
           ay(i)=ay(i)-dum*y(i)
           az(i)=az(i)-dum*z(i)
 5      CONTINUE
!     $omp end parallel do

C            ------
        IF(fixacc) CALL corracc
C                       -------
        CALL stepvel
C            -------

        CALL outstate(n)
C            --------

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE stepvel
C
C
C***********************************************************************
C
C
C     Subroutine to advance the velocities of the bodies for timestep.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================

C   Loop over all velocity components for all bodies.
C   -------------------------------------------------

!     $omp parallel do
        DO 10 p=1,nbodies
           vx(p)=vx(p)+ax(p)*dtime
           vy(p)=vy(p)+ay(p)*dtime
           vz(p)=vz(p)+az(p)*dtime
 10     CONTINUE
!     $omp end parallel do

C   Update velocity time, system time.
C   ----------------------------------
        tvel=tvel+dtime
        tnow=tvel

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stopout
C
C
C***********************************************************************
C
C
C     Subroutine to end output.
C
C
C=======================================================================
 
        INCLUDE 'scf.h'

        INTEGER i
        REAL*8 t1

C=======================================================================
 
        OPEN(UNIT=ulog,FILE=logfile,STATUS='UNKNOWN')

        DO 10 i=1,2*nsteps
           READ(ulog,120,end=20) t1
 10     CONTINUE
 
 20     CONTINUE

        BACKSPACE(ulog)
        WRITE(ulog,30) cputime1-cputime0
 30     FORMAT(//,15x,' Total cpu time used (secs) =',1pe15.7)

 120    FORMAT(15(1pe13.5))

        CLOSE(ulog)

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE terror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to terminate the program as the result of a fatal
C     error, close the output files, and dump timing information.
C
C
C=======================================================================

        INCLUDE 'scf.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER ierror
        REAL second

C=======================================================================

C   Write error message to the log file and to the terminal.
C   --------------------------------------------------------
        ierror=-1

        CALL outterm(message,ierror)
C            -------

C-----------------------------------------------------------------------
C   Stop timing, output timing data, close files, terminate the
C   simulation.
C-----------------------------------------------------------------------

        cputime1=SECOND()

        STOP
        END
C***********************************************************************
C
C
        SUBROUTINE ultrasv(nbodies,n,alpha,xi,un,unm1,ultrasp)
C
C
C***********************************************************************
C
C
C     A function to compute ultraspherical polynomials.
C
C
C=======================================================================

        INTEGER n,i,k,nbodies
        REAL*8 ultrasp,alpha,xi,unm1,un,twoalpha,c1,c2,c3

        DIMENSION xi(1),un(1),unm1(1),ultrasp(1)

        twoalpha=2.0*alpha

        IF(n.EQ.0) THEN
!     $omp parallel do
           DO 20 k=1,nbodies
              ultrasp(k)=1.0
 20        CONTINUE
!     $omp end parallel do
        ELSE
           IF(n.EQ.1) THEN
!     $omp parallel do
              DO 30 k=1,nbodies
                 ultrasp(k)=twoalpha*xi(k)
 30           CONTINUE
!     $omp end parallel do
           ELSE

!     $omp parallel do
              DO 40 k=1,nbodies
                 unm1(k)=1.0
                 un(k)=twoalpha*xi(k)
 40           CONTINUE
!     $omp end parallel do

              DO 60 i=1,n-1
                 c1=2.*(i+alpha)
                 c2=i+twoalpha-1.
                 c3=1./(i+1.)
!     $omp parallel do
                 DO 50 k=1,nbodies
                    ultrasp(k)=(c1*xi(k)*un(k)-c2*unm1(k))*c3
                    unm1(k)=un(k)
                    un(k)=ultrasp(k)
 50              CONTINUE
!     $omp end parallel do
 60           CONTINUE

           ENDIF
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
        SUBROUTINE ultrasv1(nbodies,n,alpha,xi,un,unm1,ultrasp,ultrasp1)
C
C
C***********************************************************************
C
C
C     A function to compute ultraspherical polynomials.
C
C
C=======================================================================

        INTEGER n,i,k,nbodies
        REAL*8 ultrasp,alpha,xi,unm1,un,twoalpha,c1,c2,c3,ultrasp1

        DIMENSION xi(1),un(1),unm1(1),ultrasp(1),ultrasp1(1)

        twoalpha=2.0*alpha

        IF(n.EQ.0) THEN
!     $omp parallel do
           DO 20 k=1,nbodies
              ultrasp(k)=1.0
              ultrasp1(k)=0.0
 20        CONTINUE
!     $omp end parallel do
        ELSE
           IF(n.EQ.1) THEN
!     $omp parallel do
              DO 30 k=1,nbodies
                 ultrasp(k)=twoalpha*xi(k)
                 ultrasp1(k)=1.0
 30           CONTINUE
!     $omp end parallel do
           ELSE

!     $omp parallel do
              DO 40 k=1,nbodies
                 unm1(k)=1.0
                 un(k)=twoalpha*xi(k)
 40           CONTINUE
!     $omp end parallel do

              DO 60 i=1,n-1
                 c1=2.*(i+alpha)
                 c2=i+twoalpha-1.
                 c3=1./(i+1.)
!     $omp parallel do
                 DO 50 k=1,nbodies
                    ultrasp(k)=(c1*xi(k)*un(k)-c2*unm1(k))*c3
                    unm1(k)=un(k)
                    un(k)=ultrasp(k)
 50              CONTINUE
!     $omp end parallel do
 60           CONTINUE

!     $omp parallel do
              DO 70 k=1,nbodies
                 ultrasp1(k)=((twoalpha+n-1.)*unm1(k)-n*xi(k)*
     &                       ultrasp(k))/(twoalpha*(1.-xi(k)*xi(k)))
 70           CONTINUE
!     $omp end parallel do

           ENDIF
        ENDIF

        RETURN
        END

