
C=======================================================================
C
C
C                        INCLUDE FILE scf.h
C
C
C=======================================================================
C
C
C     Parameter declarations, allocation of array storage, common
C     block definitions.
C
C
C=======================================================================

        INTEGER nbodsmax,nmax,lmax
 
cXXX	want 1,000,000, 6+, 0 for our runs

        PARAMETER(nbodsmax=100000,nmax=6,lmax=6)

        CHARACTER*50 headline
        INTEGER nsteps,noutbod,nbodies,noutlog,irsort
cXXX
	PARAMETER(nmgrpmax=20)
        parameter(rpc=3.08567802d18,year=3.155674d7)
        parameter(G0=6.668d-8)
c
	INTEGER indexm(nbodsmax),iistar(nbodsmax)
	INTEGER ndumm,nmgrp
	DOUBLE PRECISION ximf,sigma,r0,t0,tscale,refm
	DOUBLE PRECISION tlim(nmgrpmax),masslim(nmgrpmax)
	DOUBLE PRECISION vbar(nmax),vsig(nmax)

        LOGICAL selfgrav,inptcoef,outpcoef,zeroodd,zeroeven,fixacc
cXXX	IBM wants double precision, not real - change for Cray?
        DOUBLE PRECISION tnow,x,y,z,vx,vy,vz,mass,pot,dtime,G,temp1,temp2,temp3,
     &         temp4,temp5,temp6,r,phi,xi,costh,clm,dlm,elm,flm,plm,
     &         dplm,ultrasp,ultrasp1,ax,ay,az,one,pi,twoopi,onesixth,
     &         tpos,tvel,cputime0,cputime1,cputime,potext,
     &         two,zero,tiny

        COMMON/bodscom/x(nbodsmax),y(nbodsmax),z(nbodsmax),vx(nbodsmax),
     &                 vy(nbodsmax),vz(nbodsmax),mass(nbodsmax),
     &                 pot(nbodsmax),ax(nbodsmax),ay(nbodsmax),
     &                 az(nbodsmax),r(nbodsmax),phi(nbodsmax),
     &                 xi(nbodsmax),costh(nbodsmax),potext(nbodsmax)
        COMMON/tempcom/temp1(nbodsmax),temp2(nbodsmax),temp3(nbodsmax),
     &                 temp4(nbodsmax),temp5(nbodsmax),temp6(nbodsmax)
        COMMON/expcom/clm(nbodsmax),dlm(nbodsmax),elm(nbodsmax),
     &                flm(nbodsmax)
        COMMON/sortcom/irsort(nbodsmax)
        COMMON/funccom/ultrasp(nbodsmax),ultrasp1(nbodsmax),
     &                 plm(nbodsmax),dplm(nbodsmax)
        COMMON/parcomi/nbodies,nsteps,noutbod,noutlog
        COMMON/parcomr/dtime,G,one,pi,twoopi,onesixth,two,tiny,zero
        COMMON/parcomc/headline
        COMMON/parcoml/selfgrav,inptcoef,outpcoef,zeroodd,zeroeven,
     &                 fixacc
        COMMON/timecom/tpos,tnow,tvel
        COMMON/cpucom/cputime0,cputime1,cputime
cXXX
	COMMON/inmcom/ximf,sigma,r0,t0,tscale,refm
	COMMON/mascom/tlim,masslim,vbar,vsig
	COMMON/indcom/nmgrp,ndumm,indexm,iistar

C=======================================================================
C   Definitions specific to input/output.
C=======================================================================
        INTEGER uterm,upars,ulog,ubodsin,ubodsout,utermfil,uoutcoef,
     &          uincoef,ubodsel
        CHARACTER*8 parsfile,logfile,ibodfile,obodfile,termfile,
     &              outcfile,incfile,elfile

        PARAMETER(uterm=6,upars=10,ulog=11,ubodsin=12,ubodsout=13,
     &            utermfil=15,uoutcoef=16,uincoef=17,ubodsel=18)
        PARAMETER(parsfile='SCFPAR',logfile='SCFLOG',
     &            ibodfile='SCFBI',obodfile='SNAPxxxx',
     &            termfile='SCFOUT',outcfile='SCFOCOEF',
     &            incfile='SCFICOEF',elfile='SCFELxxx')
