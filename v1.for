c  program general CADSOL
c             26 Nov 13

 
C RUNS WITH CDSE01

      implicit character*1  (A-Z)
      INTEGER  NX,NY,PD,MDW,MIW,IBIG,ISMALL,SOMA,ID,IDP,ILIN1,ILIN14
      INTEGER  NK,N,MV,MNLS,MCOM,M,NDQP     
      LOGICAL LX, LY, LPOL
C
C***  PLEASE CALCULATE CAREFULLY THE DIMENSIONS OF THE ARRAYS.
C***  INCORRECT DIMENSIONS MAY LEAD TO CHAOTIC RESULTS.
C?
      PARAMETER ( Nx= 251 ,  Ny= 30 , 
     &     PD = 6 ,
     $     MDW= 498085,
     &     MIW= 121074, 
     &     IBIG=  99333040,
     &     ISMALL= 2337192,
     &     SOMA=0,
     &     ID=641,
     &     IDP=3321
     &)

      PARAMETER (
     &     ILIN1 = 0,
     &     ILIN14 = 0,
     &     NK = 5 , 
     &     N = NX*NY,
     &     MV = N, 
     &     MNLS = 100,
     &     mcom=  211920,
c acesta e parametrul ce trebuie schimbat pt o mesa mai mare
     &     M= 63545,
     &     ndqp = 4002141,
     &     LX = .FALSE.,
     &     LY = .FALSE.,
     &     LPOL = .FALSE.)

      INTEGER  MDINF,MLINF,CI1,CI2,CI3,LOUT     

      PARAMETER (
     &     MDINF=5+2*NK, MLINF=18+NK, CI1=2*PD*(2*PD+1)-1,
     &     CI2=N/MV+1, CI3=CI1*CI2, 
     &     LOUT=6)

      INTEGER  NCC,MBIG,MSMALL,IN,INR,INR2,IM,IS,NLE,NLE2,GL,PD2,RL4,
     &         VPTOT,I,IERR,I1,I2,I3,IND,MAXIT,IDOKU,NCOM,S,
     &         INDE,NS,IINFO(30),ILIN(17),IFILE(14),IWORK(MIW),
     &         anx
      INTEGER  HNX,HNY,HN,HM,HMV

      integer ncont, irun, method 
      double precision gam,x,y,pi,h1,f2,sh,ch
      double precision
     &         TIME2,DINFO(MDINF),U1(n,nk),U2(n,NK),
     &         DWORK(MDW),
     &         RE,errlim

      LOGICAL LINFO(MLINF), LLIN(M,2), LSTART
      EXTERNAL CDSU01, CDSU02, CDSU03, CDSU04, FDS500, cdse01

C
      COMMON /CDSCON/ RE
C
C     COMMON BLOCKS FOR PACKED STORING OF THE DIAGONAL PARTS
C
      COMMON /CDSINT/NCOM,INDE(MCOM)
      COMMON /NINT/NS(CI3)
C
C     COMMON BLOCKS FOR INCORE FILES
C
      double precision
     &         BIG
      INTEGER  SMALL, methoda
C
      INTEGER  NBIG,NTOTAL,IOINF1,IOINF2,IOINF3
      COMMON /FDSIO/ BIG(IBIG)
      COMMON /FDSINF/ NBIG,NTOTAL,IOINF1(99),IOINF2(99),IOINF3(99)
C
      INTEGER  NSMALL,NTOT,IOINF4,IOINF5,IOINF6
      COMMON /CDSIO/ SMALL(ISMALL)
      COMMON /CDSINF/ NSMALL,NTOT,IOINF4(99),IOINF5(99),IOINF6(99)
C
C     COMMON BLOCKS FOR PACKED STORING OF THE MATRIX
C
      INTEGER  IDQP,NIND,NPACK 
      double precision
     &         DQP
      COMMON/FDSLEN/ NPACK(ID)
      COMMON/FDSIND/ NIND,IDQP(NDQP)
      COMMON/FDSDQP/ DQP(M)

 
      integer norming

c arrays for storing values of derivs calculated by fidisol/cadsol
      double precision u0x,u0y,u0xx,u0yy,u0xy,u0
      integer m_nx,m_ny
      logical fifo_exist
      common/derivs/ u0x(n,nk),u0y(n,nk),u0xx(n,nk),u0yy(n,nk),
     *     u0xy(n,nk),u0(n,nk)
      common/dimensions/m_nx,m_ny

        double precision nr,w,con,alpha,c1,c2,c3,lambda,rh,par
  
      double precision rout(nk) 
      integer ipo(nk),ii,maxiti

      common/parae/ eng
      common nr,w,con,alpha,c1,c2,c3,lambda,rh,par
 
       
 
c      open(7,file='skyrme.fifo',err=661,status='old')
 661  continue
C
      m_nx=nx
      m_ny=ny

      NIND = NDQP
      NCOM = MCOM
C
      RE = 1.
C
      HNX = NX
      HNY = NY
      IF(LX) HNX = HNX - 1
      IF(LY) HNY = HNY - 1
      HN = HNX * HNY
      HM = HN * NK
      HMV = MV
      IF(MV.EQ.N) HMV = HN
C
C     CALCULATION OF THE REQUIRED INCORE SPACE (SEE IOINF1 AND IOINF4).
C     (THE INCORE SPACE IN IOINF1 AND IOINF4 IS CALCULATED AUTO-
C      MATICALLY. THERE IS NOTHING TO CHANGE BY THE USER.)
C     ATTENTION ! IT IS NOT POSSIBLE TO USE THE SAME ARRAY NUMBERS
C     IN IOINF1 AND IOINF4 !
C
      NCC = 2*PD*(PD+1)+1
C
      PD2  = PD + 2
      NLE  = ((PD+1)*(PD+2)) / 2
      NLE2 = ((PD2+1)*(PD2+2)) / 2
C
      IN   = HN / HMV
      IF(IN*HMV.LT.HN) IN=IN+1
      IM = HM / HMV
      IF(IM*HMV.LT.HM) IM=IM+1
      RL4 = MAX(HNX,HNY)
      IS  = HN / MNLS
      IF(IS*MNLS.LT.HN) IS=IS+1
      INR = ( MNLS*NLE ) / HMV
      IF ( INR*HMV .LT. MNLS*NLE ) INR = INR + 1
      INR2 = ( MNLS*NLE2 ) / HMV
      IF ( INR2*HMV .LT. MNLS*NLE2 ) INR2 = INR2 + 1


        PI=4.d0*dATAN(1.d0)

ccccccccccccccccccccccccccccccccccccccccccccccc
c  parameters of the problem
cccccccccccccccccccccccccccccccccccccccccccccc
 
c			 rh=0.1d0
 
c         nr = 1.d0
c           w= 0.82d0 
            Open(Unit=1,File='parameters')  
              Read(1,*)nr
              Read(1,*)w
              Read(1,*)rh
              Read(1,*)alpha
            Close(Unit=1)


c      * coefficient Z^6
          c1=1.d0

c      * coefficient Z^4
          c2= -2.d0

c      * coefficient Z^2
          c3=1.1d0
 
  
c	alpha=1.d0
 

 
               methoda = 2
                maxiti = 42520  


            write(6,*) 'n= ', nr  
            write(6,*)' w= ', w    
		  write(6,*)' rh= ', rh     
		  write(6,*)' alpha= ', alpha  
		  write(6,*)' c1= ', c1
		  write(6,*)' c2= ', c2
		  write(6,*)' c3= ', c3
		  write(6,*)' '  
 
c	go to 9999
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C  !!!!!!!!!!!!!    BEGINNING OF ITERATION     !!!!!!!!!!!!!!!!!!!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      Do 99 ii= 1, 1

	write(6,*)"i=",ii
 
c          alfa=0.02d0 + (ii-1)*0.01d0
C
      IOINF1(21) = 0
C
      IOINF1(22) = 5* IN * NLE * HMV
      IOINF1(23) = (5*NK*IN+4*IN) * HMV
      IOINF1(24) = (24*NK*(1+NK)+8) * RL4
      IOINF1(25) = 5 * NK * IN * HMV
      IOINF1(26) = MAX( IDP*IN, 5*INR2*IS ) * HMV
      IOINF1(27) = MAX( 6*NK*NK*IN, 5*INR*IS ) * HMV
      IOINF1(28) = ID * HM
      IOINF1(29) = 5 * NLE * HN
      IOINF1(30) = 5 * IN * NLE2 * HMV
      IOINF1(31) = ID * HM
      IOINF1(32) = HM
C
      IOINF4(33) = (INR*IS+IN*NLE)*HMV
      IOINF4(34) = (INR2*IS+IN*NLE2)*HMV
C
C***  IF ILIN14 > 0, THEN IOINF1(31) IS REDUCED
C***  IF ILIN14 > 1, THEN ALSO IOINF1(28) IS REDUCED
C
      IF ( (ILIN14.GT.0).AND.(SOMA.NE.0).AND.(IOINF1(31).NE.0) ) THEN
         IOINF1(31) = SOMA
         IF (ILIN14.GT.1.AND.(IOINF1(28).NE.0)) THEN
            IOINF1(28) = SOMA
         ENDIF
      ENDIF
C
      MBIG   = 0
      MSMALL = 0
      DO 1111 I=2,12
           MBIG = MBIG + IOINF1(I + 20)
 1111 CONTINUE
      DO 1112 I=13,14
           MSMALL = MSMALL + IOINF4(I + 20)
 1112 CONTINUE
C
      IF(MBIG.GT.IBIG) THEN
        WRITE(LOUT,508) IBIG,MBIG
        STOP
      ELSE IF(MBIG.LT.IBIG) THEN
        WRITE(LOUT,509) IBIG,MBIG
      ENDIF
C
      IF(MSMALL.EQ.0) THEN
        MSMALL = 1
      ENDIF
C
      IF(MSMALL.GT.ISMALL) THEN
        WRITE(LOUT,510) ISMALL,MSMALL
        STOP
      ELSE IF(MSMALL.LT.ISMALL) THEN
        WRITE(LOUT,511) ISMALL,MSMALL
      ENDIF
C
      NBIG       = IBIG
      NSMALL     = ISMALL
      NTOTAL     = 0
      NTOT       = 0
C
C***  CALCULATION OF MAIN STORAGE IN MEGABYTES
C
      IF((ILIN1.EQ.0).OR.(ILIN1.EQ.1).OR.(ILIN1.EQ.10)) THEN
         GL = 15
      ELSE IF((ILIN1.EQ.2).OR.(ILIN1.EQ.3).OR.(ILIN1.EQ.12).OR.
     &        (ILIN1.EQ.13)) THEN
         GL = 11
      ELSE IF(ILIN1.EQ.4) THEN
         GL = 10
      ELSE IF(ILIN1.EQ.11) THEN
         GL =  9
      ELSE
         WRITE(LOUT,512) ILIN1
         STOP
      ENDIF
C
      S = 2 * (MIW+NDQP) + 5*N*NK + MDW + MCOM
C
      VPTOT = INT( (IBIG + (GL+3)*M + S)/1.E6 * 8 + ISMALL/1.E6 * 4) + 1
C
      WRITE(LOUT,513) VPTOT
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***  INITIAL SOLUTION

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	

	write(6,*)
	write(6,*)"_________________________________"
	write(6,*)'start iteration ' 

            write(6,*) 'n= ', nr  
            write(6,*)' w= ', w    
		  write(6,*)' rh= ', rh     
		  write(6,*)' alpha= ', alpha  
		  write(6,*)' c1= ', c1
		  write(6,*)' c2= ', c2
		  write(6,*)' c3= ', c3
		  write(6,*)' '  

 
	write(6,*)

c    	go to 666

 	open(13,file='funct.dat')

        DO 42 I2=1,ny
        DO 41 I1=1,nx


           anx=i1+(i2-1)*nx 
           READ(13,765) dwork(anx),dwork(anx+n),(U1(anx,i3),i3=1,nk)
c         write(6,*)I2,I1
   41   CONTINUE
        read(13,*)
   42   continue
 
      Close(Unit=13)  

 
c	dwork(nx)=1.d0

      LSTART = .TRUE.

C***  UNIT NUMBERS
C
      DO 50 I=2,14
      IFILE(I) = I+20
   50 CONTINUE
C
C***  INFORMATION VECTORS
C
      idoku = 1000 
c      IDOKU = 1
c      MAXIT = 0
C
      IINFO(1)  = NK
      IINFO(2)  = NX
      IINFO(3)  = NY
      IINFO(5)  = MV
      IINFO(6)  = IDOKU
      IINFO(7)  = LOUT
      IINFO(8)  = MAXIT
      IINFO(9)  = MDW
      IINFO(10) = MIW
      IINFO(11) = MCOM
      IINFO(12) = PD
      IINFO(13) = MNLS
	 
	errlim=0.01d0
C
      DINFO(1)  = errlim
C
      LINFO(1)  = LSTART
      LINFO(2)  = LX
      LINFO(3)  = LY
      LINFO(5)  = .true.

c shnir ia .true.
c        LINFO(5) = .true.
c  da rezultate foarte diferite 
  
      LINFO(6)  = .FALSE.
      LINFO(7)  = .TRUE.
      LINFO(8)  = .FALSE.
c      LINFO(9)  = .TRUE.
      LINFO(9)  = .FALSE.
c nu vreau full output screen->   LINFO(9)  = .FALSE.
      LINFO(10) = .FALSE.  
C
      DO 100 I=1,NK
      LINFO(10+I) = .FALSE.
  100 CONTINUE
C
      LINFO(15+NK) = .FALSE.
      LINFO(16+NK) = .FALSE.
      LINFO(17+NK) = .FALSE.
      LINFO(18+NK) = LPOL
C
C
c         method
c      ILIN(1)  = 10, 2
c  cu 4 converge f rapid; dar crashes la computation xy defect
c                             la fel cu method 10

        ILIN(1)  = methoda

c  maximum iterations
      ILIN(2)  =    maxiti

      ILIN(14) = ILIN14
      ILIN(15) = 0

 
      CALL CDSE01(IINFO,DINFO,LINFO,ILIN,LLIN,IWORK,DWORK,IFILE,U1,U2,
     &            IND,IERR,CDSU01,CDSU02,CDSU03,CDSU04,FDS500)
c      CALL SECOND(TIME2)
c      TIME2 = TIME2 - TIME1
      WRITE(LOUT,501) TIME2, IND, IERR, ILIN(10)
 
 
 
 

c  U1 = F1
c  U2 = F2
c  U3 = F0
c  U4 = Z
c  U5 = W
 
  
        rewind(1)
        open(unit=1,File='funct.dat')


c first angular direction: theta =0 -- Z=0
            i2=1
         DO 1171 I1=1,Nx-1
             anx=i1+(i2-1)*nx
c        write(6,*)I1,dwork(anx)
            write(1,765)dwork(anx),0.,
     & u1(anx,1),u1(anx,2),u1(anx,3),0.,u1(anx,5) 

 1171      CONTINUE 
                I1=Nx
              anx=i1+(i2-1)*nx
c last radial point: infinity--all functions are zero
            write(1,765)1., dwork(N+anx), 0.,0.,0.,0.,0.
          write(1,*) 

        DO 170 I2=2,Ny

          DO 171 I1=1,Nx-1
             anx=i1+(i2-1)*nx
c        write(6,*)I1,dwork(anx)
            write(1,765)dwork(anx),dwork(N+anx),(u1(anx,i3),i3=1,nk)

 171      CONTINUE 
                I1=Nx
              anx=i1+(i2-1)*nx
c last radial point: infinity--all functions are zero
            write(1,765)1., dwork(N+anx), 0.,0.,0.,0.,0.
          write(1,*) 
 170    CONTINUE

             Close(Unit=1)


c       rewind(2)
        open(unit=1,File='err.dat')
        DO 270 I2=1,Ny
          DO 271 I1=1,Nx
             anx=i1+(i2-1)*nx

            write(1, 5003)dwork(anx),dwork(N+anx),(u2(anx,i3),i3=1,nk)

 271      CONTINUE 
          write(1,*) 
 270    CONTINUE

             Close(Unit=1)


	   rewind(1)
        open(unit=1,File='gridx.dat')
 
          DO 471 I1=1,Nx
             anx=i1 

            write(1,5003)dwork(anx)/(1.-dwork(anx)+1.d-12) 
 471      CONTINUE  

             Close(Unit=1)

		   rewind(2)
        open(unit=2,File='gridy.dat')
 
          DO 472 I2=1,Ny  
	y=pi/2.*( (i2-1.)/(ny-1.))
c		write(6,*)y,pi,i2 
            write(2,5003)y
 472      CONTINUE  

             Close(Unit=2)

	   rewind(1)
	  open(unit=1,File='functf.dat')
        DO 180 I2=1,Ny
          DO 181 I1=1,Nx
             anx=i1+(i2-1)*nx

            write(1,*) u1(anx,1) 
	      write(1,*) u1(anx,2)
	      write(1,*) u1(anx,3)
	      write(1,*) u1(anx,4)
	      write(1,*) u1(anx,5) 

 181      CONTINUE          
 180    CONTINUE
             Close(Unit=1)


            write(6,*) 'n= ', nr  
            write(6,*)' w= ', w    
		  write(6,*)' rh= ', rh     
		  write(6,*)' alpha= ', alpha  
		  write(6,*)' c1= ', c1
		  write(6,*)' c2= ', c2
		  write(6,*)' c3= ', c3
		  write(6,*)' '  

	      write(6,*)' '  
			      write(6,*)' ' 

 
	write(6,*)"___________________________________"
	write(6,*)" "
 
   99   continue
9999     continue

            open(unit=5,File='res.txt')
 	        WRITE(5,765)   nr,w,alpha,c1,c2,c3,rh
 	        Close(Unit=5)

C
c            PAUSE 'Press ENTER to continue'  
      STOP
C
  501 FORMAT('  TIME FOR UP CDSE01:',G15.6,' SECONDS'/
     &       '  IND,IERR,IMVM:',3I10)
  502 FORMAT(/' '/'   SOLUTIONS OF THE EXAMPLE :'/
     &       '    IX    IY',5X,'COORDINATES',15X,
     &       ' SOL. U',7X,'SOL. V',7X,'SOL. W',/' ')
  503 FORMAT(2I6,' : (',G10.3,',',G10.3,') - ',3(G11.4,2X))
  504 FORMAT(/' '/'   ESTIMATES OF THE ERRORS (RELATIVE',
     &       ' ERRORS) OF THE EXAMPLE :'/
     &       '    IX    IY',5X,'COORDINATES',15X,
     &       ' ERR. U',7X,'ERR. V',7X,'ERR. W',/' ')
  508 FORMAT(//'  ERROR IN ARRAY BIG : ',/,
     &       //'  IBIG = ',I10,' MUST BE GREATER/EQUAL:  ',I10)
  509 FORMAT(//'  IBIG = ',I10,' CAN BE REDUCED TO: ',I10)
  510 FORMAT(//'  ERROR IN ARRAY SMALL : ',/,
     &       //'  ISMALL = ',I10,' MUST BE GREATER/EQUAL:  ',I10)
  511 FORMAT(//'  ISMALL = ',I10,' CAN BE REDUCED TO: ',I10)
  512 FORMAT (//,' ILIN1 = ',I4,' THIS METHOD FOR LINEAR EQUATION ',
     &             'SOLVER DOES NOT EXIST')
  513 FORMAT(//'  NEED OF MAIN STORAGE OF YOUR PROGRAM : ',I10,' MBYTE')
  514 FORMAT(2(/1H ))
  612 FORMAT(/' MAXIMUM OF THE ESTIMATED RELATIVE ERRORS:')
  613 FORMAT(' SOLUTION COMPONENT',I2,':',G13.4,'  AT POINT',I8,
     &       ' = (',G10.3,',',G10.3,')')

 5000 FORMAT (7(e24.16/))
 5003 FORMAT(13(e24.16))
 5004 FORMAT( (e24.16))
c 765     format (9(e24.16))
  765  FORMAT (7(e24.16,1x))
C
c END OF P R O G R A M
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CX
C X
C  X
C   X
C    X
C     X
C      X
C       X
C        X
C         X
C          X
C           X
C            X
C             X
C              X
C               X
C                X
C                 X
C                  X
C                   X
C                    X
C                     X
C                      X
C                       X
C                        X
C                         X
C                          X
C                           X
C                            X
C                             X
C                              X
C                               X
C                                X
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdsu01(t,xx,y,u,ut,ux,uy,uxx,uxy,uyy,p,mt,mv,nk,nv) 
      implicit character*1 (a-z)  
      integer  mt, mv, nk, nv  
      double precision  t,xx(mv),y(mv),u(mv,nk),ut(mt,nk),ptest(nk),
     *                  ux(mv,nk),uy(mv,nk),uxx(mv,nk),uxy(mv,nk),  
     *                  uyy(mv,nk),p(mv,nk),ec11(nv),ec12(nv)   
      integer  I,k,j   
 
 

	double precision  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
      common  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
  
	Double precision 
     & U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,
     & U1x,U2x,U3x,U4x,U5x,U6x,U7x,U8x,U9x,U10x,U11x,
     & U1y,U2y,U3y,U4y,U5y,U6y,U7y,U8y,U9y,U10y,U11y,
     & U1xy,U2xy,U3xy,U4xy,U5xy,U6xy,U7xy,U8xy,U9xy,U10xy,U11xy,
     & U1xx,U2xx,U3xx,U4xx,U5xx,U6xx,U7xx,U8xx,U9xx,U10xx,U11xx,
     & U1yy,U2yy,U3yy,U4yy,U5yy,U6yy,U7yy,U8yy,U9yy,U10yy,U11yy 
 

       Double precision sn,sn2,cs,cs2,csc,ct,r,T34,T44,Ttot,ff  
       Double precision eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,e11,e12
	 Double precision r1,r2,r3,r4,r5,r6,r7,r8,r9,pat1,Z2
	 Double precision err1,err2,err3,err4,err5,err6,err7,err8,err9
	Double precision norm,rn,xn,max,rm,min,rmin,minf,rminf

	Double precision H,derH,der2H,g,Ricci,eq11,eq12,Kr

 
      err1=0.d0 
      err2=0.d0  
	norm=1.
	max=-1. 
	min=2.
	minf=2.

 

         rewind(1)
        open(unit=1,name='T44.dat')

         rewind(12)
        open(unit=12,name='Ricci.dat')

           max=0.d0

      DO 11  I = 1,NV  

	                             
	U1=U(I,1)
      U2=U(I,2)
      U3=U(I,3)
      U4=U(I,4)
      U5=U(I,5)  

      U1x=UX(I,1)
      U2x=UX(I,2) 
      U3x=UX(I,3) 
      U4x=UX(I,4) 
      U5x=UX(I,5)  

      U1y=UY(I,1)
      U2y=UY(I,2) 
      U3y=UY(I,3) 
      U4y=UY(I,4) 
      U5y=UY(I,5) 

      U1yy=UYY(I,1)
      U2yy=UYY(I,2)
      U3yy=UYY(I,3)
      U4yy=UYY(I,4)
      U5yy=UYY(I,5)  
							
	U1xx=UXX(I,1)
      U2xx=UXX(I,2) 
      U3xx=UXX(I,3) 
      U4xx=UXX(I,4) 
      U5xx=UXX(I,5)  
 
	 
         r= XX(I)/(1.-XX(I)) 

         sn=dSin(Y(I))
         cs=dCos(Y(I))
	csc=1.d0/dSin(Y(I))
         ct=dCos(Y(I))/dSin(Y(I))
	 sn2=dSin(2.*Y(I))
 	cs2=dCos(2.*Y(I))
 

	ff=1.-XX(I)
        
   
       g=r**2 + rh**2 


       H= Sqrt(g)/(Sqrt(g) + rh)

	 derH= (r*rh)/(Sqrt(g)*(Sqrt(g) + rh)**2)

	der2H=        (rh**4 - 2*r**2*rh*Sqrt(r**2 + rh**2) + 
     -    rh**3*Sqrt(r**2 + rh**2))/
     -  ((r**2 + rh**2)**1.5*(rh + Sqrt(r**2 + rh**2))**3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
 
      eq1= (0.5*derH*r*sn**2*(ff**2*U1x - 1.*ff**2*U2x))/H - 
     -  (1.*r**2*sn**2*(-1.*ff**2*U1x - 1.*ff**2*U2x + 
     -       ff**2*U3x))/g - 
     -  1.*sn**2*(-1.*r*(-2.*ff**3*U1x + ff**4*U1xx) + 
     -     ff**2*U2x*(1. + ff**2*r*U3x)) - 
     -  (1.*r*sn*(-1.*sn*U1yy + (cs + sn*U2y)*U3y))/(g*H) - 
     -  (0.25*dExp(2.*U2 - 2.*U3)*g*sn**4*
     -     (ff**4*g*H*U5x**2 + U5y**2))/(H**2*r) + 
     -  alpha**2*((-2.*dExp(2.*U1 - 2.*U2)*nr**2*
     -        r*U4**2)/(g*H) + 
     -     2.*r*sn**2*(ff**4*U4x**2 + U4y**2/(g*H)) + 
     -     (2.*dExp(2.*U1 - 2.*U3)*g*sn**2*U4**2*
     -        (nr*U5 + w)**2)/(H**2*r))    

      eq2=  (derH*ff**2*sn**2*U2x)/H + (ff**2*sn**2*U2x)/r + 
     -  (r*sn**2*(2.*ff**2*U2x + ff**2*U3x))/g + 
     -  sn**2*(-2.*ff**3*U2x + ff**4*U2x**2 + ff**4*U2xx + 
     -     ff**4*U2x*U3x) + 
     -  (sn*(2.*cs*U2y + sn*U2y**2 + sn*U2yy + 
     -       (cs + sn*U2y)*U3y))/(g*H) + 
     -  (alpha**2*U4**2*((4.*dExp(2.*U1 - 2.*U2)*
     -          nr**2)/g + 
     -       2.*dExp(2.*U1)*sn**2*
     -        (c3 + c2*U4**2 + c1*U4**4)))/H + 
     -  (0.5*dExp(2.*U2 - 2.*U3)*g*sn**4*
     -     (ff**4*g*H*U5x**2 + U5y**2))/(H**2*r**2)  


      eq3=   2.*ff**2*sn*U2x - (2.*ff**2*r**2*sn*U2x)/g + 
     -  4.*ff**2*sn*U3x + 2.*ff**4*r*sn*U2x*U3x + 
     -  2.*ff**4*r*sn*U3x**2 + 
     -  (derH*r*sn*(ff**2*U2x + 3.*ff**2*U3x))/H + 
     -  2.*r*sn*(-2.*ff**3*U3x + ff**4*U3xx) + 
     -  (2.*r*((cs + sn*U2y)*U3y + sn*U3y**2 + sn*U3yy))/(g*H) - 
     -  (1.*dExp(2.*U2 - 2.*U3)*g*sn**3*
     -     (ff**4*g*H*U5x**2 + U5y**2))/(H**2*r) + 
     -  (alpha**2*sn*U4**2*
     -     (4.*dExp(2.*U1)*r*
     -        (c3 + c2*U4**2 + c1*U4**4) - 
     -       (8.*dExp(2.*U1 - 2.*U3)*g*
     -          (nr*U5 + w)**2)/(H*r)))/H  
     

      eq4=  (-1.*dExp(2.*U1 - 2.*U2)*nr**2*U4)/H - 
     -  (1.*dExp(2.*U1)*g*sn**2*U4*
     -     (c3 + 2.*c2*U4**2 + 3.*c1*U4**4))/H + 
     -  ff**2*sn**2*(r + g*
     -      (derH/H + 1/r + ff**2*U2x + ff**2*U3x))*U4x + 
     -  g*sn**2*(-2.*ff**3*U4x + ff**4*U4xx) + 
     -  (sn*(cs + sn*U2y + sn*U3y)*U4y)/H + (sn**2*U4yy)/H + 
     -  (dExp(2.*U1 - 2.*U3)*g**2*sn**2*U4*
     -     (nr*U5 + w)**2)/(H**2*r**2)  
 

      eq5=  -1.*ff**2*r*sn**2*U5x + (5.*ff**2*r**3*sn**2*U5x)/g + 
     -  r**2*sn**2*(-2.*ff**3*U5x + 3.*ff**4*U2x*U5x - 
     -     1.*ff**4*U3x*U5x + ff**4*U5xx) + 
     -  (r**2*sn*((3.*(cs + sn*U2y) - 1.*sn*U3y)*U5y + sn*U5yy))/
     -   (g*H) - (8.*dExp(2.*U1 - 2.*U2)*alpha**2*
     -     nr*r**2*U4**2*(nr*U5 + w))/(g*H)  
 
      
	p(I,1)=eq1          
      p(I,2)=eq2  
      p(I,3)=eq3  
      p(I,4)=eq4  
      p(I,5)=eq5 


      if(Abs(U4) .le. max) go to 124
	max=Abs(U4)
	rm=r
124    continue

c      	write(6,*)r,U4,max

       T34= (-2*g*nr*U4**2*(nr*U5 + w))/(dExp(2*U3)*H*r**2)
	   
 
       T44=        -1.*c2*U4**4 - 1.*c1*U4**6 - 
     -  (1.*(ff**4*g*H*U4x**2 + U4y**2))/
     -   (dExp(2.*U1)*g) + 
     -  U4**2*(-1.*c3 - (1.*csc**2*nr**2)/
     -      (dExp(2.*U2)*g) + 
     -     (g*(nr**2*U5**2 - 1.*w**2))/
     -      (dExp(2.*U3)*H*r**2))   

c Ttot= 2Td[4, 4] - Sum[Td[i, i]]
         
	   Ttot=  2.*U4**2*(c3 + c2*U4**2 + 
     -    (dExp(2.*U3)*c1*H*r**2*U4**4 - 
     -       2.*g*w*(nr*U5 + w))/
     -     (dExp(2.*U3)*H*r**2))
 

           WRITE(1,489)T34 ,T44,Ttot

	Ricci= (0.5*(-4.*dExp(2.*U3)*H**2*r*
     -       (r + ff**2*r**2*(U1x + 2.*U2x) + 
     -         ff**2*g*(U2x + 2.*U3x - 
     -            2.*ff*r*(U1x + U2x + U3x) + 
     -            ff**2*r*(U1xx + U2x**2 + U2xx + U2x*U3x + 
     -               U3x**2 + U3xx))) - 
     -      1.*H*((2.*dExp(2.*U3)*r*
     -            (derH*sn*
     -               (r**2 + 
     -                 g*(3. + ff**2*r*(U1x + 2.*U2x + 3.*U3x)))
     -               + r*(2.*cs*(2.*U2y + U3y) + 
     -                 sn*(der2H*g + 
     -                    2.*
     -                     (-1. + U1yy + U2y**2 + U2yy + 
     -                      U2y*U3y + U3y**2 + U3yy)))))/sn - 
     -         1.*dExp(2.*U2)*ff**4*g**3*sn**2*
     -          U5x**2) + dExp(2.*U2)*g**2*sn**2*
     -       U5y**2))/
     -  (dExp(2.*(U1 + U3))*g*H*r**2)

	eq11= (0.5*(dExp(2.*U3)*derH*g*H*r*sn*
     -       (3.*r**2 + g*(-3. + 
     -            ff**2*r*(2.*U1x - 1.*U2x - 3.*U3x))) + 
     -      dExp(2.*U3)*H*r*
     -       (-4.*H*r**3*sn - 
     -         1.*g**2*sn*(der2H*r + 
     -            2.*ff**2*H*
     -             (2.*U3x - 2.*ff*r*(U2x + U3x) - 
     -               2.*U1x*(1. + ff**2*r*(U2x + U3x)) + 
     -               ff**2*r*
     -                (U2x**2 + U2xx + U3x**2 + U3xx + 
     -                  4.*alpha**2*U4x**2))) + 
     -         2.*g*r*(2.*cs*(-1.*U1y + U2y) + 
     -            sn*(-1. + U2y**2 + U2yy + 
     -               H*(3. - 1.*ff**2*r*(U2x - 3.*U3x)) + 
     -               U3y**2 - 2.*U1y*(U2y + U3y) + U3yy + 
     -               4.*alpha**2*U4y**2))) + 
     -      dExp(2.*U2)*g**3*sn**3*
     -       (ff**4*g*H*U5x**2 - 1.*U5y**2)))/
     -  (dExp(2.*(U1 + U3))*g**2*H*r**2*sn)
	
	  
	eq12= (0.5*(2.*dExp(2.*U3)*cs*ff**2*g*H*r**2*
     -       (U1x - 1.*U2x) + 
     -      sn*(dExp(2.*U3)*r*
     -          (4.*H*r**2*U3y + 
     -            g*(derH*r*(U1y - 1.*U3y) + 
     -               2.*H*(U1y*(1. + ff**2*r*(U2x + U3x)) - 
     -                  1.*U3y - 
     -                  1.*ff**2*r*
     -                   (U2xy + U2x*U2y + U3xy + U3x*U3y - 
     -                     1.*U1x*(U2y + U3y) + 
     -                     4.*alpha**2*U4x*U4y)))) + 
     -         dExp(2.*U2)*ff**2*g**3*sn**2*U5x*
     -          U5y)))/
     -  (dExp(2.*(U1 + U3))*g**2*H*r**2*sn)

	Kr =  (0.25*(4.*dExp(4.*U3)*H**2*r**2*
     -       (56.*H**2*r**6*sn**2 + 
     -         g**4*sn**2*(der2H**2*r**2 - 
     -            4.*der2H*ff**2*H*r*
     -             (U1x - 2.*U3x + 2.*ff*r*U3x + 
     -               ff**2*r*U1x*U3x - 1.*ff**2*r*U3x**2 - 
     -               1.*ff**2*r*U3xx) + 
     -            4.*ff**4*H**2*
     -             (U2x**2 + 4.*U3x**2 - 8.*ff*r*U3x**2 + 
     -               2.*U1x**2*
     -                (1. + 2.*ff**2*r*(r + U3x) + 
     -                  ff**4*r**2*(U2x**2 + U3x**2)) - 
     -               4.*ff**3*r**2*
     -                (U2x**3 + U2x*U2xx + U3x**3 + U3x*U3xx) + 
     -               ff**4*r**2*
     -                (U1xx**2 + U2x**4 + U2xx**2 + U3x**4 + 
     -                  U2x**2*(2.*U2xx + U3x**2) + 
     -                  2.*U3x**2*U3xx + U3xx**2) - 
     -               2.*U1x*
     -                (2.*U3x - 2.*ff*r*U3x + 
     -                  2.*ff**3*r**2*
     -                   (U1xx - 1.*U2x**2 - 1.*U3x**2) + 
     -                  ff**2*r*(3.*U3x**2 + U3xx) + 
     -                  ff**4*r**2*
     -                   (U2x**3 + U2x*U2xx + U3x**3 + U3x*U3xx))
     -                 + 2.*ff**2*r*
     -                (2.*r*(U2x**2 + U3x**2) + 
     -                  U3x*(U2x**2 + 2.*(U3x**2 + U3xx))))) + 
     -         8.*g*H*r**4*sn*
     -          (-1.*cs*(2.*U1y - 2.*U2y + U3y) + 
     -            sn*(-1. + 2.*U1y**2 - 1.*U1yy - 2.*U1y*U2y + 
     -               U2y**2 + U2yy + 
     -               H*(-13. + ff**2*r*(5.*U1x - 8.*U3x)) - 
     -               1.*U2y*U3y + 3.*U3y**2 - 1.*U3yy)) + 
     -         4.*g**3*(der2H*r**2*sn**2*
     -             (H*(-3. + ff**2*r*(U1x - 2.*U3x)) + U1y*U3y)
     -             + 2.*H*(-2.*ff**5*H*r**3*sn**2*
     -                (U1x**2 - 1.*U1x*U2x + 2.*U2x**2 + 
     -                  U1x*U3x - 2.*U3x**2) + 
     -               ff**6*H*r**3*sn**2*
     -                (2.*U2x**3 + 2.*U1x**2*(U2x - 1.*U3x) - 
     -                  1.*U2x**2*U3x + U2x*(2.*U2xx + U3x**2) - 
     -                  2.*U3x*(U3x**2 + U3xx) + 
     -                  U1x*
     -                   (U1xx - 2.*U2x**2 - 1.*U2xx + 
     -                     4.*U3x**2 + U3xx)) + 
     -               sn**2*(U1y - 1.*U3y)**2 - 
     -               2.*ff**3*r**2*sn*
     -                (cs*U1y*U2x + H*sn*(U1x + U2x - 3.*U3x) + 
     -                  sn*(U1x*U1yy + U1y*U2x*U2y + U1y*U3x*U3y)
     -                  ) + 
     -               ff**2*r*sn*
     -                (H*sn*(4.*U1x + U2x - 6.*U3x) + 
     -                  cs*U2x*U3y + 
     -                  sn*
     -                   (2.*U1y**2*U3x + U2x*U2y*U3y + 
     -                     2.*U3xy*U3y - 1.*U1x*U3y**2 + 
     -                     2.*U3x*U3y**2 - 
     -                     2.*U1y*(U3xy + U3x*U3y) + U1x*U3yy))
     -                + ff**4*r**2*
     -                (cs**2*(U1x - 1.*U2x)**2 + 
     -                  cs*sn*
     -                   (U1y*(-1.*U2x**2 + U2xx) + 
     -                     2.*U1x**2*U2y - 
     -                     2.*U1x*(U2xy + U2x*U2y) + 
     -                     U2x*(2.*U2xy + 2.*U2x*U2y + U3x*U3y))
     -                   + sn**2*
     -                   (U1xx*U1yy - 1.*U1x*U2x + 
     -                     U1y**2*U2x**2 - 2.*U1y*U2x*U2xy + 
     -                     U2xy**2 - 1.*U1y*U2x**2*U2y + 
     -                     U1y*U2xx*U2y - 2.*U1x*U2xy*U2y + 
     -                     2.*U2x*U2xy*U2y + U1x**2*U2y**2 - 
     -                     1.*U1x*U2x*U2y**2 + U2x**2*U2y**2 + 
     -                     U1x*U2x*U2yy + U1y**2*U3x**2 + 
     -                     H*
     -                      (-2.*U1x**2 + U1xx - 1.*U1x*U2x + 
     -                      U2xx + 9.*U1x*U3x + 2.*U2x*U3x - 
     -                      7.*U3x**2 - 3.*U3xx) - 
     -                     2.*U1y*U3x*U3xy + U3xy**2 + 
     -                     U2x*U2y*U3x*U3y - 1.*U1y*U3x**2*U3y + 
     -                     U1y*U3xx*U3y - 2.*U1x*U3xy*U3y + 
     -                     2.*U3x*U3xy*U3y + U1x**2*U3y**2 - 
     -                     1.*U1x*U3x*U3y**2 + U3x**2*U3y**2 + 
     -                     U1x*U3x*U3yy)))) + 
     -         4.*g**2*r**2*
     -          (cs**2*(2.*U1y**2 - 4.*U1y*U2y + 4.*U2y**2 + 
     -               U3y**2) + 
     -            2.*cs*sn*
     -             (2.*U1y**2*U2y + 2.*U2y**3 - 
     -               1.*U1y*
     -                (-1. + H*(-1. + ff**2*r*U2x) + 3.*U2y**2 + 
     -                  U2yy) + 
     -               H*(1. + ff**2*r*(-1.*U2x + U3x))*U3y + 
     -               U2y*(-2. + 2.*ff**2*H*r*(U1x + U2x) + 
     -                  2.*U2yy + U3y**2)) + 
     -            sn**2*(1. + U1yy**2 + 2.*U1y*U2y - 2.*U2y**2 + 
     -               2.*U1y**2*U2y**2 - 2.*U1y*U2y**3 + U2y**4 - 
     -               2.*U2yy - 2.*U1y*U2y*U2yy + 
     -               2.*U2y**2*U2yy + U2yy**2 + 
     -               H**2*(13. - 4.*ff**2*r*(4.*U1x - 7.*U3x) + 
     -                  4.*ff**3*r**2*(U1x + U2x - 3.*U3x) + 
     -                  ff**4*r**2*
     -                   (5.*U1x**2 + 2.*U1x*(U2x - 7.*U3x) - 
     -                     2.*
     -                      (U1xx - 2.*U2x**2 + U2xx + 
     -                      2.*U2x*U3x - 6.*U3x**2 - 3.*U3xx)))
     -                + 2.*U1y**2*U3y**2 + U2y**2*U3y**2 - 
     -               2.*U1y*U3y**3 + U3y**4 - 2.*U1y*U3y*U3yy + 
     -               2.*U3y**2*U3yy + U3yy**2 + 
     -               H*(3.*der2H*r**2 + 
     -                  2.*
     -                   (-2.*U1y**2 + U1yy + U1y*U2y + 
     -                     2.*U1y*U3y + U2y*U3y - 3.*U3y**2 + 
     -                     U3yy + 
     -                     ff**2*r*
     -                      (-2.*U1y*U2xy - 2.*U1y**2*U3x + 
     -                      2.*U1y*U3xy + 3.*U1y*U3x*U3y + 
     -                      U2y*U3x*U3y - 4.*U3xy*U3y - 
     -                      3.*U3x*U3y**2 + 
     -                      U2x*
     -                      (-1. + 2.*U1y**2 - 1.*U1y*U2y + 
     -                      U2y**2 + U2yy - 1.*U2y*U3y) + 
     -                      U1x*
     -                      (-1. + U1yy + U2y**2 + U2yy + 
     -                      3.*U3y**2 - 1.*U3yy) + U3x*U3yy))))))
     -        + 11.*dExp(4.*U2)*g**6*sn**6*
     -       (ff**4*g*H*U5x**2 + U5y**2)**2 + 
     -      4.*dExp(2.*U3)*derH**2*g**2*r**2*
     -       sn**2*(dExp(2.*U3)*H*
     -          (13.*H*r**4 + 
     -            g**2*H*(9. - 6.*ff**2*r*(U1x - 3.*U3x) + 
     -               ff**4*r**2*
     -                (3.*U1x**2 + 2.*U2x**2 - 6.*U1x*U3x + 
     -                  9.*U3x**2)) + 
     -            2.*g*r**2*
     -             (H*(-9. + 
     -                  ff**2*r*(5.*U1x + 2.*U2x - 9.*U3x)) + 
     -               (U1y - 1.*U3y)**2)) - 
     -         1.*dExp(2.*U2)*g**3*sn**2*U5y**2)
     -       + 4.*dExp(2.*U3)*derH*g*H*r*sn*
     -       (-2.*dExp(2.*U3)*der2H*g**2*H*r**2*
     -          sn*(3.*r**2 + g*(-3. + ff**2*r*(U1x - 3.*U3x)))
     -          + 4.*dExp(2.*U3)*H*r*
     -          (-13.*H*r**5*sn + 
     -            ff**2*g**3*H*sn*
     -             (6.*U3x - 6.*ff*r*U3x + 
     -               2.*ff**2*r*U1x**2*
     -                (1. - 1.*ff*r + ff**2*r*U3x) - 
     -               2.*ff**3*r**2*(U2x**2 + 3.*U3x**2) + 
     -               U1x*(-3. - 8.*ff**2*r*U3x + 
     -                  2.*ff**3*r**2*U3x + 
     -                  ff**4*r**2*
     -                   (U1xx - 1.*U2x**2 - 4.*U3x**2 - 1.*U3xx)
     -                  ) + 
     -               ff**2*r*(U2x**2 + 9.*U3x**2 + 3.*U3xx) + 
     -               ff**4*r**2*
     -                (U2x**3 + U2x*U2xx + U2x**2*U3x + 
     -                  3.*U3x*(U3x**2 + U3xx))) + 
     -            g*r**3*(H*sn*
     -                (22. - 1.*ff**2*r*(9.*U1x + U2x - 17.*U3x))
     -                 + cs*(U1y + U3y) + 
     -               sn*(-2.*U1y**2 + U1yy + U1y*U2y + 
     -                  2.*U1y*U3y + U2y*U3y - 3.*U3y**2 + U3yy))
     -              + g**2*r*
     -             (H*sn*(-9. + 
     -                  3.*ff**2*r*(4.*U1x + U2x - 7.*U3x) - 
     -                  2.*ff**3*r**2*(U1x + U2x - 3.*U3x) + 
     -                  ff**4*r**2*
     -                   (-1.*U1x**2 + U1xx + 2.*U2x**2 + U2xx - 
     -                     2.*U1x*(U2x - 5.*U3x) + 2.*U2x*U3x - 
     -                     9.*U3x**2 - 3.*U3xx)) + 
     -               cs*ff**2*r*U2x*(U1y + U3y) + 
     -               sn*(2.*U1y**2*(1. + ff**2*r*U3x) + 
     -                  2.*U3y**2 + 
     -                  U1y*
     -                   (-1.*U3y + 
     -                     ff**2*r*
     -                      (U2x*U2y - 2.*U3xy - 1.*U3x*U3y)) + 
     -                  ff**2*r*
     -                   (U3y*(U2x*U2y + 2.*U3xy + 2.*U3x*U3y) + 
     -                     U1x*(U1yy - 1.*U3y**2 + U3yy))))) + 
     -         dExp(2.*U2)*g**3*sn**3*
     -          (ff**4*g*H*(-9.*g + 11.*r**2)*U5x**2 + 
     -            ff**6*g**2*H*r*(3.*U1x + 2.*U2x - 9.*U3x)*
     -             U5x**2 + 2.*(-2.*g + 3.*r**2)*U5y**2 + 
     -            ff**2*g*r*U5y*
     -             (4.*U1y*U5x - 4.*U3y*U5x + 2.*U5xy - 
     -               5.*U1x*U5y + 7.*U2x*U5y - 4.*U3x*U5y))) - 
     -      4.*dExp(2.*(U2 + U3))*g**3*H*sn**2*
     -       (6.*cs**2*g*r**2*
     -          (2.*ff**4*g*H*U5x**2 + 3.*U5y**2) - 
     -         2.*cs*g*r**2*sn*
     -          (-6.*ff**2*H*r*U5x*U5y + 
     -            ff**4*g*H*U5x*
     -             (7.*U1y*U5x - 12.*U2y*U5x + 7.*U3y*U5x - 
     -               6.*U5xy - 2.*U1x*U5y - 4.*U2x*U5y) + 
     -            U5y*(5.*U1y*U5y - 16.*U2y*U5y + 7.*U3y*U5y - 
     -               6.*U5yy)) + 
     -         sn**2*(8.*ff**5*g**2*H**2*r*(g - 4.*r**2)*
     -             U5x**2 + 
     -            4.*ff**7*g**3*H**2*r**2*U5x*
     -             (2.*U1x*U5x - 5.*U2x*U5x - 1.*U3x*U5x - 
     -               2.*U5xx) + 
     -            2.*ff**8*g**3*H**2*r**2*
     -             (2.*U1x**2*U5x**2 + 8.*U2x**2*U5x**2 - 
     -               1.*U2xx*U5x**2 + 4.*U3x**2*U5x**2 + 
     -               3.*U3xx*U5x**2 - 2.*U3x*U5x*U5xx + 
     -               U5xx**2 - 
     -               1.*U1x*U5x*
     -                (5.*U2x*U5x + U3x*U5x + 2.*U5xx) + 
     -               U2x*U5x*(-7.*U3x*U5x + 6.*U5xx)) + 
     -            2.*ff**6*g**2*H**2*r*U5x*
     -             (g*(4.*r*U5x - 1.*U1x*U5x - 7.*U2x*U5x + 
     -                  8.*U3x*U5x - 2.*U5xx) + 
     -               r**2*(-2.*U1x*U5x + 23.*U2x*U5x - 
     -                  15.*U3x*U5x + 8.*U5xx)) - 
     -            8.*ff**3*g**2*H*r**2*U1y*U5x*U5y - 
     -            2.*ff**2*g*H*r*
     -             (g*U5y*(6.*U1y*U5x - 4.*U3y*U5x + 2.*U5xy - 
     -                  5.*U1x*U5y + 7.*U2x*U5y - 4.*U3x*U5y) - 
     -               2.*r**2*
     -                (5.*U1y*U5x*U5y + 3.*U2y*U5x*U5y - 
     -                  5.*U3y*U5x*U5y + 2.*U5xy*U5y - 
     -                  4.*U1x*U5y**2 + 6.*U2x*U5y**2 - 
     -                  3.*U3x*U5y**2 + U5x*U5yy)) + 
     -            ff**4*g*H*
     -             (56.*H*r**4*U5x**2 + 
     -               g**2*(2.*H + 3.*der2H*r**2)*U5x**2 - 
     -               2.*g*r**2*
     -                (19.*H*U5x**2 - 2.*U1y**2*U5x**2 + 
     -                  7.*U1y*U2y*U5x**2 - 6.*U2y**2*U5x**2 - 
     -                  5.*U1y*U3y*U5x**2 + 7.*U2y*U3y*U5x**2 - 
     -                  2.*U3y**2*U5x**2 + 4.*U1y*U5x*U5xy - 
     -                  6.*U2y*U5x*U5xy + 2.*U3y*U5x*U5xy - 
     -                  2.*U5xy**2 - 2.*U1y*U2x*U5x*U5y + 
     -                  2.*U2xy*U5x*U5y - 2.*U1x*U2y*U5x*U5y - 
     -                  4.*U2x*U2y*U5x*U5y + 
     -                  6.*U1y*U3x*U5x*U5y - 6.*U3xy*U5x*U5y + 
     -                  6.*U1x*U3y*U5x*U5y - 
     -                  4.*U3x*U3y*U5x*U5y - 2.*U1y*U5xx*U5y + 
     -                  4.*U1x*U5xy*U5y - 6.*U2x*U5xy*U5y + 
     -                  2.*U3x*U5xy*U5y - 2.*U1x**2*U5y**2 + 
     -                  7.*U1x*U2x*U5y**2 - 6.*U2x**2*U5y**2 - 
     -                  5.*U1x*U3x*U5y**2 + 7.*U2x*U3x*U5y**2 - 
     -                  2.*U3x**2*U5y**2 - 2.*U1x*U5x*U5yy)) + 
     -            2.*(2.*g**2*H*U5y**2 + 5.*H*r**4*U5y**2 + 
     -               g*r**2*
     -                (-1.*
     -                   (-1. + 6.*H - 2.*U1y**2 + 5.*U1y*U2y - 
     -                     8.*U2y**2 + U2yy + U1y*U3y + 
     -                     7.*U2y*U3y - 4.*U3y**2 - 3.*U3yy)*
     -                   U5y**2 - 
     -                  2.*(U1y - 3.*U2y + U3y)*U5y*U5yy + 
     -                  U5yy**2))))))/
     -  (dExp(4.*(U1 + U3))*g**4*H**2*r**4*sn**2)

           WRITE(12,489)Ricci,eq11,eq12,Kr

 11        CONTINUE

	        Close(Unit=1)
	        Close(Unit=12)


c            write(6,*)" "
c     	        write(6,*)"max err1 at=", r1,t1
c  	     	write(6,*)"max err2 at=", r2,t2 
   	  write(6,*)" " 
   
c          WRITE(6,*)lambda,alpha

c	        open(unit=5,File='res.txt')
c	        WRITE(5,765)  eta,nr,w,g,lambda1,lambda2,eta1,const 
c	        Close(Unit=5)
 
	   	   	write(6,*)"Max Z... at (r) =", max,rm
  	  write(6,*)" "

        open(unit=1,name='derx.dat')
         rewind(1)
        Close(Unit=1) 

         rewind(1)
        open(unit=1,name='derx.dat')
      DO 110  I = 1,NV
          WRITE(1,789)  (1.-XX(I))**2*UX(I,1) 
          WRITE(1,789)  (1.-XX(I))**2*UX(I,2) 
          WRITE(1,789)  (1.-XX(I))**2*UX(I,3) 
	    WRITE(1,789)  (1.-XX(I))**2*UX(I,4) 
	    WRITE(1,789)  (1.-XX(I))**2*UX(I,5) 	  		   		    
  110   CONTINUE
c        Close(Unit=1)
                


        open(unit=1,name='dery.dat')
         rewind(1)
        Close(Unit=1)
                                    
         rewind(1)
        open(unit=1,name='dery.dat')                      
      DO 111  I = 1,NV
        WRITE(1,789)  UY(I,1)
        WRITE(1,789)  UY(I,2) 
        WRITE(1,789)  UY(I,3)
        WRITE(1,789)  UY(I,4)
        WRITE(1,789)  UY(I,5) 
  111   CONTINUE  


         write(6,*)" "
       do 1200 i=1,nk
        ptest(i)=0.d0
 	      do 1300 j=1,nv
 	        ptest(i)=ptest(i)+p(j,i)**2
 1300   continue 
        write(6,*)"res=",i, sqrt(ptest(i))
 1200 continue 

           write(6,*)" "

c  	write(6,*)"here"   

  789     format (e24.16)
  489     format (6(e24.16))
  765     format (8(e24.16))  
c         write(6,*)"HERE@#"                     
      return  
      end   
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CX
C X
C  X
C   X
C    X
C     X
C      X
C       X
C        X
C         X
C          X
C           X
C            X
C             X
C              X
C               X
C                X
C                 X
C                  X
C                   X
C                    X
C                     X
C                      X
C                       X
C                        X
C                         X
C                          X
C                           X
C                            X
C                             X
C                              X
C                               X
C                                X
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
      SUBROUTINE CDSU02(IRAND,T,XX,Y,U,UT,UX,UY,UXX,UXY,UYY,
     *                  P,MT,MV,NK,NB)
      implicit character*1 (a-z)
      INTEGER  IRAND, MT, MV, NK, NB
      double precision  T,XX(MV),Y(MV),U(MV,NK),UT(MT,NK),
     &       UX(MV,NK),UY(MV,NK),UXX(MV,NK),UXY(MV,NK),
     &       UYY(MV,NK),P(MV,NK)
      double precision eng 
      INTEGER  I,J,k
      logical fifo_exist
 
      integer m_nx,m_ny
      double precision  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
 
      common   nr,w,con,alpha,c1,c2,c3,lambda,rh,par

      common/derivs/ u0x(1505,8),u0y(1505,8),u0xx(1505,8),u0yy(1505,8),
     *     u0xy(1505,8),u0(1505,8)
 
      common/dimensions/m_nx,m_ny

      common /parae/ eng
c      DO 5 J=1,NK
c      DO 5 I=1,NB
c         P(I,J) = 0.
c    5 CONTINUE

        
			rewind(1)
			rewind(2)
			rewind(3)
			rewind(4)
			rewind(5)
			rewind(26)
			rewind(7)
			rewind(8)
			rewind(9)
			rewind(10)
			rewind(11)
			rewind(12)
			rewind(13)
			rewind(14)
			rewind(15)
			rewind(17)
			rewind(18)
			rewind(19)
			rewind(20)
			rewind(21) 

	   Open(Unit=1,File='f-0.txt') 
	   Open(Unit=2,File='fx-0.txt')
	   Open(Unit=3,File='fxx-0.txt') 
	   Open(Unit=4,File='fy-0.txt') 
	   Open(Unit=5,File='fyy-0.txt')
c
	   Open(Unit=26,File='f-inf.txt') 
	   Open(Unit=7,File='fx-inf.txt') 
	   Open(Unit=8,File='fxx-inf.txt') 
	   Open(Unit=9,File='fy-inf.txt') 
	   Open(Unit=10,File='fyy-inf.txt')
c
	   Open(Unit=11,File='f-t0.txt')
	   Open(Unit=12,File='fx-t0.txt')
	   Open(Unit=13,File='fxx-t0.txt')
	   Open(Unit=14,File='fy-t0.txt')
	   Open(Unit=15,File='fyy-t0.txt')
c
c
	   Open(Unit=16,File='f-tPi2.txt')
	   Open(Unit=17,File='fx-tPi2.txt')
	   Open(Unit=18,File='fxx-tPi2.txt')
	   Open(Unit=19,File='fy-tPi2.txt')
	   Open(Unit=20,File='fyy-tPi2.txt')


 
      GO TO (10,20,30,40) IRAND

c  U1 = f1
c  U2 = f2
c  U3 = f0
c  U4 = Z
c  U5 = W
 
 
c I treat separate the i=1 case ( u(1,4) is registered in the wrong way)

   10 continue

          p(1,1)=Ux(1,1) 
  	    p(1,2)=Ux(1,2) 
	    p(1,3)=ux(1,3)  
	    p(1,4)=ux(1,4)  
	    p(1,5)=u(1,5) +w/nr

      write(1,765)y(1),u(1,1),u(1,2),u(1,3),0.d0,u(1,5)  
      write(2,765)y(1),ux(1,1),ux(1,2),ux(1,3),ux(1,4),ux(1,5)  
      write(3,765)y(1),uxx(1,1),uxx(1,2),uxx(1,3),uxx(1,4),uxx(1,5) 
      write(4,765)y(1),uy(1,1),uy(1,2),uy(1,3),uy(1,4),uy(1,5) 
      write(5,765)y(1),uyy(1,1),uyy(1,2),uyy(1,3),uyy(1,4),uyy(1,5) 

      
       DO 11  I = 2,NB 
          p(i,1)=Ux(I,1) 
  	    p(i,2)=Ux(I,2) 
	    p(i,3)=ux(i,3)  
	    p(i,4)=ux(i,4)  
	    p(i,5)=u(i,5) +w/nr

c      	   write(6,*)U(I,4) 

      write(1,765)y(I),u(i,1),u(i,2),u(i,3),u(i,4),u(i,5)  
      write(2,765)y(I),ux(i,1),ux(i,2),ux(i,3),ux(i,4),ux(i,5)  
      write(3,765)y(I),uxx(i,1),uxx(i,2),uxx(i,3),uxx(i,4),uxx(i,5) 
      write(4,765)y(I),uy(i,1),uy(i,2),uy(i,3),uy(i,4),uy(i,5) 
      write(5,765)y(I),uyy(i,1),uyy(i,2),uyy(i,3),uyy(i,4),uyy(i,5) 


   11  CONTINUE

                 Close(Unit=1)
                 Close(Unit=2)
                 Close(Unit=3)
                 Close(Unit=4)
                 Close(Unit=5)
 
      GO TO 50


c       conditii la infinit 
   20 DO 26  I = 1,NB         
          p(i,1)=u(i,1) 
	    p(i,2)=u(i,2) 
	    p(i,3)=u(i,3) 
	    p(i,4)=u(i,4)  
	    p(i,5)=u(i,5) 

      write(26,765)y(I),u(i,1),u(i,2),u(i,3),u(i,4),u(i,5)     
      write(7,765)y(I),ux(i,1),ux(i,2),ux(i,3),0.d0,ux(i,5)   
      write(8,765)y(I),uxx(i,1),uxx(i,2),uxx(i,3),uxx(i,4),uxx(i,5)  
      write(9,765)y(I),uy(i,1),uy(i,2),uy(i,3),uy(i,4),uy(i,5)  
      write(10,765)y(I),uyy(i,1),uyy(i,2),uyy(i,3),uyy(i,4),uyy(i,5)  

   26   CONTINUE

   
                 Close(Unit=26)
                 Close(Unit=7)
                 Close(Unit=8)
                 Close(Unit=9)
                 Close(Unit=10)
      GO TO 50 

c       conditii theta=0
   30   DO 37  i = 1,NB
c          p(i,1)=u(i,1)-U(I,2) 
          p(i,1)=uy(i,1)
	    p(i,2)=uy(i,2) 
	    p(i,3)=uy(i,3) 
	    p(i,4)=u(i,4)  
	    p(i,5)=uy(i,5) 

 
      write(11,765)y(I),u(i,1),u(i,2),u(i,3),0.d0,u(i,5)  
      write(12,765)y(I),ux(i,1),ux(i,2),ux(i,3),ux(i,4),ux(i,5)  
      write(13,765)y(I),uxx(i,1),uxx(i,2),uxx(i,3),uxx(i,4),uxx(i,5) 
      write(14,765)y(I),uy(i,1),uy(i,2),uy(i,3),uy(i,4),uy(i,5) 
      write(15,765)y(I),uyy(i,1),uyy(i,2),uyy(i,3),uyy(i,4),uyy(i,5) 
   37   CONTINUE

                 Close(Unit=11)
                 Close(Unit=12)
                 Close(Unit=13)
                 Close(Unit=14)
                 Close(Unit=15)
      GO TO 50

c       conditii theta=Pi/2 
   40   DO 47  i = 1,NB      
          p(i,1)=uy(i,1) 
	    p(i,2)=uy(i,2) 
	    p(i,3)=uy(i,3) 
	    p(i,4)=uy(i,4) 
	    p(i,5)=uy(i,5)

      write(16,765)y(I),u(i,1),u(i,2),u(i,3),u(i,4),u(i,5)  
      write(17,765)y(I),ux(i,1),ux(i,2),ux(i,3),ux(i,4),ux(i,5)  
      write(18,765)y(I),uxx(i,1),uxx(i,2),uxx(i,3),uxx(i,4),uxx(i,5) 
      write(19,765)y(I),uy(i,1),uy(i,2),uy(i,3),uy(i,4),uy(i,5) 
      write(20,765)y(I),uyy(i,1),uyy(i,2),uyy(i,3),uyy(i,4),uyy(i,5)  
	    
 47   CONTINUE

                 Close(Unit=16)
                 Close(Unit=17)
                 Close(Unit=18)
                 Close(Unit=19)
                 Close(Unit=20)
      GO TO 50 
 

 50   CONTINUE 

 765     format (9(e24.16)) 
      R E T U R N
 111  FORMAT(44(e24.16,1x))
c 112  FORMAT('ux:',13(e24.16,1x))
c 113  FORMAT('uy:',13(e24.16,1x))
 114  FORMAT('uxx:',13(e24.16,1x))
 115  FORMAT('uxy:',13(e24.16,1x))
 116  FORMAT('uyy:',13(e24.16,1x))
      END


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CX
C X
C  X
C   X
C    X
C     X
C      X
C       X
C        X
C         X
C          X
C           X
C            X
C             X
C              X
C               X
C                X
C                 X
C                  X
C                   X
C                    X
C                     X
C                      X
C                       X
C                        X
C                         X
C                          X
C                           X
C                            X
C                             X
C                              X
C                               X
C                                X
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      SUBROUTINE CDSU04(IRAND,IEQU,ICOM,T,X,Y,U,UT,UX,UY,UXX,UXY,UYY,
     *                  PU,PUT,PUX,PUY,PUXX,PUXY,PUYY,MT,MV,NK,NB)
c      implicit double precision (a-h,o-z)
      implicit character*1 (a-z)
      INTEGER  IRAND, IEQU, ICOM, MT, MV, NK, NB
      double precision T,X(MV),Y(MV),U(MV,NK),UT(MT,NK),UX(MV,NK),
     *                  UY(MV,NK),UXX(MV,NK),UXY(MV,NK),UYY(MV,NK),
     *                  PU(MV),PUT(MT),PUX(MV),PUY(MV),PUXX(MV),
     *                  PUXY(MV),PUYY(MV)
	double precision rH,OmegaH,alfa,mu
	 common  rH,OmegaH,alfa,mu


      INTEGER  I,J
          DO 50 J=1,Nb
      PU(J)=0.
      PUX(J)=0.
      PUY(J)=0.
      PUXX(J)=0. 
      PUYY(J)=0.
   50 continue
  
      GO TO (1000,2000,3000,4000) IRAND

c boundary conditions r=0
 1000 CONTINUE
      DO 1111 i=1,Nb

c         p(i,1)=U(I,1)+rh/(rh+par)/2.*Ux(I,1)
c	    p(i,2)=U(I,2)+rh/(rh+par)/2.*Ux(I,2)  
c	    p(i,3)=u(i,3) 
c	    p(i,4)=u(i,4)-OmegaH 
c	    p(i,5)=ux(i,5) 

c          write(6,*)rh,par

 	 IF(IEQU.EQ.1.AND.ICOM.EQ.1) pux(i)=1. 
 	 IF(IEQU.EQ.2.AND.ICOM.EQ.2) pux(i)=1.    
       IF(IEQU.EQ.3.AND.ICOM.EQ.3) pux(i)=1.
       IF(IEQU.EQ.4.AND.ICOM.EQ.4) pux(i)=1.  
       IF(IEQU.EQ.5.AND.ICOM.EQ.5) pu(i)=1. 

 1111 continue
      GOTO 7000
 

c boundary conditions r=infinity
 2000 CONTINUE
      DO 2100 I=1,NB 
	 IF(IEQU.EQ.1.AND.ICOM.EQ.1) pu(i)=1.
       IF(IEQU.EQ.2.AND.ICOM.EQ.2) pu(i)=1. 
       IF(IEQU.EQ.3.AND.ICOM.EQ.3) pu(i)=1.
       IF(IEQU.EQ.4.AND.ICOM.EQ.4) pu(i)=1.
       IF(IEQU.EQ.5.AND.ICOM.EQ.5) pu(i)=1. 
 2100 CONTINUE
      GO TO 7000


c boundary conditions theta=0
 3000 CONTINUE
      DO 3100 I=1,NB
c         p(i,1)=uy(i,1) 
c	    p(i,2)=uy(i,2) 
c	    p(i,3)=uy(i,3) 
c	    p(i,4)=u(i,4) 
c	    p(i,5)=uy(i,5) 
	 IF(IEQU.EQ.1.AND.ICOM.EQ.1) puy(i)=1.
       IF(IEQU.EQ.2.AND.ICOM.EQ.2) puy(i)=1. 
       IF(IEQU.EQ.3.AND.ICOM.EQ.3) puy(i)=1.
       IF(IEQU.EQ.4.AND.ICOM.EQ.4) pu(i)=1.
       IF(IEQU.EQ.5.AND.ICOM.EQ.5) puy(i)=1. 
 3100 CONTINUE
      GO TO 7000
 

c boundary conditions theta=Pi/2
 4000 CONTINUE
      DO 4100 I=1,NB
c         p(i,1)=uy(i,1) 
c	    p(i,2)=uy(i,2) 
c	    p(i,3)=uy(i,3) 
c	    p(i,4)=uy(i,4) 
c	    p(i,5)=uy(i,5)  
	 IF(IEQU.EQ.1.AND.ICOM.EQ.1) puy(i)=1.
       IF(IEQU.EQ.2.AND.ICOM.EQ.2) puy(i)=1. 
       IF(IEQU.EQ.3.AND.ICOM.EQ.3) puy(i)=1.
       IF(IEQU.EQ.4.AND.ICOM.EQ.4) puy(i)=1.
       IF(IEQU.EQ.5.AND.ICOM.EQ.5) puy(i)=1. 

 4100 CONTINUE
      GO TO 7000 
 7000 CONTINUE        
 
      RETURN
      END




c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CX
C X
C  X
C   X
C    X
C     X
C      X
C       X
C        X
C         X
C          X
C           X
C            X
C             X
C              X
C               X
C                X
C                 X
C                  X
C                   X
C                    X
C                     X
C                      X
C                       X
C                        X
C                         X
C                          X
C                           X
C                            X
C                             X
C                              X
C                               X
C                                X
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      subroutine cdsu03(iequ,icom,t,xx,y,u,ut,ux,uy,uxx,uxy,uyy,
     *                  pu,put,pux,puy,puxx,puxy,puyy,
     *                  mt,mv,nk,nv)
      implicit character*1 (a-z)
      integer  iequ, icom, mt, mv, nk, nv,I
      double precision  t,xx(mv),y(mv),u(mv,nk),ut(mt,nk),ux(mv,nk),
     *                  uy(mv,nk),uxx(mv,nk),uxy(mv,nk),uyy(mv,nk),
     *                  pu(mv),puxy(mv),
     *                  put(mt),pux(mv),puxx(mv),puy(mv),puyy(mv)  
  
	Double precision 
     & U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,
     & U1x,U2x,U3x,U4x,U5x,U6x,U7x,U8x,U9x,U10x,U11x,
     & U1y,U2y,U3y,U4y,U5y,U6y,U7y,U8y,U9y,U10y,U11y,
     & U1xy,U2xy,U3xy,U4xy,U5xy,U6xy,U7xy,U8xy,U9xy,U10xy,U11xy,
     & U1xx,U2xx,U3xx,U4xx,U5xx,U6xx,U7xx,U8xx,U9xx,U10xx,U11xx,
     & U1yy,U2yy,U3yy,U4yy,U5yy,U6yy,U7yy,U8yy,U9yy,U10yy,U11yy 
 

	 Double precision 
     & Vxx11,Vxx12,Vxx13,Vxx14,Vxx15,Vxx16,Vxx17,Vxx18,Vxx19,
     & Vxx110,Vxx111,Vxx112,
     & Vxx21,Vxx22,Vxx23,Vxx24,Vxx25,Vxx26,Vxx27,Vxx28,Vxx29,
     & Vxx210,Vxx211,Vxx212,
     & Vxx31,Vxx32,Vxx33,Vxx34,Vxx35,Vxx36,Vxx37,Vxx38,Vxx39,
     & Vxx310,Vxx311,Vxx312,
     & Vxx41,Vxx42,Vxx43,Vxx44,Vxx45,Vxx46,Vxx47,Vxx48,Vxx49,
     & Vxx410,Vxx411,Vxx412,
     & Vxx51,Vxx52,Vxx53,Vxx54,Vxx55,Vxx56,Vxx57,Vxx58,Vxx59,
     & Vxx510,Vxx511,Vxx512,
     & Vxx61,Vxx62,Vxx63,Vxx64,Vxx65,Vxx66,Vxx67,Vxx68,Vxx69,
     & Vxx610,Vxx611,Vxx612,
     & Vxx71,Vxx72,Vxx73,Vxx74,Vxx75,Vxx76,Vxx77,Vxx78,Vxx79,
     & Vxx710,Vxx711,Vxx712,
     & Vxx81,Vxx82,Vxx83,Vxx84,Vxx85,Vxx86,Vxx87,Vxx88,Vxx89,
     & Vxx810,Vxx811,Vxx812,
     & Vxx91,Vxx92,Vxx93,Vxx94,Vxx95,Vxx96,Vxx97,Vxx98,Vxx99,
     & Vxx910,Vxx911,Vxx912,
     & Vxx101,Vxx102,Vxx103,Vxx104,Vxx105,Vxx106,Vxx107,
     & Vxx108,Vxx109,Vxx1010,Vxx1011,Vxx1012,
     & Vxx11n1,Vxx11n2,Vxx11n3,Vxx11n4,Vxx11n5,Vxx11n6,Vxx11n7, 
     & Vxx11n8,Vxx11n9,Vxx11n10,Vxx11n11,Vxx11n12,     
     & Vxx12n1,Vxx12n2,Vxx12n3,Vxx12n4,Vxx12n5,Vxx12n6,Vxx12n7, 
     & Vxx12n8,Vxx12n9,Vxx12n10,Vxx12n11,Vxx12n12


		 Double precision 
     & Vyy11,Vyy12,Vyy13,Vyy14,Vyy15,Vyy16,Vyy17,Vyy18,Vyy19,
     & Vyy110,Vyy111,Vyy112,
     & Vyy21,Vyy22,Vyy23,Vyy24,Vyy25,Vyy26,Vyy27,Vyy28,Vyy29,
     & Vyy210,Vyy211,Vyy212,
     & Vyy31,Vyy32,Vyy33,Vyy34,Vyy35,Vyy36,Vyy37,Vyy38,Vyy39,
     & Vyy310,Vyy311,Vyy312,
     & Vyy41,Vyy42,Vyy43,Vyy44,Vyy45,Vyy46,Vyy47,Vyy48,Vyy49,
     & Vyy410,Vyy411,Vyy412,
     & Vyy51,Vyy52,Vyy53,Vyy54,Vyy55,Vyy56,Vyy57,Vyy58,Vyy59,
     & Vyy510,Vyy511,Vyy512,
     & Vyy61,Vyy62,Vyy63,Vyy64,Vyy65,Vyy66,Vyy67,Vyy68,Vyy69,
     & Vyy610,Vyy611,Vyy612,
     & Vyy71,Vyy72,Vyy73,Vyy74,Vyy75,Vyy76,Vyy77,Vyy78,Vyy79,
     & Vyy710,Vyy711,Vyy712,
     & Vyy81,Vyy82,Vyy83,Vyy84,Vyy85,Vyy86,Vyy87,Vyy88,Vyy89,
     & Vyy810,Vyy811,Vyy812,
     & Vyy91,Vyy92,Vyy93,Vyy94,Vyy95,Vyy96,Vyy97,Vyy98,Vyy99,
     & Vyy910,Vyy911,Vyy912,
     & Vyy101,Vyy102,Vyy103,Vyy104,Vyy105,Vyy106,Vyy107,
     &                 Vyy108,Vyy109,Vyy1010,Vyy1011,Vyy1012,
     & Vyy11n1,Vyy11n2,Vyy11n3,Vyy11n4,Vyy11n5,Vyy11n6,Vyy11n7, 
     &       Vyy11n8,Vyy11n9,Vyy11n10,Vyy11n11,Vyy11n12,     
     & Vyy12n1,Vyy12n2,Vyy12n3,Vyy12n4,Vyy12n5,Vyy12n6,Vyy12n7, 
     &      Vyy12n8,Vyy12n9,Vyy12n10,Vyy12n11,Vyy12n12



		 Double precision 
     & Vxy11,Vxy12,Vxy13,Vxy14,Vxy15,Vxy16,Vxy17,Vxy18,
     &                                           Vxy19,Vxy110,Vxy111, 
     & Vxy21,Vxy22,Vxy23,Vxy24,Vxy25,Vxy26,Vxy27,Vxy28,
     &                                           Vxy29,Vxy210,Vxy211,                                        
     & Vxy31,Vxy32,Vxy33,Vxy34,Vxy35,Vxy36,Vxy37,Vxy38,
     &                                           Vxy39,Vxy310,Vxy311,   
     & Vxy41,Vxy42,Vxy43,Vxy44,Vxy45,Vxy46,Vxy47,Vxy48,
     &                                           Vxy49,Vxy410,Vxy411,    
     & Vxy51,Vxy52,Vxy53,Vxy54,Vxy55,Vxy56,Vxy57,Vxy58,
     &                                           Vxy59,Vxy510,Vxy511, 
     & Vxy61,Vxy62,Vxy63,Vxy64,Vxy65,Vxy66,Vxy67,Vxy68,
     &                                           Vxy69,Vxy610,Vxy611, 
     & Vxy71,Vxy72,Vxy73,Vxy74,Vxy75,Vxy76,Vxy77,Vxy78,
     &                                           Vxy79,Vxy710,Vxy711, 
     & Vxy81,Vxy82,Vxy83,Vxy84,Vxy85,Vxy86,Vxy87,Vxy88,
     &                                           Vxy89,Vxy810,Vxy811,  
     & Vxy91,Vxy92,Vxy93,Vxy94,Vxy95,Vxy96,Vxy97,Vxy98,
     &                                           Vxy99,Vxy910,Vxy911,
     & Vxy101,Vxy102,Vxy103,Vxy104,Vxy105,Vxy106,Vxy107,Vxy108,
     &                                        Vxy109,Vxy1010,Vxy1011,
     & Vxy11n1,Vxy11n2,Vxy11n3,Vxy11n4,Vxy11n5,Vxy11n6,Vxy11n7,Vxy11n8,
     &                                Vxy11n9,Vxy11n10,Vxy11n11

       Double precision 
     & Vx11,Vx12,Vx13,Vx14,Vx15,Vx16,Vx17,Vx18,Vx19,
     & Vx110,Vx111,Vx112,
     & Vx21,Vx22,Vx23,Vx24,Vx25,Vx26,Vx27,Vx28,Vx29,
     & Vx210,Vx211,Vx212,
     & Vx31,Vx32,Vx33,Vx34,Vx35,Vx36,Vx37,Vx38,Vx39,
     & Vx310,Vx311,Vx312,
     & Vx41,Vx42,Vx43,Vx44,Vx45,Vx46,Vx47,Vx48,Vx49,
     & Vx410,Vx411,Vx412,
     & Vx51,Vx52,Vx53,Vx54,Vx55,Vx56,Vx57,Vx58,Vx59,
     & Vx510,Vx511,Vx512,
     & Vx61,Vx62,Vx63,Vx64,Vx65,Vx66,Vx67,Vx68,Vx69,
     & Vx610,Vx611,Vx612,
     & Vx71,Vx72,Vx73,Vx74,Vx75,Vx76,Vx77,Vx78,Vx79,
     & Vx710,Vx711,Vx712,
     & Vx81,Vx82,Vx83,Vx84,Vx85,Vx86,Vx87,Vx88,Vx89,
     & Vx810,Vx811,Vx812,
     & Vx91,Vx92,Vx93,Vx94,Vx95,Vx96,Vx97,Vx98,Vx99,
     & Vx910,Vx911,Vx912,
     & Vx101,Vx102,Vx103,Vx104,Vx105,Vx106,Vx107,Vx108,Vx109,
     & Vx1010,Vx1011,Vx1012,
     & Vx11n1,Vx11n2,Vx11n3,Vx11n4,Vx11n5,Vx11n6,Vx11n7,Vx11n8,Vx11n9,
     & Vx11n10,Vx11n11,Vx11n12,
     & Vx12n1,Vx12n2,Vx12n3,Vx12n4,Vx12n5,Vx12n6,Vx12n7,Vx12n8,Vx12n9,
     & Vx12n10,Vx12n11,Vx12n12
c
      Double precision 
     & Vy11,Vy12,Vy13,Vy14,Vy15,Vy16,Vy17,Vy18,Vy19,Vy110,Vy111,Vy112,
     & Vy21,Vy22,Vy23,Vy24,Vy25,Vy26,Vy27,Vy28,Vy29,Vy210,Vy211,Vy212,
     & Vy31,Vy32,Vy33,Vy34,Vy35,Vy36,Vy37,Vy38,Vy39,Vy310,Vy311,Vy312,
     & Vy41,Vy42,Vy43,Vy44,Vy45,Vy46,Vy47,Vy48,Vy49,Vy410,Vy411,Vy412,
     & Vy51,Vy52,Vy53,Vy54,Vy55,Vy56,Vy57,Vy58,Vy59,Vy510,Vy511,Vy512,
     & Vy61,Vy62,Vy63,Vy64,Vy65,Vy66,Vy67,Vy68,Vy69,Vy610,Vy611,Vy612,
     & Vy71,Vy72,Vy73,Vy74,Vy75,Vy76,Vy77,Vy78,Vy79,Vy710,Vy711,Vy712,
     & Vy81,Vy82,Vy83,Vy84,Vy85,Vy86,Vy87,Vy88,Vy89,Vy810,Vy811,Vy812,
     & Vy91,Vy92,Vy93,Vy94,Vy95,Vy96,Vy97,Vy98,Vy99,Vy910,Vy911,Vy912,
     & Vy101,Vy102,Vy103,Vy104,Vy105,Vy106,Vy107,
     & Vy108,Vy109,Vy1010,Vy1011,Vy1012,
     & Vy11n1,Vy11n2,Vy11n3,Vy11n4,Vy11n5,Vy11n6,Vy11n7,
     & Vy11n8,Vy11n9,Vy11n10,Vy11n11,Vy11n12,
     & Vy12n1,Vy12n2,Vy12n3,Vy12n4,Vy12n5,Vy12n6,Vy12n7,
     & Vy12n8,Vy12n9,Vy12n10,Vy12n11,Vy12n12
c
       Double precision 
     & V11,V12,V13,V14,V15,V16,V17,V18,V19,V110,V111,V112,
     & V21,V22,V23,V24,V25,V26,V27,V28,V29,V210,V211,V212,
     & V31,V32,V33,V34,V35,V36,V37,V38,V39,V310,V311,V312,
     & V41,V42,V43,V44,V45,V46,V47,V48,V49,V410,V411,V412,
     & V51,V52,V53,V54,V55,V56,V57,V58,V59,V510,V511,V512,
     & V61,V62,V63,V64,V65,V66,V67,V68,V69,V610,V611,V612,
     & V71,V72,V73,V74,V75,V76,V77,V78,V79,V710,V711,V712,
     & V81,V82,V83,V84,V85,V86,V87,V88,V89,V810,V811,V812,
     & V91,V92,V93,V94,V95,V96,V97,V98,V99,V910,V911,V912,
     & V101,V102,V103,V104,V105,V106,V107,V108,V109,V1010,V1011,V1012, 
     & V11n1,V11n2,V11n3,V11n4,V11n5,V11n6,V11n7,V11n8,
     & V11n9,V11n10,V11n11,V11n12,
     & V12n1,V12n2,V12n3,V12n4,V12n5,V12n6,V12n7,V12n8,
     & V12n9,V12n10,V12n11,V12n12
   
        Double precision H,derH,der3H,g 

      double precision ff, r,z, sn,cs,sn2,cs2,csc ,ct 

	double precision nr,w,con,alpha,c1,c2,c3,lambda,rh,par
	
	common  nr,w,con,alpha,c1,c2,c3,lambda,rh,par

c_______________________________________________________________
c_______________________________________________________________
c_______________________________________________________________


      do 50 I=1,nv  
      pu(I)=0.  
      pux(I)=0.  
      puy(I)=0.  
      puxx(I)=0.  
      puxy(I)=0. 
      puyy(I)=0.  
   50 continue  

 
      do 9999 I=1,nv   
	
	U1=U(I,1)
      U2=U(I,2)
      U3=U(I,3)
      U4=U(I,4)
      U5=U(I,5) 

      U1x=UX(I,1)
      U2x=UX(I,2) 
      U3x=UX(I,3) 
      U4x=UX(I,4) 
      U5x=UX(I,5) 

      U1y=UY(I,1)
      U2y=UY(I,2) 
      U3y=UY(I,3) 
      U4y=UY(I,4) 
      U5y=UY(I,5) 

      U1yy=UYY(I,1)
      U2yy=UYY(I,2)
      U3yy=UYY(I,3)
      U4yy=UYY(I,4)
      U5yy=UYY(I,5) 


      U1xy=UXY(I,1)
      U2xy=UXY(I,2)
      U3xy=UXY(I,3)
      U4xy=UXY(I,4)
      U5xy=UXY(I,5) 
							
	U1xx=UXX(I,1)
      U2xx=UXX(I,2) 
      U3xx=UXX(I,3) 
      U4xx=UXX(I,4) 
      U5xx=UXX(I,5) 

       r= XX(I)/(1.-XX(I)) 

         sn=dSin(Y(I))
         cs=dCos(Y(I))
	csc=1.d0/dSin(Y(I))
         ct=dCos(Y(I))/dSin(Y(I))
	 sn2=dSin(2.*Y(I))
 	cs2=dCos(2.*Y(I))
 

	ff=1.-XX(I) 
       g=r**2 + rh**2 

       H= Sqrt(g)/(Sqrt(g) + rh)

	 derH= (r*rh)/(Sqrt(g)*(Sqrt(g) + rh)**2)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
  
       Vxx11=ff**4*r*sn**2
       Vxy11=0.
       Vyy11=(r*sn**2)/(g*H)
       Vx11=(0.5*ff**2*r*(derH*g + 2.*H*(-2.*ff*g + r))*sn**2)/(g*H)
       Vy11=0.
       V11=
     & (4.*dExp(2.*U1)*alpha**2*U4**2*
     &     ((-1.*H*nr**2*r**2)/dExp(2.*U2) + 
     &       (g**2*sn**2*(nr*U5 + w)**2)/dExp(2.*U3)))/
     &   (g*H**2*r)
     &
       Vxx12=0.
       Vxy12=0.
       Vyy12=0.
       Vx12=
     & (-0.5*ff**2*sn**2*(-2.*H*r**2 + 
     &       g*(derH*r + 2.*H*(1. + ff**2*r*U3x))))/(g*H)
     &
       Vy12=(-1.*r*sn**2*U3y)/(g*H)
       V12=
     & (4.*dExp(2.*U1 - 2.*U2)*alpha**2*nr**2*r*
     &      U4**2)/(g*H) - (0.5*dExp(2.*U2 - 2.*U3)*g*
     &      sn**4*(ff**4*g*H*U5x**2 + U5y**2))/(H**2*r)
     &
       Vxx13=0.
       Vxy13=0.
       Vyy13=0.
       Vx13=(-1.*ff**2*r*sn**2*(r + ff**2*g*U2x))/g
       Vy13=(-1.*r*sn*(cs + sn*U2y))/(g*H)
       V13=
     & (0.5*g*sn**2*(dExp(2.*U2)*sn**2*
     &        (ff**4*g*H*U5x**2 + U5y**2) - 
     &       8.*dExp(2.*U1)*alpha**2*U4**2*
     &        (nr*U5 + w)**2))/(dExp(2.*U3)*H**2*r)
     &
       Vxx14=0.
       Vxy14=0.
       Vyy14=0.
       Vx14=4.*alpha**2*ff**4*r*sn**2*U4x
       Vy14=(4.*alpha**2*r*sn**2*U4y)/(g*H)
       V14=
     & (4.*dExp(2.*U1)*alpha**2*U4*
     &     ((-1.*H*nr**2*r**2)/dExp(2.*U2) + 
     &       (g**2*sn**2*(nr*U5 + w)**2)/dExp(2.*U3)))/
     &   (g*H**2*r)
     &
       Vxx15=0.
       Vxy15=0.
       Vyy15=0.
       Vx15=
     & (-0.5*dExp(2.*U2 - 2.*U3)*ff**4*g**2*sn**4*U5x)/
     &   (H*r)
     &
       Vy15=
     & (-0.5*dExp(2.*U2 - 2.*U3)*g*sn**4*U5y)/(H**2*r)
     &
       V15=
     & (4.*dExp(2.*U1 - 2.*U3)*alpha**2*g*nr*sn**2*
     &     U4**2*(nr*U5 + w))/(H**2*r)
     &
       Vxx21=0.
       Vxy21=0.
       Vyy21=0.
       Vx21=0.
       Vy21=0.
       V21=
     & (4.*dExp(2.*U1)*alpha**2*U4**2*
     &     ((2.*nr**2)/(dExp(2.*U2)*g) + 
     &       sn**2*(c3 + c2*U4**2 + c1*U4**4)))/H
     &
       Vxx22=ff**4*sn**2
       Vxy22=0.
       Vyy22=sn**2/(g*H)
       Vx22=
     & (ff**2*sn**2*(2.*H*r**2 + 
     &       g*(derH*r + H*(1. - 2.*ff*r + ff**2*r*(2.*U2x + U3x)))))/
     &   (g*H*r)
     &
       Vy22=(sn*(2.*cs + sn*(2.*U2y + U3y)))/(g*H)
       V22=
     & (-8.*dExp(2.*U1 - 2.*U2)*alpha**2*nr**2*U4**2)/
     &    (g*H) + (dExp(2.*U2 - 2.*U3)*g*sn**4*
     &      (ff**4*g*H*U5x**2 + U5y**2))/(H**2*r**2)
     &
       Vxx23=0.
       Vxy23=0.
       Vyy23=0.
       Vx23=(ff**2*sn**2*(r + ff**2*g*U2x))/g
       Vy23=(sn*(cs + sn*U2y))/(g*H)
       V23=
     & (-1.*dExp(2.*U2 - 2.*U3)*g*sn**4*
     &     (ff**4*g*H*U5x**2 + U5y**2))/(H**2*r**2)
     &
       Vxx24=0.
       Vxy24=0.
       Vyy24=0.
       Vx24=0.
       Vy24=0.
       V24=
     & (4.*dExp(2.*U1 - 2.*U2)*alpha**2*U4*
     &     (2.*nr**2 + dExp(2.*U2)*g*sn**2*
     &        (c3 + 2.*c2*U4**2 + 3.*c1*U4**4)))/(g*H)
     &
       Vxx25=0.
       Vxy25=0.
       Vyy25=0.
       Vx25=
     & (dExp(2.*U2 - 2.*U3)*ff**4*g**2*sn**4*U5x)/
     &   (H*r**2)
     &
       Vy25=
     & (dExp(2.*U2 - 2.*U3)*g*sn**4*U5y)/(H**2*r**2)
     &
       V25=0.
       Vxx31=0.
       Vxy31=0.
       Vyy31=0.
       Vx31=0.
       Vy31=0.
       V31=
     & (8.*dExp(2.*U1)*alpha**2*sn*U4**2*
     &     (r**2*(c3 + c2*U4**2 + c1*U4**4) - 
     &       (2.*g*(nr*U5 + w)**2)/(dExp(2.*U3)*H)))/
     &   (H*r)
     &
       Vxx32=0.
       Vxy32=0.
       Vyy32=0.
       Vx32=
     & (ff**2*sn*(-2.*H*r**2 + g*(derH*r + 2.*H*(1. + ff**2*r*U3x))))/
     &   (g*H)
     &
       Vy32=(2.*r*sn*U3y)/(g*H)
       V32=
     & (-2.*dExp(2.*U2 - 2.*U3)*g*sn**3*
     &     (ff**4*g*H*U5x**2 + U5y**2))/(H**2*r)
     &
       Vxx33=2.*ff**4*r*sn
       Vxy33=0.
       Vyy33=(2.*r*sn)/(g*H)
       Vx33=
     & (ff**2*sn*(3.*derH*r + 2.*H*
     &        (2. - 2.*ff*r + ff**2*r*(U2x + 2.*U3x))))/H
     &
       Vy33=(2.*r*(cs + sn*(U2y + 2.*U3y)))/(g*H)
       V33=
     & (2.*g*sn*(dExp(2.*U2)*sn**2*
     &        (ff**4*g*H*U5x**2 + U5y**2) + 
     &       8.*dExp(2.*U1)*alpha**2*U4**2*
     &        (nr*U5 + w)**2))/(dExp(2.*U3)*H**2*r)
     &
       Vxx34=0.
       Vxy34=0.
       Vyy34=0.
       Vx34=0.
       Vy34=0.
       V34=
     & (8.*dExp(2.*U1 - 2.*U3)*alpha**2*sn*U4*
     &     (dExp(2.*U3)*c3*H*r**2 + 
     &       2.*dExp(2.*U3)*c2*H*r**2*U4**2 + 
     &       3.*dExp(2.*U3)*c1*H*r**2*U4**4 - 
     &       2.*g*nr**2*U5**2 - 4.*g*nr*U5*w - 2.*g*w**2))/(H**2*r)
     &
       Vxx35=0.
       Vxy35=0.
       Vyy35=0.
       Vx35=
     & (-2.*dExp(2.*U2 - 2.*U3)*ff**4*g**2*sn**3*U5x)/
     &   (H*r)
     &
       Vy35=
     & (-2.*dExp(2.*U2 - 2.*U3)*g*sn**3*U5y)/(H**2*r)
     &
       V35=
     & (-16.*dExp(2.*U1 - 2.*U3)*alpha**2*g*nr*sn*
     &     U4**2*(nr*U5 + w))/(H**2*r)
     &
       Vxx41=0.
       Vxy41=0.
       Vyy41=0.
       Vx41=0.
       Vy41=0.
       V41=
     & (2.*dExp(2.*U1)*U4*
     &     ((-1.*H*nr**2)/dExp(2.*U2) - 
     &       1.*g*H*sn**2*(c3 + 2.*c2*U4**2 + 3.*c1*U4**4) + 
     &       (g**2*sn**2*(nr*U5 + w)**2)/
     &        (dExp(2.*U3)*r**2)))/H**2
     &
       Vxx42=0.
       Vxy42=0.
       Vyy42=0.
       Vx42=ff**4*g*sn**2*U4x
       Vy42=(sn**2*U4y)/H
       V42=(2.*dExp(2.*U1 - 2.*U2)*nr**2*U4)/H
       Vxx43=0.
       Vxy43=0.
       Vyy43=0.
       Vx43=ff**4*g*sn**2*U4x
       Vy43=(sn**2*U4y)/H
       V43=
     & (-2.*dExp(2.*U1 - 2.*U3)*g**2*sn**2*U4*
     &     (nr*U5 + w)**2)/(H**2*r**2)
     &
       Vxx44=ff**4*g*sn**2
       Vxy44=0.
       Vyy44=sn**2/H
       Vx44=
     & ff**2*sn**2*(-2.*ff*g + r + 
     &     g*(derH/H + 1/r + ff**2*(U2x + U3x)))
     &
       Vy44=(sn*(cs + sn*(U2y + U3y)))/H
       V44=
     & (dExp(2.*U1 - 2.*(U2 + U3))*
     &     (-1.*dExp(2.*U3)*H*nr**2*r**2 - 
     &       1.*dExp(2.*(U2 + U3))*g*H*r**2*sn**2*
     &        (c3 + 6.*c2*U4**2 + 15.*c1*U4**4) + 
     &       dExp(2.*U2)*g**2*sn**2*(nr*U5 + w)**2))/
     &   (H**2*r**2)
     &
       Vxx45=0.
       Vxy45=0.
       Vyy45=0.
       Vx45=0.
       Vy45=0.
       V45=
     & (2.*dExp(2.*U1 - 2.*U3)*g**2*nr*sn**2*U4*
     &     (nr*U5 + w))/(H**2*r**2)
     &
       Vxx51=0.
       Vxy51=0.
       Vyy51=0.
       Vx51=0.
       Vy51=0.
       V51=
     & (-16.*dExp(2.*U1 - 2.*U2)*alpha**2*nr*r**2*
     &     U4**2*(nr*U5 + w))/(g*H)
     &
       Vxx52=0.
       Vxy52=0.
       Vyy52=0.
       Vx52=3.*ff**4*r**2*sn**2*U5x
       Vy52=(3.*r**2*sn**2*U5y)/(g*H)
       V52=
     & (16.*dExp(2.*U1 - 2.*U2)*alpha**2*nr*r**2*
     &     U4**2*(nr*U5 + w))/(g*H)
     &
       Vxx53=0.
       Vxy53=0.
       Vyy53=0.
       Vx53=-1.*ff**4*r**2*sn**2*U5x
       Vy53=(-1.*r**2*sn**2*U5y)/(g*H)
       V53=0.
       Vxx54=0.
       Vxy54=0.
       Vyy54=0.
       Vx54=0.
       Vy54=0.
       V54=
     & (-16.*dExp(2.*U1 - 2.*U2)*alpha**2*nr*r**2*U4*
     &     (nr*U5 + w))/(g*H)
     &
       Vxx55=ff**4*r**2*sn**2
       Vxy55=0.
       Vyy55=(r**2*sn**2)/(g*H)
       Vx55=
     & (ff**2*r*sn**2*(5.*r**2 + 
     &       g*(-1. - 2.*ff*r + ff**2*r*(3.*U2x - 1.*U3x))))/g
     &
       Vy55=(r**2*sn*(3.*cs + 3.*sn*U2y - 1.*sn*U3y))/(g*H)
       V55=
     & (-8.*dExp(2.*U1 - 2.*U2)*alpha**2*nr**2*r**2*
     &     U4**2)/(g*H)
    

    

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                  
c        write(6,*)r,V22
      GO TO (100,200,300,400,500  ) IEQU 
      

  100 continue
c          write(6,*)'ec.1'
      GO TO (110,120,130,140,150  ) ICOM 

  110 continue 
      PUXX(I)=Vxx11
      PUYY(I)=Vyy11
      PUX(I)=Vx11
      PUY(I)=Vy11
      PU(I)=V11 
      GO TO 9999
      
C           derivate in raport cu H2
  120  continue

     
      PUXX(I)=Vxx12
      PUYY(I)=Vyy12
      PUX(I)=Vx12
      PUY(I)=Vy12
      PU(I)=V12 
	GO TO 9999

  130  continue
      PUXX(I)=Vxx13
      PUYY(I)=Vyy13
      PUX(I)=Vx13
      PUY(I)=Vy13
      PU(I)=V13 
      GO TO 9999

  140  continue    
      PUXX(I)=Vxx14
      PUYY(I)=Vyy14
      PUX(I)=Vx14
      PUY(I)=Vy14
      PU(I)=V14 
      GO TO 9999

  150  continue     
      PUXX(I)=Vxx15
      PUYY(I)=Vyy15
      PUX(I)=Vx15
      PUY(I)=Vy15
      PU(I)=V15 
	GO TO 9999 
 

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C            a doua ecuatie 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  200 continue
c          write(6,*)'ec.1'
      GO TO (210,220,230,240,250 ) ICOM 

  210 continue 
      PUXX(I)=Vxx21
      PUYY(I)=Vyy21
      PUX(I)=Vx21
      PUY(I)=Vy21
      PU(I)=V21 

c            write(6,*)Vxx11,Vyy11,Vx11
      GO TO 9999
      
C           derivate in raport cu H2
  220  continue

     
      PUXX(I)=Vxx22
      PUYY(I)=Vyy22
      PUX(I)=Vx22
      PUY(I)=Vy22
      PU(I)=V22 
      
	GO TO 9999

  230  continue

     
      PUXX(I)=Vxx23
      PUYY(I)=Vyy23
      PUX(I)=Vx23
      PUY(I)=Vy23
      PU(I)=V23 

      GO TO 9999

  240  continue

     
      PUXX(I)=Vxx24
      PUYY(I)=Vyy24
      PUX(I)=Vx24
      PUY(I)=Vy24
      PU(I)=V24 

      GO TO 9999

  250  continue     
      PUXX(I)=Vxx25
      PUYY(I)=Vyy25
      PUX(I)=Vx25
      PUY(I)=Vy25
      PU(I)=V25 

	GO TO 9999 

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C            a TREIA ecuatie 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  300 continue
 
      GO TO (310,320,330,340,350  ) ICOM

c   prima ecuatie si prima componenta: U(I,1)=H1

  310 continue 
      PUXX(I)=Vxx31
      PUYY(I)=Vyy31
      PUX(I)=Vx31
      PUY(I)=Vy31
      PU(I)=V31 
c            write(6,*)Vxx11,Vyy11,Vx11
      GO TO 9999
      
C           derivate in raport cu H2
  320  continue

     
      PUXX(I)=Vxx32
      PUYY(I)=Vyy32
      PUX(I)=Vx32
      PUY(I)=Vy32
      PU(I)=V32 
      
	GO TO 9999

  330  continue 
      PUXX(I)=Vxx33
      PUYY(I)=Vyy33
      PUX(I)=Vx33
      PUY(I)=Vy33
      PU(I)=V33 

      GO TO 9999

  340  continue
      PUXX(I)=Vxx34
      PUYY(I)=Vyy34
      PUX(I)=Vx34
      PUY(I)=Vy34
      PU(I)=V34 
      GO TO 9999

  350  continue
      PUXX(I)=Vxx35
      PUYY(I)=Vyy35
      PUX(I)=Vx35
      PUY(I)=Vy35
      PU(I)=V35 
	GO TO 9999 
   

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C            a PATRA ecuatie 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  400 continue
 
      GO TO (410,420,430,440,450 ) ICOM

c   prima ecuatie si prima componenta: U(I,1)=H1

  410 continue 
      PUXX(I)=Vxx41
      PUYY(I)=Vyy41
      PUX(I)=Vx41
      PUY(I)=Vy41
      PU(I)=V41 
c            write(6,*)Vxx11,Vyy11,Vx11
      GO TO 9999
      
C           derivate in raport cu H2
 420  continue
      PUXX(I)=Vxx42
      PUYY(I)=Vyy42
      PUX(I)=Vx42
      PUY(I)=Vy42
      PU(I)=V42 
      
	GO TO 9999

 430  continue
      PUXX(I)=Vxx43
      PUYY(I)=Vyy43
      PUX(I)=Vx43
      PUY(I)=Vy43
      PU(I)=V43 

      GO TO 9999

  440  continue
      PUXX(I)=Vxx44
      PUYY(I)=Vyy44
      PUX(I)=Vx44
      PUY(I)=Vy44
      PU(I)=V44 
      GO TO 9999

  450  continue  
      PUXX(I)=Vxx45
      PUYY(I)=Vyy45
      PUX(I)=Vx45
      PUY(I)=Vy45
      PU(I)=V45 
	GO TO 9999 

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C            a CINCEA ecuatie 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  500 continue
 
      GO TO (510,520,530,540,550 ) ICOM

c   prima ecuatie si prima componenta: U(I,1)=H1

  510 continue 
      PUXX(I)=Vxx51
      PUYY(I)=Vyy51
      PUX(I)=Vx51
      PUY(I)=Vy51
      PU(I)=V51 

c            write(6,*)Vxx11,Vyy11,Vx11
      GO TO 9999
      
C           derivate in raport cu H2
 520  continue 
      PUXX(I)=Vxx52
      PUYY(I)=Vyy52
      PUX(I)=Vx52
      PUY(I)=Vy52
      PU(I)=V52 
      
	GO TO 9999

 530  continue 
      PUXX(I)=Vxx53
      PUYY(I)=Vyy53
      PUX(I)=Vx53
      PUY(I)=Vy53
      PU(I)=V53 

      GO TO 9999

  540  continue    
      PUXX(I)=Vxx54
      PUYY(I)=Vyy54
      PUX(I)=Vx54
      PUY(I)=Vy54
      PU(I)=V54 
      GO TO 9999

  550  continue 
      PUXX(I)=Vxx55
      PUYY(I)=Vyy55
      PUX(I)=Vx55
      PUY(I)=Vy55
      PU(I)=V55 
	GO TO 9999  
9999   continue
  


 1111 continue 
      return 
      end 

    
 



