c             25.02.13
 

      PARAMETER (NX= 251 , NY= 30  , NK= 5 , N=NX*NY ,LOUT=6, NXMAX=200,
     &           NOX= 4 , NOY= 4 , MNOX= 8 , MNOY= 8,
     &           M=NXMAX*NXMAX*NK , MV=5*N,
     &           M10=6*NK+(2+6*NK)*MV , M20=1+NK+(7+5*NK)*MV,
     &           MD=2*M-1 , MLI=10+NK,                                          
     &           MIW=3*MV, MDI=2*NXMAX+5+2*NK,MDW=1885163  )
C                                                                               
      INTEGER     RL4                                                           
C                                                                               
C****            METHOD ....   METHOD FOR THE LINEAR EQUATION SOLVER            
C****            IPACK  ....   SEE ILIN(14) IN YOUR DOCUMENTATION               
C****            PACK   ....   SPARSITY OF THE MATRIX                           
C****            IDP    ....   NUMBER OF DIAGONALPARTS
C****            ID     ....   NUMBER OF DIAGONALS                              
C                                                                               
      PARAMETER (METHOD=2  , IPACK= 2 , PACK = 0.4,
C    &           IDP  = NK*NK*(2.d0*(MNOX+MNOY)+5),                                
     &           IDP  = 173,                                                    
C    &           ID   = (2.d0*MNOX+2)*NK+(2.d0*MNOY+3)*(2.d0*NK-1),                      
     &           ID   = 83,
C    &           NDQP = ID * M * 0.67 ,
     &           NDQP = 378600,                                                  
C                                                                               
C****            FILE SIZES
C
C****            IZ     ....   NUMBER OF PERIODIC DIRECTIONS
     &           RL4=NXMAX , IZ=0.d0 , IN=(N+MV-1)/MV,                             
     &           NF1 = 40 + 2*RL4 + 2*NK + 2*M+2,                                 
     &           NF2 =2*( ((MNOX+MNOY)*4+16)*IN-IZ*6*IN+14 ) * (MV+10),
     &           NF3 = 7*IN*NK * MV+2,
     &           NF4 = 30*NK*(1+NK) * RL4+2 ,
     &           NF5 = (4+IN)*NK * MV+2 ,
C**** SHOULD BE EQUAL MAX0(4.d0IN,4+IN)*NK * MV
     &           NF6 = 80*IDP   * N+2 ,
C**** SHOULD BE EQUAL MAX0(2+5*NK,IDP)*IN * MV
     &           NF7= 5*NK*NK*IN * MV+2,                                          
     &           NF89=  ID*M*(2+IPACK*(PACK-1.0))+2 ,                             
     &           NF10=M+2 ,
     &           MBIG=NF1+NF2+NF3+NF4+NF5+NF6+NF7+ NF89 + NF10)                 
C
      INTEGER  I,I1,I2,IFILE(10),IND,IER,IW(MIW),IINFO(35),                     
     &         ILIN(17),GL,ii,maxiti,methoda
C    

c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
c           for FDS01                                                                           
        DOUBLE PRECISION DINFO(MDI),DW(MDW),U1(NX,NY,NK),
     &         U2(NX,NY,NK),BYTES(10) 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU



 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
c           for FDS02    
C        DOUBLE PRECISION DINFO(MDI),DW(MDW),U1(NXMAX,NXMAX,NK),
C     &         U2(NXMAX,NXMAX,NK),BYTES(10) 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU       
	 
	DOUBLE PRECISION   eta,nr,w,con,alpha,c1,c2,c3,lambda,
     & x,y,k1,k2,k3,k4,fi1,fi2,a1,r,s,rh,par,g
 

      LOGICAL  LINFO(MLI),LLIN(MD)                                              
      EXTERNAL  U01E1,U02E1,U03E1,U04E1,FDS500
C
      COMMON /FDSCON/  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
                                                                                   
C
      DOUBLE PRECISION BIG
      COMMON/FDSIO/ BIG(MBIG)				                                         
C                                                                               
      INTEGER NBIG,NTOTAL,IOINF1,IOINF2,IOINF3                                  
      COMMON/FDSINF/ NBIG,NTOTAL,IOINF1(99),IOINF2(99),IOINF3(99)
C
C***  THE COMMON-BLOCKS FDSLEN ,FDSIND AND FDSDQP ARE ONLY USED FOR             
C     PACKED STORING OF THE NORMALIZED MATRIX                                   
C                                                                               
      INTEGER  IDQP,NIND,NPACK
      DOUBLE PRECISION    DQP                                                              
      COMMON/FDSLEN/ NPACK(IDP)
      COMMON/FDSIND/ NIND,IDQP(NDQP)                                            
      COMMON/FDSDQP/ DQP(M)   
               
 
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
c  extremal BH: 0.81<w<0.95
 
 			 rh=0.1d0
 
          nr = 1.d0
				
	    w= 0.954d0 

c literature
c          c1=1.d0
c          c2=-2.d0
c          c3=1.1d0

c      * coefficient Z^6
          c1=1.d0

c      * coefficient Z^4
          c2= -2.d0

c      * coefficient Z^2
          c3=1.1d0
 
  
	alpha=1.d0
 
               methoda = 2
                maxiti =    42520  


            write(6,*) 'n= ', nr  
            write(6,*)'w= ', w 
            write(6,*)'c1= ', c1
            write(6,*)'c2= ', c2 
            write(6,*)'c3= ', c3 		  	
		  write(6,*)'alpha= ', alpha  
 		  write(6,*)' '          

		  write(6,*)' '   

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
          Open(Unit=1,File='gridx.dat')  
                                          
      DO 17 I=1,NX
       Read(1,789)s                                                                 
                           DINFO(I)=s/(s+1. ) 
c            write(6,*)'grid x',s,DINFO(I)                                       
  17  CONTINUE 

                  Close(Unit=1)     
           
		 DINFO(1)=0.d0
		 DINFO(nx)=1.d0
        
  	   Open(Unit=1,File='gridy.dat')  

      DO 18 I=1,Ny
       Read(1,789)s                                                                 
                            DINFO(nx+I)=s 
c			    write(6,*)'grid y',I,DINFO(nx+I)                                    
  18  CONTINUE
                   
                   Close(Unit=1)   
  
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C  !!!!!!!!!!!!!    BEGINNING OF ITERATION     !!!!!!!!!!!!!!!!!!!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

      Do 99 ii= 1, 2


c 	      c1= 0.6d0 + (ii-1)*0.1d0
c   	      c2= -1.6d0 - (ii-1)*0.1d0
c           nr= 2.05d0 +(ii-1)*0.05d0 
c           w= 0.84d0-(ii-1)*0.01d0 
c           rh= 0.11d0 +(ii-1)*0.01d0 


      IOINF1(21) = NF1
      IOINF1(22) = NF2
      IOINF1(23) = NF3
      IOINF1(24) = NF4
      IOINF1(25) = NF5
      IOINF1(26) = NF6
      IOINF1(27) = NF7
      IOINF1(28) = 0
c       write(6,*)Ipack
c      IF( IPACK .LE. 1 ) THEN
c        IOINF1(28) = ID*M
c      ELSE
c        IOINF1(28) = ID*M*PACK
c      ENDIF

      IF( IPACK .EQ. 0 ) THEN
        IOINF1(29) = ID*M
      ELSE
        IOINF1(29) = 2*ID*M*PACK
      ENDIF
      IOINF1(30) = NF10
      NIND       = NDQP
      NBIG       = MBIG
      NTOTAL     = 0
C
      IF ( NF5 .LT. (MAX0(4*IN,4+IN)*NK * MV) )
     &   WRITE(6,1010) NF5, (MAX0(4*IN,4+IN)*NK * MV)
      IF ( NF6 .LT.  (MAX0(2+5*NK,IDP)*IN * MV) )
     &   WRITE(6,1020) NF6, (MAX0(2+5*NK,IDP)*IN * MV)
C
C***  PARAMETERS FOR FILES
C
      IF     ( (METHOD .LE. 1) .OR. (METHOD .EQ. 10) ) THEN
         GL = 15
      ELSEIF ((METHOD .EQ. 2) .OR. (METHOD .EQ. 3) ) THEN
         GL = 11
      ELSEIF  (METHOD .EQ. 4) THEN
         GL = 10
      ELSE
         WRITE(6,999)
      ENDIF
C
      DO 111 J=1,10
         BYTES(J)  = DBLE(IOINF1(20+J))* 8. * 1.D-6
  111 CONTINUE
C
      WRITE(6,1030) N,M,NDQP,MV,MBIG,IDP,ID
      WRITE(6,1040)
      WRITE(6,1050) (J, IOINF1(20+J),BYTES(J) ,J=1,10 )
C
C***  UNIT - NUMBERS
C
      IFILE(1)=21
      IFILE(2)=22
      IFILE(3)=23
      IFILE(4)=24
      IFILE(5)=25
      IFILE(6)=26
      IFILE(7)=27
      IFILE(8)=28
      IFILE(9)=29
      IFILE(10)=30
C
C***INFORMATION-VECTORS :



      IINFO(1)  = NK
      IINFO(2)  = NX
      IINFO(3)  = NY
      IINFO(5)  = MV
      IINFO(6)  = 1000
      IINFO(7)  = LOUT
      IINFO(8)  =  0
      IINFO(9)  = MDW
      IINFO(10) = MIW
      IINFO(11) = NOX
      IINFO(12) = NOY
      IINFO(14) = MNOX
      IINFO(15) = MNOX
      IINFO(17) = NXMAX
c
      LINFO(1) = .TRUE.
      LINFO(2) = .FALSE.
      LINFO(3) = .FALSE.
c true - se ia step 0
c       LINFO(5) = .FALSE.
      LINFO(5) = .true.
      LINFO(6) = .FALSE.
      LINFO(7) = .TRUE.
      LINFO(8) = .FALSE.
      LINFO(9) = .false.
      LINFO(10) = .false.
      LINFO(11) = .FALSE.
      LINFO(12) = .FALSE.
      LINFO(13) = .FALSE.                                                       
      LINFO(14) = .FALSE.                                                       
      LINFO(15) = .FALSE.   
      LINFO(16) = .FALSE.            
 
      ILIN(13) = 3
      ILIN(14) = 0.
      ILIN(15) = 0.

c               PARAMETRI:
C
C
C     toleranta: desired relative accuracy for the
c                computation of the solution

       DINFO(NX+NY+1) =  0.1d0
c      ILIN(1) = METHOD
C          10 <=> Multialgoritm;
 
        ILIN(1)  = methoda

c  maximum iterations
         ILIN(2)  =    maxiti

C***  PARAMETERS OF THE PROBLEM


c       write(6,*)'NK=',NK
       write(6,*)'111111111111111111111111111111111111111111111'
       write(6,*)'  '
            write(6,*)'start iteration for i= ', ii
            write(6,*)' ' 
            write(6,*)'n= ', nr 
            write(6,*)'rh= ',rh
            write(6,*)'w= ', w 
            write(6,*)'c1= ', c1
            write(6,*)'c2= ', c2 
            write(6,*)'c3= ', c3 		  	
		  write(6,*)'alpha= ', alpha  
 		  write(6,*)' '   

       write(6,*)'  '
       write(6,*)'  '
	 write(6,*)"i=",ii
       write(6,*)'  '

      Open(Unit=1,File='functf.dat')
  
                                                                               
      DO 50 I2=1,NY                                                             
      DO 40 I1=1,NX  
c  U1 = f1
c  U2 = f2
c  U3 = f0
c  U4 = Z
c  U5 = W
        
	  r= DINFO(I1)/(1.-DINFO(I1)+1.d-12)+1.d-12
	g=r**2+rh**2
c	write(6,*) r
	      
                               Read(1,789)a1                                                      
         U1(I1,I2,1) = a1 
                               Read(1,789)a1                                                            
         U1(I1,I2,2) = a1 
 	                         Read(1,789)a1                                                            
         U1(I1,I2,3) = a1   
      		                 Read(1,789)a1  
	   U1(I1,I2,4) = a1    
       		                 Read(1,789)a1       
 	   U1(I1,I2,5) = a1   
     		                     
  
  40  CONTINUE

  50  CONTINUE

  
      Close(Unit=1)
C
C***CALL OF FIDISOL
C

c      CALL SECOND(TIME)
c      WRITE(6,1900) TIME
                    
      CALL FDSE01(IINFO,DINFO,LINFO,ILIN,LLIN,IW,DW,IFILE,
     &         U1,U2,IND,IER,U01E1,U02E1,U03E1,U04E1,FDS500)
 


        rewind(1)
        open(unit=1,File='functf.dat')
        do 91 I2=1,NY
           DO 90 I1=1,NX
        WRITE(1,*)  U1(I1,I2,1)
        WRITE(1,*)  U1(I1,I2,2) 
        WRITE(1,*)  U1(I1,I2,3) 
        WRITE(1,*)  U1(I1,I2,4) 
        WRITE(1,*)  U1(I1,I2,5)  	 	  	   
   90  CONTINUE
   91  continue
            


         rewind(1)
        open(unit=1,File='errf.dat')
        do 191 I2=1,NY
           DO 190 I1=1,NX
        WRITE(1,789)  U2(I1,I2,1)
        WRITE(1,789)  U2(I1,I2,2)  
        WRITE(1,789)  U2(I1,I2,3)  
        WRITE(1,789)  U2(I1,I2,4)  
        WRITE(1,789)  U2(I1,I2,5) 	   
  190  CONTINUE
  191  continue



       write(6,*)'  '
       write(6,*)'  '
            write(6,*)'end iteration for i= ', ii
            write(6,*)' ' 
            write(6,*)'n= ', nr 
            write(6,*)'rh= ',rh
            write(6,*)'w= ', w 
            write(6,*)'c1= ', c1
            write(6,*)'c2= ', c2 
            write(6,*)'c3= ', c3 		  	
		  write(6,*)'alpha= ', alpha  
 		  write(6,*)' '    
       write(6,*)'____________________________________________'

 
              Close(Unit=1)
          
 
        rewind(1)
        open(unit=1,File='funct.dat')
        DO 1170 I2=1,Ny
          DO 1171 I1=1,Nx
           

      write(1,765) DINFO(I1), DINFO(nx+I2), 
     &U1(I1,I2,1),U1(I1,I2,2),U1(I1,I2,3),U1(I1,I2,4),U1(I1,I2,5) 

 1171      CONTINUE 
          write(1,*) 
 1170    CONTINUE
 		 

         open(unit=1,File='err.dat')
        DO 270 I2=1,Ny
          DO 271 I1=1,Nx
      write(1,765) DINFO(I1), DINFO(nx+I2), 
     &U2(I1,I2,1),U2(I1,I2,2),U2(I1,I2,3),U2(I1,I2,4),U2(I1,I2,5)

 271      CONTINUE 
          write(1,*) 
 270    CONTINUE

             Close(Unit=1)

            open(unit=5,File='res.txt')
 	        WRITE(5,765)   nr,w,alpha,c1,c2,c3,rh
 	        Close(Unit=5)


   99  continue
              
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C  !!!!!!!!!!!!!    END ITERATION !!!!!!!!!!!!!!!!!!!
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       write(6,*)'  FINAL   '
       WRITE(6,*)' winding number n= ', nr


  999 FORMAT (//,' ****THIS PROGRAM DOES NOT EXIST YET FOR THE METHOD '
     *,' CHOSEN',/,' ****WE WILL BE GLAD IF YOU COULD DO IT')
 1010 FORMAT (1X,'WARNING!! SIZE OF FILE 5 ',I10,' SHOULD BE >= ',I10)
 1020 FORMAT (1X,'WARNING!! SIZE OF FILE 6 ',I10,' SHOULD BE >= ',I10)
 1030 FORMAT (1X,'N       : ',I12,/
     *   ,1X,'M       : ',I12,8X,'IDQP    : ',I12,' WORDS',/
     *   ,1X,'MV      : ',I12,8X,'MBIG    : ',I12,' WORDS',/
     *   ,1X,'IDP     : ',I12,8X,'ID      : ',I12,/)
 1040 FORMAT (19X,' WORDS (= 8 BYTES )',9X,'MBYTES')
 1050 FORMAT (10(1X,'SIZE FOR FILE : ',I2,3X,I12,4X,F12.2,/))
 1900 FORMAT(' TIME:',G15.6)
 1910 FORMAT(' TIME:',G15.6,'  IND,IER,IMVM:',3I10)
C
 2001 FORMAT(//'  SOLUTION-COMPONENT:',I4//'  GRID:  /',
     &    5(F8.4,6X)/1X,79('-'))
 2002 FORMAT(1X,F7.4,' /',5(E12.4,2X))
c 2003 FORMAT(//'  ERROR-ESTIMATES:')
C
  789     format (e24.16)
  766     format (8(e24.16))
c  765     format (7(e24.16))
  765  FORMAT (7(e24.16,1x))
           S T O P


C
       E N D
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
C    23.08.99 -INCEPUT
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C**********************************************************************
C**                                                             ***
C**                                                             ***
      SUBROUTINE  U01E1(T,X,Y,U,UT,UX,UXX,UY,UYY,P,MT,MV,NK,NV)
 
C**                                                             ***
C**                                                             ***
C**********************************************************************
C**                                                             ***
C**  F D S U 0 1   PROGRAM FRAME FOR THE USER-SUBROUTINE,       ***
C**  IN WHICH THE PDE-SYSTEM AT INNER GRIDPOINTS IS PRESCRIBED. ***
C**                                                             ***
C**********************************************************************
C**                                                             ***
C**  FORMAL PARAMETERS :                                        ***
C**  -------------------                                        ***
C**                                                             ***
      INTEGER  MT, MV, NK, NV
      DOUBLE PRECISION  T,X(MV),Y(MV),P(MV,NK),U(MV,NK),UT(MT,NK),
     &               UX(MV,NK),UXX(MV,NK),UY(MV,NK),UYY(MV,NK),
     &               En1(MV),En2(MV),En3(MV), ptest(nk)

       Double precision U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,U12
       Double precision U1x,U2x,U3x,U4x,U5x,U6x,U7x,U8x,
     & U9x,U10x,U11x,U12x
	 Double precision U1xx,U2xx,U3xx,U4xx,U5xx,U6xx,U7xx,
     & U8xx,U9xx,U10xx,U11xx,U12xx
       Double precision U1y,U2y,U3y,U4y,U5y,U6y,U7y,U8y,U9y,
     & U10y,U11y,U12y
	 Double precision U1yy,U2yy,U3yy,U4yy,U5yy,U6yy,U7yy,
     & U8yy,U9yy,U10yy,U11yy,U12yy
       Double precision sn,cs,csc,ct,r,T34,T35a,T53a,ff  
       Double precision eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,e11,e12
	 Double precision r1,r2,r3,r4,r5,r6,r7,r8,r9,pat1,Z2
	 Double precision err1,err2,err3,err4,err5,err6,err7,err8,err9
	Double precision norm,rn,xn,max,rm,min,rmin,minf,rminf

	Double precision H,derH,der2H,g 
c
C**                                                             ***
C**-----------------------------------------------------------------***
C**                                                             ***
C** LIST OF FORMAL PARAMETERS : (SEE DOCUMENTATION 1.4 A)         ***
C** ---------------------------                                   ***
C**                                                             ***
C**-----------------------------------------------------------------***
C**                                                             ***
C**    LOCAL PARAMETERS :                                       ***
C**    ----------------                                         ***
C**                                                             ***
      INTEGER  I
      DOUBLE PRECISION  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
      COMMON /FDSCON/  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
C**                                                             ***
C**-----------------------------------------------------------------***
C
C
C****START OF CALCULATION :
C     ---------------------
C 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     lista variabilelor:
c     x=xbara=x(r)/(1+x(r))
c     Y(I)=teta
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C              MATTER EQUATIONS
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

      err1=0.d0 
      err2=0.d0  
	norm=1.
	max=-1. 
	min=2.
	minf=2.
c	write(6,*)eta,nr,w,con,g,lambda1,lambda2,eta1,const
c                       write(6,*)"TH*IS HERE"



         rewind(1)
        open(unit=1,name='T44.dat')


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
 
	 
         r= X(I)/(1.-X(I)) 

         sn=dSin(Y(I))
         cs=dCos(Y(I))
	csc=1.d0/dSin(Y(I))
         ct=dCos(Y(I))/dSin(Y(I))
	 sn2=dSin(2.*Y(I))
 	cs2=dCos(2.*Y(I))
 

	ff=1.-X(I)
        
   
       g=r**2 + rh**2 

       H= rh + Sqrt(g)

	 derH= r/Sqrt(g)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
 
      eq1=(0.25*(-8.*dExp(2.*(U1 + U3))*alpha**2*g*
     -       H**3*nr**2*r**4*U4**2 + 
     -      4.*dExp(2.*(U2 + U3))*g*r**4*sn*
     -       (-2.*ff**3*g*H*r**2*sn*U1x + 
     -         ff**2*r*sn*(g*
     -             (H*(U1x - 2.*U2x) + derH*r*(-1.*U1x + U2x)) + 
     -            H*r**2*(U1x + U2x - 1.*U3x)) + 
     -         ff**4*g*H*r**2*sn*
     -          (U1xx - 1.*U2x*U3x + 2.*alpha**2*U4x**2) + 
     -         H**3*(-1.*cs*U3y + 
     -            sn*(U1yy - 1.*U2y*U3y + 2.*alpha**2*U4y**2)))
     -       - 1.*dExp(4.*U2)*H**3*sn**4*
     -       (r**2*(-2.*r*U5 + ff**2*g*U5x)**2 + g*H**2*U5y**2)
     -       + 8.*dExp(2.*(U1 + U2))*alpha**2*g*
     -       H**5*sn**2*U4**2*(nr*U5 + g*w)**2))/
     -  (dExp(2.*(U2 + U3))*g**2*H*r**4) 

      eq2=  (0.5*(2.*sn*(H*r**2*rh**2*sn + 
     -         ff**2*g**2*r*sn*
     -          (-2.*derH*r*U2x + 
     -            H*(ff**2*r*U2x**2 + ff**2*r*U2xx + 
     -               U2x*(3. + ff*r*(-2. + ff*U3x)))) + 
     -         g*(-2.*derH*r**3*sn + 
     -            H*r**2*sn*(3. + ff**2*r*(2.*U2x + U3x)) + 
     -            H**3*(cs*(2.*U2y + U3y) + 
     -               sn*(-1. + U2yy + U2y*(U2y + U3y))))) + 
     -      4.*dExp(2.*U1)*alpha**2*g*H**3*U4**2*
     -       ((2.*nr**2)/dExp(2.*U2) + 
     -         g*sn**2*(c3 + c2*U4**2 + c1*U4**4)) + 
     -      (dExp(2.*U2 - 2.*U3)*H**3*sn**4*
     -         (r**2*(-2.*r*U5 + ff**2*g*U5x)**2 + g*H**2*U5y**2)
     -         )/r**4))/(g**2*H)  

      eq3=  (0.5*sn*(2.*dExp(2.*U3)*r**4*
     -       (-1.*H**2*r**2*(r**2 + 2.*rh**2)*sn + 
     -         g**2*sn*(3.*derH**2*r**2 - 
     -            1.*derH*ff**2*H*r**2*(U2x + 3.*U3x) + 
     -            H**2*(4. + 
     -               ff**2*r*
     -                (U2x*(2. + ff**2*r*U3x) + 
     -                  U3x*(5. + ff*r*(-2. + ff*U3x)) + 
     -                  ff**2*r*U3xx))) + 
     -         g*H*(-2.*derH*r*(3.*r**2 + 4.*rh**2)*sn - 
     -            1.*ff**2*H*r**3*sn*U2x + 
     -            H**3*(U3y*(cs + sn*(U2y + U3y)) + sn*U3yy))) + 
     -      4.*dExp(2.*(U1 + U3))*alpha**2*g**2*
     -       H**4*r**4*sn*U4**2*(c3 + c2*U4**2 + c1*U4**4) - 
     -      1.*dExp(2.*U2)*H**4*sn**3*
     -       (r**2*(-2.*r*U5 + ff**2*g*U5x)**2 + g*H**2*U5y**2)
     -       - 8.*dExp(2.*U1)*alpha**2*g*H**6*sn*
     -       U4**2*(nr*U5 + g*w)**2))/
     -  (dExp(2.*U3)*g**2*H**2*r**4)
      
	eq4=  (-1.*dExp(2.*U1 - 2.*U2)*H**2*nr**2*U4)/
     -   g - 1.*dExp(2.*U1)*H**2*sn**2*U4*
     -   (c3 + 2.*c2*U4**2 + 3.*c1*U4**4) + 
     -  (sn*(-2.*ff**3*g*H*r**2*sn*U4x + 
     -       ff**2*r*(3.*g*H - 2.*derH*g*r + H*r**2)*sn*U4x + 
     -       ff**4*g*H*r**2*sn*((U2x + U3x)*U4x + U4xx) + 
     -       H**3*((cs + sn*(U2y + U3y))*U4y + sn*U4yy)))/(g*H)
     -   + (dExp(2.*U1 - 2.*U3)*H**4*sn**2*U4*
     -     (nr*U5 + g*w)**2)/(g*r**4)  
 

      eq5=  (-1.*dExp(2.*U2)*sn*
     -     (2.*r**4*sn*U5 + 
     -       ff**2*g*r**3*sn*(6.*U2x*U5 - 2.*U3x*U5 - 1.*U5x) + 
     -       ff**2*g**2*r*sn*U5x + 
     -       ff**3*g**2*r**2*sn*
     -        ((2. + ff*(-3.*U2x + U3x))*U5x - 1.*ff*U5xx) - 
     -       1.*g*H**2*(3.*cs*U5y + 
     -          sn*(3.*U2y*U5y - 1.*U3y*U5y + U5yy))) - 
     -    8.*dExp(2.*U1)*alpha**2*g*H**2*nr*U4**2*
     -     (nr*U5 + g*w))/(dExp(2.*U2)*g**3) 
 
      
	p(I,1)=eq1          
      p(I,2)=eq2  
      p(I,3)=eq3  
      p(I,4)=eq4  
      p(I,5)=eq5 


      if(Abs(U4) .le. max) go to 124
	max=Abs(U4)
	rm=r
124    continue

c       	write(6,*)r,eq5 

       T34= (-2.*H**2*nr*U4**2*(nr*U5 + g*w))/
     -  (dExp(2.*U3)*r**4)
	   
 
       T44=  -1.*c2*U4**4 - 1.*c1*U4**6 - 
     -  (1.*(ff**4*g*r**2*U4x**2 + H**2*U4y**2))/
     -   (dExp(2.*U1)*g*H**2) + 
     -  U4**2*(-1.*c3 - (1.*csc**2*nr**2)/
     -      (dExp(2.*U2)*g) - 
     -     (1.*g*H**2*((-1.*nr**2*U5**2)/g**2 + w**2))/
     -      (dExp(2.*U3)*r**4))   

c Ttot= 2Td[4, 4] - Sum[Td[i, i]]
         
	 Ttot=  2.*U4**2*(c3 + c2*U4**2 + 
     -    (dExp(2.*U3)*c1*r**4*U4**4 - 
     -       2.*H**2*w*(nr*U5 + g*w))/
     -     (dExp(2.*U3)*r**4))

           WRITE(1,489)T34 ,T44,Ttot

 11        CONTINUE

	        Close(Unit=1)



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

	
       write(6,*)" "
       do 1200 i=1,nk
        ptest(i)=0.d0
 	      do 1300 j=1,nv
 	        ptest(i)=ptest(i)+p(j,i)**2
 1300   continue 
        write(6,*)"res=",i, sqrt(ptest(i))
 1200 continue 
       write(6,*)" "



         rewind(1)
        open(unit=1,name='derx.dat')
      DO 110  I = 1,NV
          WRITE(1,789)  (1.-X(I))**2*UX(I,1) 
          WRITE(1,789)  (1.-X(I))**2*UX(I,2) 
          WRITE(1,789)  (1.-X(I))**2*UX(I,3) 
	    WRITE(1,789)  (1.-X(I))**2*UX(I,4) 
	    WRITE(1,789)  (1.-X(I))**2*UX(I,5) 	  		   		    
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
  
c  	write(6,*)"here"   

  789     format (e24.16)
  489     format (6(e24.16))
  765     format (8(e24.16))  
c         write(6,*)"HERE@#"                     
C
C****END OF CALCULATION
C     ------------------
C
      R E T U R N
C-----END OF FDSU01----------------------------------------------------
      E    N    D

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
C                                 X
C                                  X
C                                   X
C                                    X      
C                                     X
C                                      X
C                                       X
C                                        X
C                                         X
C                                          X
C                                           X
C                                            X      
C                                             X
C                                              X
C                                               X
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C    23.08.99 -INCEPUT
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C**                                                               ***       
C**                                                               ***
      SUBROUTINE  U02E1(IRAND,T,X,Y,U,UT,UX,UXX,UY,UYY,                         
     &                   P,MT,MV,NK,NB)                                          
C**                                                               ***       
C**                                                               ***       
C**********************************************************************
C**                                                               ***       
C**    F D S U 0 2   PROGRAM FRAME FOR THE USER-SUBROUTINE,       ***
C**    IN WHICH THE BOUNDARY-CONDITIONS ON THE 4 BOUNDARY-LINES   ***       
C**    ARE PRESCRIBED.                                            ***       
C**                                                               ***
C**********************************************************************       
C**                                                               ***       
C**    FORMAL PARAMETERS :                                        ***       
C**    -------------------                                        ***       
C**                                                               ***       
      INTEGER  IRAND,MT, MV, NK, NB 
      DOUBLE PRECISION  T,X(MV),Y(MV),P(MV,NK),U(MV,NK),UT(MT,NK),              
     &       UX(MV,NK),UXX(MV,NK),UY(MV,NK),UYY(MV,NK)               

 
	DOUBLE PRECISION val,OmegaH 
      DOUBLE PRECISION  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
      COMMON /FDSCON/  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
C**                                                               ***
C**-----------------------------------------------------------------***
C**                                                               ***
C** LIST OF FORMAL PARAMETERS : (SEE DOCUMENTATION 1.4 B)         ***
C** ---------------------------                                   ***       
C**                                                               ***
C**-----------------------------------------------------------------***
C**                                                               ***       
C**      LOCAL PARAMETERS :                                       ***       
C**      ----------------                                         ***
C**                                                               ***       
      INTEGER  I                                                                
C**                                                               ***       
C**-----------------------------------------------------------------***
C                                                                               
C                                                                               
C****START OF CALCULATION :
C     ---------------------                                                     
C                                                                               
C
C***SAME BOUNDARY-CONDITIONS FOR BOUNDARY 3 AND 4/                            
C                                                                               
c      IRAND = 1:  X = X(1) =0 (r=0)
C     IRAND = 2 : X = X(NX)=1 (r=infinity)  
C     IRAND = 3 : Y = Y(1) =0 
C     IRAND = 4 : Y = Y(NY)=Pi/2
c
C     NB = NUMBER OF GRIDPOINTS WITHIN THE BOUNDARY-LINE
c      write(6,*)'se apeleaza boundary conditions'
        
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
	    p(1,4)=u(1,4)  
	    p(1,5)=u(1,5) +w/nr*rh**2

      write(1,765)y(1),u(1,1),u(1,2),u(1,3),0.d0,u(1,5)  
      write(2,765)y(1),ux(1,1),ux(1,2),ux(1,3),ux(1,4),ux(1,5)  
      write(3,765)y(1),uxx(1,1),uxx(1,2),uxx(1,3),uxx(1,4),uxx(1,5) 
      write(4,765)y(1),uy(1,1),uy(1,2),uy(1,3),uy(1,4),uy(1,5) 
      write(5,765)y(1),uyy(1,1),uyy(1,2),uyy(1,3),uyy(1,4),uyy(1,5) 

      
       DO 11  I = 2,NB 
          p(i,1)=Ux(I,1) 
  	    p(i,2)=Ux(I,2) 
	    p(i,3)=ux(i,3)  
	    p(i,4)=u(i,4)  
	    p(i,5)=u(i,5) +w/nr*rh**2

c       	   write(6,*)p(i,4),P(I,5),u(i,5)

      write(1,765)y(I),u(i,1),u(i,2),u(i,3),0.,u(i,5)  
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
C
C
C****END OF CALCULATION
C     ------------------
C
      R E T U R N
C-----END OF FDSU02----------------------------------------------------
      E    N    D
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C    23.08.99 -INCEPUT
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C**                                                               ***
C**                                                               ***
      SUBROUTINE  U04E1(IRAND,IEQU,ICOM,T,X,Y,U,UT,UX,UXX,UY,UYY,
     &                   PU,PUT,PUX,PUXX,PUY,PUYY,MT,MV,NK,NB)
C**                                                               ***       
C**                                                               ***       
C**********************************************************************       
C**                                                               ***       
C**    F D S U 0 4       PROGRAM FRAME FOR THE USER-SUBROUTINE,   ***       
C**    IN WHICH THE JACOBIAN MATRICES ON THE 4 BOUNDARY-LINES     ***       
C**    ARE PRESCRIBED. ONLY NONZERO ELEMENTS MUST BE DEFINED.     ***       
C**                                                               ***       
C**********************************************************************
C**                                                               ***       
C**    FORMAL PARAMETERS :                                        ***
C**    -------------------                                        ***
C*IRAND   NUMBER OF THE BOUNDARY-LINES                           
C         IRAND = 1 : X = X(1)
C         IRAND = 2 : X = X(NX)                                  
C         IRAND = 3 : Y = Y(1)                                   
C         IRAND = 4 : Y = Y(NY)
C IEQU    NUMBER OF THE ACTUAL EQUATION OF THE PDE-SYSTEM  
C ICOM    NUMBER OF THE SOLUTION COMPONENT                                                                                                   ***       
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C  semnificatie:
C 
      INTEGER  IRAND, IEQU, ICOM, MT, MV, NK, NB                                
      DOUBLE PRECISION  T,X(MV),Y(MV),U(MV,NK),UT(MT,NK),UX(MV,NK),             
     &                   UXX(MV,NK),UY(MV,NK),UYY(MV,NK),PU(MV),                 
     &                   PUT(MT),PUX(MV),PUXX(MV),PUY(MV),PUYY(MV)
     
           DOUBLE PRECISION  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
      COMMON /FDSCON/  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
	              
C**                                                               ***
C**-----------------------------------------------------------------***       
C**                                                               ***       
C** LIST OF FORMAL PARAMETERS : (SEE DOCUMENTATION 1.4 D)         ***       
C** ---------------------------                                   ***       
C**                                                               ***       
C**-----------------------------------------------------------------***       
C**                                                               ***       
C**      LOCAL PARAMETERS :                                       ***       
C**      ----------------                                         ***       
C**                                                               ***       
      INTEGER  I
C**                                                               ***
C**-----------------------------------------------------------------***
C
C
C****START OF CALCULATION :                                                    
C     ---------------------                                                     
C                    
C
c       write(6,*)'se apeleaza boundary-jacobian'

c       write(6,*)"HERE+++++++++++++++"
  
      GO TO (1000,2000,3000,4000) IRAND

c boundary conditions r=0
 1000 CONTINUE
      DO 1111 i=1,Nb

c    

c          write(6,*)rh,par

 	 IF(IEQU.EQ.1.AND.ICOM.EQ.1) pux(i)=1. 
 	 IF(IEQU.EQ.2.AND.ICOM.EQ.2) pux(i)=1.    
       IF(IEQU.EQ.3.AND.ICOM.EQ.3) pux(i)=1.
       IF(IEQU.EQ.4.AND.ICOM.EQ.4) pu(i)=1.  
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
C-----END OF FDSU04----------------------------------------------------         
      E    N    D
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C    23.08.99 -INCEPUT
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C**                                                             ***
C**                                                             ***     
      SUBROUTINE  U03E1(IEQU,ICOM,T,X,Y,U,UT,UX,UXX,UY,UYY,
     &               PU,PUT,PUX,PUXX,PUY,PUYY,MT,MV,NK,NV)                   
C**                                                             ***
C**                                                             ***     
C**********************************************************************     
C**                                                             ***     
C**  F D S U 0 3       PROGRAM FRAME FOR THE USER-SUBROUTINE,   ***     
C**  IN WHICH THE JACOBIAN MATRICES AT INNER GRIDPOINTS ARE     ***     
C**  PRESCRIBED. ONLY NONZERO ELEMENTS MUST BE DEFINED.         ***
C**                                                             ***     
C**********************************************************************
C**                                                             ***
C**  FORMAL PARAMETERS :                                        ***     
C**  -------------------                                        ***     
C**                                                             ***     
      INTEGER  IEQU, ICOM, MT, MV, NK, NV                                       
      DOUBLE PRECISION  T,X(MV),Y(MV),U(MV,NK),UT(MT,NK),UX(MV,NK),
     &               UXX(MV,NK),UY(MV,NK),UYY(MV,NK),PU(MV), puxy(mv),                 
     &               PUT(MT),PUX(MV),PUXX(MV),PUY(MV),PUYY(MV),max 
       Double precision sn,cs,csc,ct,r,z,ff,f

       Double precision 
     & U1,U2,U3,
     & U1x,U2x,U3x, 
     & U1y,U2y,U3y 

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

C**                                                             ***
C**-----------------------------------------------------------------***
C**                                                             ***
C** LIST OF FORMAL PARAMETERS :  (SEE DOCUMENTATION 1.4 C)        ***
C** ---------------------------                                   ***
C**                                                             ***
C**-----------------------------------------------------------------***
C**                                                             ***     
C**    LOCAL PARAMETERS :                                       ***
C**    ----------------                                         ***
C**                                                             ***
      INTEGER  I
      DOUBLE PRECISION   H,derH,der2H,g 
      DOUBLE PRECISION  nr,w,con,alpha,c1,c2,c3,lambda,rh,par
      COMMON /FDSCON/  nr,w,con,alpha,c1,c2,c3,lambda,rh,par                          
C**                                                             ***
C**-----------------------------------------------------------------***
C
C
C****START OF CALCULATION :
C     ---------------------
C
C
C     IEQU  NUMBER OF THE ACTUAL EQUATION OF THE PDE-SYSTEM
c     ICOM  NUMBER OF THE SOLUTION COMPONENT
             
         


c	write(6,*)eta,nr,w,con,g,lambda1,lambda2,eta1,const

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
							
	U1xx=UXX(I,1)
      U2xx=UXX(I,2) 
      U3xx=UXX(I,3) 
      U4xx=UXX(I,4) 
      U5xx=UXX(I,5) 

       r= X(I)/(1.-X(I)) 

         sn=dSin(Y(I))
         cs=dCos(Y(I))
	csc=1.d0/dSin(Y(I))
         ct=dCos(Y(I))/dSin(Y(I))
	 sn2=dSin(2.*Y(I))
 	cs2=dCos(2.*Y(I))
 

	ff=1.-X(I) 
       g=r**2 + rh**2 

       H= rh + Sqrt(g)

	 derH= r/Sqrt(g)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       Vxx11=ff**4*r**2*sn**2
       Vxy11=0.
       Vyy11=(H**2*sn**2)/g
       Vx11=
     & (ff**2*r*(H*r**2 + g*(H -derH*r-2.*ff*H*r))*sn**2)/(g*H)
     &
       Vy11=0.
       V11=
     & (4.*dExp(2.*U1)*alpha**2*H**2*U4**2*
     &     ((-1.*nr**2)/dExp(2.*U2) + 
     &       (H**2*sn**2*(nr*U5 + g*w)**2)/
     &        (dExp(2.*U3)*r**4)))/g
     &
       Vxx12=0.
       Vxy12=0.
       Vyy12=0.
       Vx12=
     & (ff**2*r*sn**2*(H*r**2 + 
     &       g*(derH*r - 1.*H*(2. + ff**2*r*U3x))))/(g*H)
     &
       Vy12=(-1.*H**2*sn**2*U3y)/g
       V12=
     & (0.5*H**2*(8.*dExp(2.*U1)*alpha**2*g*nr**2*
     &        U4**2 - (1.*dExp(4.*U2 - 2.*U3)*sn**4*
     &          ((-2.*r*U5 + ff**2*g*U5x)**2 + (g*H**2*U5y**2)/r**2)
     &)/r**2))/(dExp(2.*U2)*g**2)
     &
       Vxx13=0.
       Vxy13=0.
       Vyy13=0.
       Vx13=(-1.*ff**2*r**2*sn**2*(r + ff**2*g*U2x))/g
       Vy13=(-1.*H**2*sn*(cs + sn*U2y))/g
       V13=
     & (0.5*H**2*sn**2*(dExp(2.*U2)*sn**2*
     &        (4.*r**4*U5**2 - 4.*ff**2*g*r**3*U5*U5x + 
     &          ff**4*g**2*r**2*U5x**2 + g*H**2*U5y**2) - 
     &       8.*dExp(2.*U1)*alpha**2*g*H**2*U4**2*
     &        (nr*U5 + g*w)**2))/
     &   (dExp(2.*U3)*g**2*r**4)
     &
       Vxx14=0.
       Vxy14=0.
       Vyy14=0.
       Vx14=4.*alpha**2*ff**4*r**2*sn**2*U4x
       Vy14=(4.*alpha**2*H**2*sn**2*U4y)/g
       V14=
     & (4.*dExp(2.*U1)*alpha**2*H**2*U4*
     &     ((-1.*nr**2)/dExp(2.*U2) + 
     &       (H**2*sn**2*(nr*U5 + g*w)**2)/
     &        (dExp(2.*U3)*r**4)))/g
     &
       Vxx15=0.
       Vxy15=0.
       Vyy15=0.
       Vx15=
     & (-0.5*dExp(2.*U2 - 2.*U3)*ff**2*H**2*sn**4*
     &     (-2.*r*U5 + ff**2*g*U5x))/(g*r**2)
     &
       Vy15=
     & (-0.5*dExp(2.*U2 - 2.*U3)*H**4*sn**4*U5y)/
     &   (g*r**4)
     &
       V15=
     & (H**2*sn**2*(-1.*dExp(2.*U2)*r**3*sn**2*
     &        (2.*r*U5 - 1.*ff**2*g*U5x) + 
     &       4.*dExp(2.*U1)*alpha**2*g*H**2*nr*U4**2*
     &        (nr*U5 + g*w)))/(dExp(2.*U3)*g**2*r**4)
     &
       Vxx21=0.
       Vxy21=0.
       Vyy21=0.
       Vx21=0.
       Vy21=0.
       V21=
     & (4.*dExp(2.*U1 - 2.*U2)*alpha**2*H**2*U4**2*
     &     (2.*nr**2 + dExp(2.*U2)*g*sn**2*
     &        (c3 + c2*U4**2 + c1*U4**4)))/g
     &
       Vxx22=ff**4*r**2*sn**2
       Vxy22=0.
       Vyy22=(H**2*sn**2)/g
       Vx22=
     & (ff**2*r*sn**2*(2.*H*r**2 + 
     &       g*(-2.*derH*r + H*(3. - 2.*ff*r + 
     &             ff**2*r*(2.*U2x + U3x)))))/(g*H)
     &
       Vy22=(H**2*sn*(2.*cs + sn*(2.*U2y + U3y)))/g
       V22=
     & (H**2*(-8.*dExp(2.*U1)*alpha**2*g*nr**2*
     &        U4**2 + (dExp(4.*U2 - 2.*U3)*sn**4*
     &          ((-2.*r*U5 + ff**2*g*U5x)**2 + (g*H**2*U5y**2)/r**2)
     &)/r**2))/(dExp(2.*U2)*g**2)
     &
       Vxx23=0.
       Vxy23=0.
       Vyy23=0.
       Vx23=(ff**2*r**2*sn**2*(r + ff**2*g*U2x))/g
       Vy23=(H**2*sn*(cs + sn*U2y))/g
       V23=
     & (-1.*dExp(2.*U2 - 2.*U3)*H**2*sn**4*
     &     ((-2.*r*U5 + ff**2*g*U5x)**2 + (g*H**2*U5y**2)/r**2))/
     &   (g**2*r**2)
     &
       Vxx24=0.
       Vxy24=0.
       Vyy24=0.
       Vx24=0.
       Vy24=0.
       V24=
     & (4.*dExp(2.*U1 - 2.*U2)*alpha**2*H**2*U4*
     &     (2.*nr**2 + dExp(2.*U2)*g*sn**2*
     &        (c3 + 2.*c2*U4**2 + 3.*c1*U4**4)))/g
     &
       Vxx25=0.
       Vxy25=0.
       Vyy25=0.
       Vx25=
     & (dExp(2.*U2 - 2.*U3)*ff**2*H**2*sn**4*
     &     (-2.*r*U5 + ff**2*g*U5x))/(g*r**2)
     &
       Vy25=
     & (dExp(2.*U2 - 2.*U3)*H**4*sn**4*U5y)/(g*r**4)
     &
       V25=
     & (-2.*dExp(2.*U2 - 2.*U3)*H**2*sn**4*
     &     (-2.*r*U5 + ff**2*g*U5x))/(g**2*r)
     &
       Vxx31=0.
       Vxy31=0.
       Vyy31=0.
       Vx31=0.
       Vy31=0.
       V31=
     & 4.*dExp(2.*U1)*alpha**2*H**2*sn**2*U4**2*
     &   (c3 + c2*U4**2 + c1*U4**4 - 
     &     (2.*g*H**2*((nr*U5)/g + w)**2)/
     &      (dExp(2.*U3)*r**4))
     &
       Vxx32=0.
       Vxy32=0.
       Vyy32=0.
       Vx32=
     & (ff**2*r*sn**2*(-1.*H*r**2 + 
     &       g*(-1.*derH*r + H*(2. + ff**2*r*U3x))))/(g*H)
     &
       Vy32=(H**2*sn**2*U3y)/g
       V32=
     & (-1.*dExp(2.*U2 - 2.*U3)*H**2*sn**4*
     &     ((-2.*r*U5 + ff**2*g*U5x)**2 + (g*H**2*U5y**2)/r**2))/
     &   (g**2*r**2)
     &
       Vxx33=ff**4*r**2*sn**2
       Vxy33=0.
       Vyy33=(H**2*sn**2)/g
       Vx33=
     & (ff**2*r*sn**2*(-3.*derH*r + 
     &       H*(5. - 2.*ff*r + ff**2*r*(U2x + 2.*U3x))))/H
     &
       Vy33=(H**2*sn*(cs + sn*(U2y + 2.*U3y)))/g
       V33=
     & (H**2*sn**2*(dExp(2.*U2)*sn**2*
     &        (4.*r**4*U5**2 - 4.*ff**2*g*r**3*U5*U5x + 
     &          ff**4*g**2*r**2*U5x**2 + g*H**2*U5y**2) + 
     &       8.*dExp(2.*U1)*alpha**2*g*H**2*U4**2*
     &        (nr*U5 + g*w)**2))/
     &   (dExp(2.*U3)*g**2*r**4)
     &
       Vxx34=0.
       Vxy34=0.
       Vyy34=0.
       Vx34=0.
       Vy34=0.
       V34=
     & (4.*dExp(2.*U1 - 2.*U3)*alpha**2*H**2*sn**2*U4*
     &     (dExp(2.*U3)*c3*g*r**4 + 
     &       2.*dExp(2.*U3)*c2*g*r**4*U4**2 + 
     &       3.*dExp(2.*U3)*c1*g*r**4*U4**4 - 
     &       2.*H**2*nr**2*U5**2 - 4.*g*H**2*nr*U5*w - 
     &       2.*g**2*H**2*w**2))/(g*r**4)
     &
       Vxx35=0.
       Vxy35=0.
       Vyy35=0.
       Vx35=
     & (-1.*dExp(2.*U2 - 2.*U3)*ff**2*H**2*sn**4*
     &     (-2.*r*U5 + ff**2*g*U5x))/(g*r**2)
     &
       Vy35=
     & (-1.*dExp(2.*U2 - 2.*U3)*H**4*sn**4*U5y)/
     &   (g*r**4)
     &
       V35=
     & (-2.*H**2*sn**2*(dExp(2.*U2)*r**3*sn**2*
     &        (2.*r*U5 - 1.*ff**2*g*U5x) + 
     &       4.*dExp(2.*U1)*alpha**2*g*H**2*nr*U4**2*
     &        (nr*U5 + g*w)))/(dExp(2.*U3)*g**2*r**4)
     &
       Vxx41=0.
       Vxy41=0.
       Vyy41=0.
       Vx41=0.
       Vy41=0.
       V41=
     & 2.*dExp(2.*U1)*H**2*U4*
     &   ((-1.*nr**2)/(dExp(2.*U2)*g) - 
     &     1.*sn**2*(c3 + 2.*c2*U4**2 + 3.*c1*U4**4) + 
     &     (g*H**2*sn**2*((nr*U5)/g + w)**2)/
     &      (dExp(2.*U3)*r**4))
     &
       Vxx42=0.
       Vxy42=0.
       Vyy42=0.
       Vx42=ff**4*r**2*sn**2*U4x
       Vy42=(H**2*sn**2*U4y)/g
       V42=(2.*dExp(2.*U1 - 2.*U2)*H**2*nr**2*U4)/g
       Vxx43=0.
       Vxy43=0.
       Vyy43=0.
       Vx43=ff**4*r**2*sn**2*U4x
       Vy43=(H**2*sn**2*U4y)/g
       V43=
     & (-2.*dExp(2.*U1 - 2.*U3)*g*H**4*sn**2*U4*
     &     ((nr*U5)/g + w)**2)/r**4
     &
       Vxx44=ff**4*r**2*sn**2
       Vxy44=0.
       Vyy44=(H**2*sn**2)/g
       Vx44=
     & (ff**2*r*sn**2*(H*r**2 + 
     &       g*(-2.*derH*r + H*(3. - 2.*ff*r + ff**2*r*(U2x + U3x))))
     &)/(g*H)
     &
       Vy44=(H**2*sn*(cs + sn*(U2y + U3y)))/g
       V44=
     & (-1.*dExp(2.*U1 - 2.*U2 - 2.*U3)*H**2*
     &     (dExp(2.*U3)*nr**2*r**4 + 
     &       dExp(2.*(U2 + U3))*g*r**4*sn**2*
     &        (c3 + 6.*c2*U4**2 + 15.*c1*U4**4) - 
     &       1.*dExp(2.*U2)*H**2*sn**2*
     &        (nr*U5 + g*w)**2))/(g*r**4)
     &
       Vxx45=0.
       Vxy45=0.
       Vyy45=0.
       Vx45=0.
       Vy45=0.
       V45=
     & (2.*dExp(2.*U1 - 2.*U3)*H**4*nr*sn**2*U4*
     &     ((nr*U5)/g + w))/r**4
     &
       Vxx51=0.
       Vxy51=0.
       Vyy51=0.
       Vx51=0.
       Vy51=0.
       V51=
     & (-16.*dExp(2.*U1 - 2.*U2)*alpha**2*H**2*nr*
     &     U4**2*(nr*U5 + g*w))/g**2
     &
       Vxx52=0.
       Vxy52=0.
       Vyy52=0.
       Vx52=(3.*ff**2*r**2*sn**2*(-2.*r*U5 + ff**2*g*U5x))/g**2
       Vy52=(3.*H**2*sn**2*U5y)/g**2
       V52=
     & (16.*dExp(2.*U1 - 2.*U2)*alpha**2*H**2*nr*
     &     U4**2*(nr*U5 + g*w))/g**2
     &
       Vxx53=0.
       Vxy53=0.
       Vyy53=0.
       Vx53=(-1.*ff**2*r**2*sn**2*(-2.*r*U5 + ff**2*g*U5x))/g**2
       Vy53=(-1.*H**2*sn**2*U5y)/g**2
       V53=0.
       Vxx54=0.
       Vxy54=0.
       Vyy54=0.
       Vx54=0.
       Vy54=0.
       V54=
     & (-16.*dExp(2.*U1 - 2.*U2)*alpha**2*H**2*nr*U4*
     &     (nr*U5 + g*w))/g**2
     &
       Vxx55=(ff**4*r**2*sn**2)/g
       Vxy55=0.
       Vyy55=(H**2*sn**2)/g**2
       Vx55=
     & (ff**2*r*sn**2*(r**2 + g*
     &        (-1. - 2.*ff*r + ff**2*r*(3.*U2x - 1.*U3x))))/g**2
     &
       Vy55=(H**2*sn*(3.*cs + 3.*sn*U2y - 1.*sn*U3y))/g**2
       V55=
     & (-2.*(dExp(2.*U2)*r**3*sn**2*
     &        (r + ff**2*g*(3.*U2x - 1.*U3x)) + 
     &       4.*dExp(2.*U1)*alpha**2*g*H**2*nr**2*
     &        U4**2))/(dExp(2.*U2)*g**3)
   


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
  
  
C                                                                               
C****END OF CALCULATION                                                        
C     ------------------                                                        
C                                                                               

        RETURN
                                                                    
C-----END OF FDSU03----------------------------------------------------
        
      E    N    D                                                               

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

