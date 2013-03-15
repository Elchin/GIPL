
! Notes: read new thermo and then put it in the old structure 
! File thermo and input1. txt had been corrected
! no oldId, Csnow and nbound declared inside of the code
! line 347 YYS(i,1)=0!YM(1) is not correct
! spinup added corrected version of perm4
! modified 11/03/09 :      gt_zone_code(:) was added
! modified 11/06/09
! modified 11/18/09 : SnowFix is added ( no snow in fall when temp>0)
! modified 11/21/09 : new reading procedure added ak_gtzone_18km instead of initial.txt
! modified 11/23/09 : reads xt(1) from bound; first day of calculations and makes tini=xt(1)
! modified 02/14/10
! modified 04/24/10 : nmean is corrected for 366 case, allocated(numsl) removed
! modified 09/04/10 : modified for correct nmean
! modified 09/15/10 : added vegetation file read, 1 pros initializes all the necess input variables
! modified 06/20/11 : getrid of auxilary function references
! modified 06/27/11 : fully serial
! modified 07/06/11 : fully serial
! NOTE!!! MODIFIE START READING ROUTINNE
      MODULE BND
	INTEGER NXT
	real*8,allocatable::  XT(:)
	real*8,allocatable::  YT(:,:)

      	real*8,allocatable::  XTOUT(:)
	real*8 ,allocatable:: YTOUT(:,:)
      	real*8 TINIR,TINI
	real*8  ,allocatable:: XXSNOW1(:),XXSNOW2(:,:)
	real*8 ,allocatable:: SNOWYL1(:),SNOWYL2(:,:)
	real*8 ,allocatable:: YSNOUT (:,:)
	real*8 ,allocatable:: RSNOUT (:,:)
      	INTEGER nxxsnow,NLsnow
	integer lbound
	
    	integer :: ist      !cmd
	real*8 	:: STEP,TAUM,TMIN
	real*8  :: itfin,itfinl
	real*8  :: E1,UWK
	real*8  :: perl!,SLEV
	integer :: itmax,NMEAN,NFRMAX
	real*8  :: YFMAX,YFMIN,WFAZ

      END MODULE BND

      MODULE THERMO
      	real*8 QPHASE,SLEV
	real*8,allocatable:: VARNOD(:,:),ACLV(:,:),BCLV(:,:),TFPV(:,:),EE(:,:)
	real*8,allocatable:: CAPFR(:,:),CAPTH(:,:),XKFR(:,:),XKTH(:,:)

	real*8,allocatable:: VARNOD1(:,:),ACLV1(:,:),BCLV1(:,:),EE1(:,:)
	real*8,allocatable:: CAPFR1(:,:),CAPTH1(:,:),XKFR1(:,:),XKTH1(:,:)
	real*8 CSNOW

	real*8, allocatable :: U(:,:)
	real, allocatable:: YYS(:,:)
	integer k0
	
        integer, allocatable :: snow_code(:),veg_code(:)     !input params
	integer, allocatable :: geo_code(:),gt_zone_code(:)
	real*8, allocatable  :: GRAD(:)

	real*8 ,allocatable:: RES(:,:)

	
      END MODULE THERMO

      MODULE GRDG
	INTEGER NX,NY,NMY
      	real*8,allocatable:: YM(:),HY(:)
	INTEGER,allocatable:: NLX(:,:),MY(:),NUMSL(:)
	real*8, allocatable :: XTI(:),YTI(:,:)!initial
	integer :: nxti,nnpg
      END MODULE GRDG

      MODULE AL
      	integer,allocatable::NFRNT(:,:)
	integer,allocatable::IMEAN(:)
	real*8 ,allocatable::YFRNT(:,:,:)
      END MODULE AL


      PROGRAM permgis
	USE BND
      	USE THERMO
      	USE GRDG
	USE AL
      implicit none
	
      	EXTERNAL BOUND,GLKY,SUNWATER,GLC

      	real*8 UNWATER,DUNWATER,GLCQ,GLKY,GLC,BOUND,SLEVEL
	real*8 RSNOW,SUNWATER

      	integer I,II,J,JJ,nout
      	integer NMEANS,NNN

      	real*8 ,allocatable:: RESTMP(:,:)
      	real*8, allocatable:: TTT(:),tlet(:)
      	real*8 ,allocatable:: YYFRNT(:)

      	real*8 tfin,tfinl
      	real tau1,FREEZUPD,FREEZUP

        character(64)  :: startf,resultf1,MEANF1
	character(210) :: FMT1,FMT2
! MPI varriables----------------------------------------------
      integer   my_rank  ! My process rank.
      integer   mype     ! The number of processes.
      integer   nx1,nx2  ! partitioning with dimensions nx2-nx1
      integer   subset ! NX/mype
      integer   ierr
!------------------------------------------------------------


!my_rank=3
!mype=4
!	OPEN(60,FILE='process.txt')
!   	  READ(60,*)my_rank !process number 
!          READ(60,*)mype    !number processors used
!	CLOSE(60)
	my_rank=0
	mype=1

	call initialize


	
	subset = nx/mype      ! subset length
	nx1 = my_rank*subset+1
	nx2 = (my_rank+1)*subset 
	if ((my_rank+1).eq.mype) nx2=(my_rank+1-mype)*subset + nx

	allocate(RESTMP(NMY+3,NX),STAT=IERR)
      	allocate(YYFRNT(NMEAN),STAT=IERR)
      	allocate(TTT(NX),STAT=IERR)
      	allocate(TLET(NX),STAT=IERR)
        allocate(RES(NMEAN,NMY+3),STAT=IERR)
!________________INITIALIZATION________________________

!      NOUT=NMEAN*ITFINL+3
	NOUT=NMEAN+2
      	TFIN=STEP*DBLE(NMEAN*ITFIN)
	TFINL=STEP*DBLE(NMEAN*ITFINL)
	TTT=0.0D0
	TINIR=0.0D0
!-----------------------------------------------initial data read
	call init_cond(ist,my_rank,nx1,nx2)
	
	write(MEANF1,'(A10,I3.3,A4)') 'dump/mean_',my_rank+1,'.txt'
	write(resultf1,'(A12,I3.3,A4)') 'dump/result_',my_rank+1,'.txt'
	write(startf,'(A11,I3.3,A4)') 'dump/start_',my_rank+1,'.txt'
	open(1,file=resultf1,STATUS='unknown')
	open(2,file=MEANF1,STATUS='unknown')
      	open(3,file=startf,STATUS='unknown')
	write(FMT1,'(A30,I0,A12)')'(1x,I10,1x,F12.3,2(1x,F16.12),',NMY,'(1x,F16.12))'
	write(FMT2,'(A28,I0,A40)')'(1x,I10,1x,F12.3,2(1x,F8.3),',NMY,  '(1x,F8.3),(1x,F8.3,1x,F12.3),(1x,F12.3))'

      	allocate(XTOUT(NOUT),STAT=IERR)
	allocate(YTOUT(NOUT,NX),STAT=IERR)
	allocate(YSNOUT(NOUT,NX),STAT=IERR)
	allocate(RSNOUT(NOUT,NX),STAT=IERR)
	
!****************************************************************************	
      do I=1,NMEAN+2
	 XTOUT(I)=TINI+DBLE(I-1)*STEP
      enddo
      IMEAN=1

      do I=nx1,nx2
	 if (lbound.EQ.2)GRAD(i)=GRAD(i)*ym(ny) 
         do ii=1,NUMSL(i)
	      TFPV(ii,i)=-(varnod(ii,i)/aclv(ii,i))**(1.d0/bclv(ii,i))
	 enddo
	 call LININTRP(XT,YT(:,I),NXT,XTOUT,YTOUT(:,I),NOUT)
      	 call LININTRP(XXSNOW1,XXSNOW2(:,I),NXXSNOW,XTOUT,YSNOUT(:,I),NOUT)
	 call SnowFix(YTOUT(:,I),YSNOUT(:,I),NOUT)
      	 call LININTRP(SNOWYL1,SNOWYL2(:,I),NLSNOW,XTOUT,RSNOUT(:,I),NOUT)
	 call active_layer(I)
      enddo

65    CONTINUE
      do I=nx1,nx2
          TLET(I)=TTT(I)+TINI
	  call save_results(I,TLET(I),TTT(I))
          6666  continue

	  call stefan1D(U(I,:),NY,HY,TTT(I),TAUM,TMIN,STEP,I,NLX(I,:) &
	    ,ITMAX,E1,UWK,lbound,GRAD(I),BOUND,GLC,GLKY,SUNWATER)
	
      !--------------------------------------------
      !--------------------------------------------
!	    do j=1,ny			! WRITTING RESULTS
!	       write(1,'3(f12.3)') U(I,j),ym(j),GLKY(U(I,j),I,J,tlet(I))
!	    enddo
	
	
	
	  TTT(I)=TTT(I)+STEP
          TLET(I)=TTT(I)+TINI
          if(IMEAN(I).LT.NMEAN)  then
            IMEAN(I)=IMEAN(I)+1
	    call save_results(I,TLET(I),TTT(I))
	    call active_layer(I)

            GOTO 6666
          endif
          if(TFIN.LT.TFINL.AND.TTT(nx1).GT.TFIN)then
	    do II=1,NMEAN			! WRITTING RESULTS
	       write(1,FMT1) I, (RES(II,JJ),JJ=1,NMY+3)
	    enddo
	  endif
      	  do jj=1,nmy+3
		  RESTMP(jj,i)=sum((res(:,JJ)))
	  enddo
      enddo
      
      IMEAN=1
      do I=nx1,nx2
           FREEZUP=-7777.D0
	   FREEZUPD=FREEZUP
           do J=2,NMEAN
              if((NFRNT(J,I)-NFRNT(J-1,I)).EQ.-2)then
        	if(YFRNT(J-1,NFRNT(J-1,I),I).GE.YFMIN) FREEZUP=SNGL(RES(J,1))
              endif
           enddo
	   if(FREEZUP.GT.0.0)then
	       FREEZUPD=AMOD(FREEZUP,REAL(NMEAN))
	       if(FREEZUPD.EQ.0.0)FREEZUPD=REAL(NMEAN)
           endif
           NMEANS=NMEAN
           YYFRNT=YFRNT(:,1,I)

  	   call save_results(I,TLET(I),TTT(I))

	   call active_layer(I)

	   !____WRITTING MEAN
	   write(2,FMT2) I,(RESTMP(JJ,I)/DBLE(NMEAN),JJ=1,NMY+3), &
 		   YYFRNT(NMEAN),FREEZUP,FREEZUPD	
	   do j=1,NOUT
	     XTOUT(j)=TLET(nx1)+DBLE(j-1)*STEP
	   enddo
      	   call LININTRP(XT,YT(:,I),NXT,XTOUT,YTOUT(:,I),NOUT)
      	   call LININTRP(XXSNOW1,XXSNOW2(:,I),NXXSNOW,XTOUT,YSNOUT(:,I),NOUT)
	   call SnowFix(YTOUT(:,I),YSNOUT(:,I),NOUT)
	   call LININTRP(SNOWYL1,SNOWYL2(:,I),NLSNOW,XTOUT,RSNOUT(:,I),NOUT)
      enddo

      rewind(3) ! -------------start file writting begin     
      write(3, * ) tini
!      write(3, * ) TLET(1)
      do J=1,NY
         write (3,* ) ( U(II,J),II=nx1,nx2)
      enddo     ! -------------start file writting end     

      TINIR=TTT(nx1)
!************************************************************
6444    if(TTT(nx1).LT.TFINL)GO TO 65
	close(1);close(2);close(3)
!   call MPI_FINALIZE(ierr)

END

subroutine initialize
	USE BND
      	USE THERMO
      	USE GRDG
	USE AL
implicit none 

	integer IREAD,ierr
	integer :: i,j,k,z_num

	real*8, allocatable:: Y(:)
      	real*8 glm,hcscale
	
	real*8 ,allocatable ::gtzone(:,:)
  	character*64 stdummy

        character*64 finput,boundf,snowf,rsnowf,inif
        character*64 cmdf,cetkaf,fvegetation,fgeology

	real*8,allocatable:: A1(:,:),A2(:,:),A3(:,:),A4(:,:),A5(:,:)
	real*8,allocatable:: A6(:,:),A7(:,:),A8(:,:),A9(:,:),A10(:,:)
        integer, allocatable :: veg_class(:), num_vl(:)
	integer :: vln

	real*8,allocatable:: B1(:,:),B2(:,:),B3(:,:),B4(:,:),B5(:,:)
	real*8,allocatable:: B6(:,:),B7(:,:),B8(:,:)
        integer, allocatable :: geo_class(:), num_gl(:)
	real*8 :: layer_thick
	integer :: gln


      finput='in/input1.txt'
      boundf='in/bound.txt'
      snowf='in/snow.txt'
      rsnowf='in/rsnow.txt'
      inif='in/initial.txt'      
      cmdf='in/cmd.txt'      
      cetkaf='in/grid.txt'
      fvegetation='in/vegetation.txt'            
      fgeology='in/geo.txt'
      
      	call filexist(finput)
	call filexist(boundf)
	call filexist(snowf)
	call filexist(rsnowf)
	call filexist(cetkaf)
	call filexist(inif)
	call filexist(fgeology)
	call filexist(fvegetation)
	call filexist(cmdf)

      open(60,FILE=finput)
   	read(60,*)NX
        allocate(snow_code(nx),STAT=IERR)
        allocate(veg_code(nx),STAT=IERR)
        allocate(geo_code(nx),STAT=IERR)
	allocate(gt_zone_code(nx),STAT=IERR)
        allocate(GRAD(nx),STAT=IERR)
        do i=1,nx
          read(60,*) IREAD,snow_code(i),veg_code(i),geo_code(i),&
	  gt_zone_code(i),GRAD(i)
        enddo
      close(60)
!      print*, trim(finput),' has been read'
      
      open(60,file=BOUNDF)
      	read(60,*)nxt
	allocate(XT(NXT),STAT=IERR)
	allocate(YT(NXT,NX),STAT=IERR)
	do i=1,nxt
          read(60,*) XT(I),(YT(I,J),J=1,NX)
	enddo
      close(60)
!      print*,trim(boundf),' has been read'

      open(60,file=RSNOWF)
         read(60,*)NLsnow
	allocate(SNOWYL1(NLSNOW),STAT=IERR)
	allocate(SNOWYL2(NLSNOW,NX),STAT=IERR)
	do i=1,nlsnow
      	   read(60,*) snowyl1(i),(snowyl2(i,J),J=1,NX)
	enddo
      close(60)
!      print*,trim(rsnowf),' has been read'

      open(60,file=SNOWF)
      	read(60,*)nxxsnow
	allocate(XXSNOW1(NXXSNOW),STAT=IERR)
	allocate(XXSNOW2(NXXSNOW,NX),STAT=IERR)
	do I=1,nxxsnow
          read(60,*) xxsnow1(i),(xxsnow2(i,J),J=1,NX)
	enddo
      close(60)
!      print*,trim(snowf),' has been read' 

      open(60,file=inif,action='read')
        read(60,*)z_num,nxti!,TINI
	allocate(XTI(NXTI),STAT=IERR)
	allocate(YTI(NXTI,NX),STAT=IERR)
        allocate(gtzone(NXTI,z_num+1),STAT=IERR)
   	read(60,*)stdummy
        do i=1,nxti
	  read(60,*) (gtzone(i,j),j=1,z_num+1)
        enddo
      close(60)
!      print*,trim(inif),'has been read'	
      
      TINI=xt(1)
      XTI(:)=gtzone(:,1)
      do i=1,nx
	k=gt_zone_code(i)
	YTI(:,I)=gtzone(:,k+1)
      enddo

      open(60,FILE=cmdf)
	read(60,*)IST
	read(60,*)STEP,TAUM,TMIN
      	read(60,*)ITFIN,ITFINL
      	read(60,*)E1,UWK,itmax
      	read(60,*)perl,NMEAN
	read(60,*)SLEV,NFRMAX
        read(60,*)YFMIN,YFMAX
      	read(60,*)WFAZ
      close(60)
!      print*,trim(cmdf),' has been read'


      open(60,file=cetkaf)
	read(60,*)NY
	allocate(YM(NY),STAT=IERR)
      	do i=1,NY
         read(60,*) YM(i)
      	enddo
	Read(60,*)NMY
	allocate(MY(NMY),STAT=IERR)
	do j=1,NMY
 	   Read(60,*)MY(j)
	enddo
	close(60)
!      print*,trim(cetkaf),' has been read'

! note: that all max NUMSL layers has to read ro it will a give segmantation error
      NNPG=10!MAXVAL(NUMSL)      
!----------------------------------------------------      
      open (60, file=fvegetation)
	read(60,*) vln ! reads numbers of  classes
	allocate(A1(NNPG,vln),STAT=IERR) ! varnod
	allocate(A2(NNPG,vln),STAT=IERR) ! aclv
	allocate(A3(NNPG,vln),STAT=IERR) ! bclv
	allocate(A4(NNPG,vln),STAT=IERR) ! capfr
	allocate(A5(NNPG,vln),STAT=IERR)  !capth 
	allocate(A6(NNPG,vln),STAT=IERR)  !xkfr
	allocate(A7(NNPG,vln),STAT=IERR)  !xkth
	allocate(A8(vln,NNPG),STAT=IERR)  !bot_cond
	allocate(veg_class(vln),STAT=IERR) !veg_class
	allocate(num_vl(vln),STAT=IERR)  !num_vl number of vegetation layers
	
	do I = 1,vln
      	  read(60,*)veg_class(i),num_vl(i)
 	  do j=1,num_vl(i)
	  read(60,*)A1(J,I),A2(J,I),A3(J,I), &
		  A4(J,I),A5(J,I),A6(J,I),A7(J,I),A8(I,J)
          enddo
	enddo
      close(60)
!      print*,trim(fvegetation),' has been read'
      
      open (60, file='in/geo.txt')
	read(60,*) gln ! reads numbers of  classes
	allocate(B1(NNPG,gln),STAT=IERR) ! varnod
	allocate(B2(NNPG,gln),STAT=IERR) ! aclv
	allocate(B3(NNPG,gln),STAT=IERR) ! bclv
	allocate(B4(NNPG,gln),STAT=IERR) ! capfr
	allocate(B5(NNPG,gln),STAT=IERR)  !capth 
	allocate(B6(NNPG,gln),STAT=IERR)  !xkfr
	allocate(B7(NNPG,gln),STAT=IERR)  !xkth
	allocate(B8(gln,NNPG),STAT=IERR)  !bot_cond
	allocate(geo_class(gln),STAT=IERR) !geo_class
	allocate(num_gl(gln),STAT=IERR)  !num_vl number of lithologic layers
	do I = 1,gln
      	  read(60,*)geo_class(i),num_gl(i)
 	  do j=1,num_gl(i)
  		read(60,*)B1(J,I),B2(J,I),B3(J,I), &
		  B4(J,I),B5(J,I),B6(J,I),B7(J,I),B8(I,J)
          enddo
	enddo
      close(60)
!      print*,trim(fgeology),' has been read'


      allocate(VARNOD(NNPG,nx),STAT=IERR)
      allocate(ACLV(NNPG,nx),STAT=IERR)
      allocate(BCLV(NNPG,nx),STAT=IERR)
      allocate(EE(NNPG,nx),STAT=IERR)
      allocate(CAPFR(NNPG,nx),STAT=IERR)
      allocate(CAPTH(NNPG,nx),STAT=IERR)
      allocate(XKFR(NNPG,nx),STAT=IERR)
      allocate(XKTH(NNPG,nx),STAT=IERR)
      allocate(numsl(nx),STAT=IERR)
      allocate(YYS(nx,NNPG+1),STAT=IERR)

      do i = 1,nx
	layer_thick=0
	YYS(i,1)=layer_thick
	  	layer_thick=0
	YYS(i,1)=layer_thick
 	do j=1,num_vl(veg_code(i))
	   VARNOD(J,I)=A1(j,veg_code(i));
	   ACLV(J,I)=A2(j,veg_code(i));
	   BCLV(J,I)=A3(j,veg_code(i));
	   CAPTH(J,I)=A4(j,veg_code(i));
	   CAPFR(J,I)=A5(j,veg_code(i));
	   XKTH(J,I)=A6(j,veg_code(i));
	   XKFR(J,I)=A7(j,veg_code(i));
	   if (j.eq.1) then 
	      layer_thick=A8(veg_code(i),j)
	   else
	      layer_thick=layer_thick+A8(veg_code(i),j);
	   endif
	   YYS(i,j+1)=layer_thick
	   EE(J,I)=0
!	     write(*,'(3(f8.3),2(f12.1),3(f8.3))') VARNOD(J,I),ACLV(J,I),BCLV(J,I), &
!		  CAPTH(J,I),CAPFR(J,I),XKTH(J,I),XKFR(J,I),YYS(i,j+1)
        enddo
	k=1
	NUMSL(I)=num_vl(veg_code(i))+num_gl(geo_code(i))
 	do j=num_vl(veg_code(i))+1,NUMSL(I)
	   VARNOD(J,I)=B1(k,geo_code(i));
	   ACLV(J,I)=B2(k,geo_code(i));
	   BCLV(J,I)=B3(k,geo_code(i));
	   CAPTH(J,I)=B4(k,geo_code(i));
	   CAPFR(J,I)=B5(k,geo_code(i));
	   XKTH(J,I)=B6(k,geo_code(i));
	   XKFR(J,I)=B7(k,geo_code(i));
	   EE(J,I)=0
           layer_thick=layer_thick+B8(geo_code(i),k);
	   YYS(i,j+1)=layer_thick!B8(geo_code(i),j)
	   k=k+1
        enddo
	   YYS(i,NUMSL(I)+1)=YM(NY)
      enddo

	csnow=840000.0
	lbound=2 !1 Dirichlet, 2 heat flux condition at the bottom boundary

	allocate(Y(NY),STAT=IERR)
	allocate(HY(NY),STAT=IERR)
	glm=ym(ny) 
	Y=YM/GLM
	do i=2,ny
	 hy(i)=y(i)-y(i-1)
	enddo
	HCSCALE=glm*glm/perl
	CAPFR=CAPFR*HCSCALE
	CAPTH=CAPTH*HCSCALE
  	CSNOW=CSNOW*HCSCALE
	QPHASE=HCSCALE*333.2*1.D+6

	allocate(U(NX,NY),STAT=IERR)
	allocate(NLX(NX,NY),STAT=IERR)
	allocate(IMEAN(NX),STAT=IERR)
	allocate(YFRNT(NMEAN,NFRMAX,NX),STAT=IERR)
      	allocate(NFRNT(NMEAN,NX),STAT=IERR)
	allocate(TFPV(NNPG,nx),STAT=IERR)	
	
	call  NSLOJS(NNPG,NUMSL,NX,NY,YM,YYS,NLX)
	

end subroutine initialize

subroutine init_cond(q,curr,first,last)
	USE BND
      	USE THERMO
      	USE GRDG
implicit none
   integer q,curr,first,last
   integer i,j
   character*64 inif
   
	if(q.EQ.1)then !ist=1 means reading initial data from 
          do I=first,last
    	    call LININTRP(XTI,YTI(:,I),NXTI,YM,U(I,:),NY)
	  enddo
	elseif(IST.EQ.0)then  			!ist=0 enbales spinup
     	  write(INIF,'(A11,I3.3,A4)') 'dump/start_',curr+1,'.txt'
!	  write(INIF,'(A9,I3.3,A4)') 'in/start_',my_rank+1,'.txt'
      	  open(60,file=INIF,action='READ')
	  read(60,*)TINI
	  do J=1,NY
      		read (60,* ) ( U(i,j),i=first,last)
	  enddo
	  close(60)
	endif

end subroutine init_cond

subroutine active_layer(k)
	USE BND
      	USE THERMO
      	USE GRDG
	USE AL
implicit none

integer :: k,j,jj
real*8 GA,GB,YFRON,GX,GY
real*8 SUNWATER

	  YFRNT(IMEAN(k),:,k)=SLEV
	  NFRNT(IMEAN(k),k)=0
	  do 1329 JJ=1,NY-1
             J=NY-JJ
	     if (YM(J).GE.SLEV.AND.YM(J+1).LE.YFMAX)then
               GA=SUNWATER(U(k,J),NLX(k,J),k)
               GB=SUNWATER(U(k,J+1),NLX(k,J+1),k)
        	  if((GA-WFAZ)*(GB-WFAZ).LE.0.D0) then
        	     GY=(GA-GB)/(YM(J)-YM(J+1))
        	     GX=(GA+GB-GY*(YM(J)+YM(J+1)))/2.D0
        	     if(GY.EQ.0.D0) then
                       YFRON=(YM(J)+YM(J+1))/2.D0
        	     else
                       YFRON=(WFAZ-GX)/GY
        	     endif
        	  else
        	     GOTO 1329
        	  endif
        	  if(NFRNT(IMEAN(k),k).LT.NFRMAX)then
        	    NFRNT(IMEAN(k),k)=NFRNT(IMEAN(k),k)+1
        	    YFRNT(IMEAN(k),NFRNT(IMEAN(k),k),k)=YFRON
        	  endif
	     endif
1329      CONTINUE

end subroutine active_layer

subroutine save_results(k,time1,time2)
USE THERMO
USE GRDG
USE AL
implicit none
	integer :: k,j
     	real*8  :: time1,time2
      	real*8  :: BOUND,SLEVEL

	  RES(IMEAN(k),1)=time1
          RES(IMEAN(k),2)=bound(time2,k)
          RES(IMEAN(k),3)=slevel(k,time2)
	  do  J=1,NMY
            RES(IMEAN(k),J+3)=U(k,MY(J))
	  enddo

end subroutine save_results

!________________________________________________
!__________________FUNCTIONS_____________________
!________________________________________________
      REAL*8 FUNCTION UNWATER(V,NNN,I)
	USE THERMO
        IMPLICIT REAL*8(A-H,O-Z)
     
	TFP=TFPV(NNN,I) ! change I to k0 everywhere except TFP
	E=EE(NNN,I)
    	X=VARNOD(NNN,I)
      	ACL=ACLV(NNN,I)
      	BCL=BCLV(NNN,I)
	
        IF(V.LE.TFP-E)THEN
	       UNWATER=ACL*((DABS(V))**BCL)
	ELSEIF(V.GT.TFP)THEN
	       UNWATER=X
      	ELSE
	       UNWATER=ACL*((DABS(TFP-E))**BCL)
	       UNWATER=UNWATER+(X-UNWATER)*(V+E-TFP)/E
	ENDIF
      RETURN
      END
!-----------------------------------------------
      REAL*8 FUNCTION SUNWATER(V,NNN,I)!Saturated unforzen water
      USE THERMO
      IMPLICIT REAL*8(A-H,O-Z)
      
	TFP=TFPV(NNN,I)
	E=EE(NNN,I)
    	X=VARNOD(NNN,I)
      	ACL=ACLV(NNN,I)
      	BCL=BCLV(NNN,I)
        IF(V.LE.TFP-E)THEN
           SUNWATER=ACL*((DABS(V))**BCL)
        ELSEIF(V.GT.TFP)THEN
           SUNWATER=X
        ELSE
           SUNWATER=ACL*((DABS(TFP-E))**BCL)
           SUNWATER=SUNWATER+(X-SUNWATER)*(V+E-TFP)/E
        ENDIF
	SUNWATER=SUNWATER/X
      RETURN
      END
!-----------------------------------------------
      REAL*8 FUNCTION DUNWATER(V,NNN,I)
      USE THERMO
      IMPLICIT REAL*8(A-H,O-Z)

	TFP=TFPV(NNN,I)
	E=EE(NNN,I)
	X=VARNOD(NNN,I)
	ACL=ACLV(NNN,I)
	BCL=BCLV(NNN,I)

	IF(V.LE.TFP-E)THEN
		DUNWATER=-BCL*ACL*((DABS(V))**(BCL-1.0D0))
	ELSEIF(V.GT.TFP)THEN
      		DUNWATER=0.0D0
	ELSE
		DUNWATER=ACL*((DABS(TFP-E))**BCL)
		DUNWATER=(X-DUNWATER)/E
	ENDIF
      RETURN
      END
!----------------------------------------
      REAL*8 FUNCTION BOUND(T,I)
      USE BND
      REAL*8 T
      INTEGER I,II
	II=1+IDINT((T-TINIR)/STEP)
	BOUND=YTOUT(II,I)+(T+TINI-XTOUT(II)) &
	  *(YTOUT(II+1,I)-YTOUT(II,I))/(XTOUT(II+1)-XTOUT(II))
      RETURN
      END
!----------------------------------------
Subroutine SnowFix(air_temp,sn_depth,n)

real*8, intent (in)  :: air_temp(n)
real*8, intent (out) :: sn_depth(n)
integer :: n

   if(air_temp(1).gt.0.and.sn_depth(1).gt.0)sn_depth(1)=0 
   do i=2,n 
      if(air_temp(i).gt.0.and.sn_depth(i).gt.0)then 
	if (sn_depth(i-1).eq.0)sn_depth(i)=0 ! puts zeros only at the begining of the year
      endif
   enddo

return
end subroutine SnowFix

!----------------------------------------
SUBROUTINE LININTRP(XIN,YIN,NIN,XOUT,YOUT,NOUT)
	! Linear interpolation
      REAL*8 XIN(NIN),YIN(NIN)
      REAL*8 XOUT(NOUT),YOUT(NOUT)
      INTEGER NIN,NOUT
      DO  I=1,NOUT
	 IF(XOUT(I).LE.XIN(1))THEN
	   YOUT(I)=YIN(1)
	   GOTO 1
	 ELSEIF(XOUT(I).GT.XIN(NIN))THEN
           YOUT(I)=YIN(NIN)
	   GOTO 1
         ELSE
           DO J=1,NIN-1
            IF (XIN(J).LT.XOUT(I).AND.XOUT(I).LE.XIN(J+1))THEN
            YOUT(I)=YIN(J)+(XOUT(I)-XIN(J))*(YIN(J+1)-YIN(J))/(XIN(J+1)-XIN(J))
	    GOTO 1
            END IF
           ENDDO
	  ENDIF
1        CONTINUE
      ENDDO
      RETURN
      END
!----------------------------------------
      SUBROUTINE NSLOJS(NNPG,NUMSL,NX,NY,YM,YYS,NLX)
	!assigns correspond layer to the point with with depth
	!starting from surface to the bottom
        INTEGER NX,NY,NNPG,NLX,NUMSL
	DIMENSION NLX(NX,NY),NUMSL(NX)
	REAL*8 YM(ny)
	real YYS(NX,NNPG+1)

      do J=1,NX
	do 6 I=1,NY
	NLX(J,I)=NUMSL(J)
      	  do K=1,NUMSL(J)-1
             IF ( YYS(J,K).LE.YM(I).AND.YM(I).LT.YYS(J,K+1))THEN
	     NLX(J,I)=K
	     GOTO 6
	  ENDIF
	  enddo
6       CONTINUE
      enddo
      RETURN
      END
!----------------------------------------
      REAL*8 FUNCTION SLEVEL(I,t)
	USE BND
	REAL*8 T
	INTEGER I,II
      II=1+IDINT((T-TINIR)/STEP)
      SLEVEL=YSNOUT(II,I)+(T+TINI-XTOUT(II))* &
	(YSNOUT(II+1,I)-YSNOUT(II,I))/(XTOUT(II+1)-XTOUT(II))
	RETURN
      END
!-----------------------------------------------

      REAL*8 FUNCTION GLKY(V,I,J,TLET)
      use bnd
      USE GRDG
      USE THERMO
      IMPLICIT REAL*8(A-H,O-Z)
      integer :: II
      
       XSNOW=SLEV
       dsnow=SLEV-slevel(i,tlet)
       NS=NLX(I,J)
       IF(YM(j).le.dsnow)THEN  		!atmosphere
	      GLKY=1.d4
       ELSEIF (YM(j).Lt.XSNOW)THEN	!snow
              II=1+IDINT((tlet-TINIR)/STEP)
	      glky=RSNOUT(II,I)+(tlet+TINI-XTOUT(II))* &
 			(RSNOUT(II+1,I)-RSNOUT(II,I))/(XTOUT(II+1)-XTOUT(II))
       ELSE				!ground
              WC=UNWATER(V,NS,I)/VARNOD(NS,I)
	      GLKY=(XKTH(NS,I)**WC)*(XKFR(NS,I)**(1.0-WC))
       ENDIF
   	RETURN
      END
!----------------------------------------
	REAL*8 FUNCTION GLCQ(V,NNUS,I)
	USE THERMO
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION NNUS(2),V(2),CL(2)
        
	H=1/(V(1)-V(2)) 
        IF(DABS(V(1)-V(2)).LT.1.D-6) THEN
           GLCQ=0.5d0*(DUNWATER(V(1),NNUS(1),I)+DUNWATER(V(2),NNUS(2),I))
	else
	  if (nnus(1).ne.nnus(2))THEN
      	   GLCQ=0.5D0*( H*(UNWATER(V(1),NNUS(1),I)-UNWATER(V(2),NNUS(1),I))+ &
	   	        H*(UNWATER(V(1),NNUS(2),I)-UNWATER(V(2),NNUS(2),I))    )
	  ELSE
	   GLCQ=H*(UNWATER(V(1),NNUS(1),I)-UNWATER(V(2),NNUS(2),I))
          ENDIF
	endif
        GLCQ=QPHASE*DABS(GLCQ)
      RETURN
      END FUNCTION GLCQ
!----------------------------------------
!----------------------------------------
      REAL*8 FUNCTION GLC(V,I,J)       ! Apparent heat capacity
      USE THERMO
      USE GRDG
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(NY),WW(2),NN(2)
	NS=NLX(I,J)
        XSNOW=SLEV
	if(YM(J).lE.XSNOW)then
 	  GLC=CSNOW            ! heat capacity for snow
	else		    
	   WC=UNWATER(V(J),NS,I)/VARNOD(NS,I)
	   GLC=CAPTH(NS,I)*WC+CAPFR(NS,I)*(1.0-WC)
!	   GLC=GLCUW(V(J),NS,I)	   
!             write(2,*)GLC
	   if(J.GT.(1).AND.J.LT.NY)then
	       WW(1)=(V(J-1)+V(J))/2.D0
	       NN(1)=NLX(I,J-1)
	       WW(2)=V(J)
	       NN(2)=NLX(I,J)
	       GLC=GLC+GLCQ(WW,NN,I)*HY(J)/(HY(J+1)+HY(J))
	       WW(1)=V(J)
	       NN(1)=NLX(I,J)
	       WW(2)=(V(J+1)+V(J))/2.D0
    	       NN(2)=NLX(I,J+1)
	       GLC=GLC+GLCQ(WW,NN,I)*HY(J+1)/(HY(J+1)+HY(J))
	    elseif(J.EQ.1)then
	       WW(1)=V(J)
	       NN(1)=NLX(I,J)
	       WW(2)=(V(J+1)+V(J))/2.D0
    	       NN(2)=NLX(I,J+1)
	       GLC=GLC+GLCQ(WW,NN,I)
	    elseif(J.EQ.NY)then
	       WW(1)=(V(J-1)+V(J))/2.D0
	       NN(1)=NLX(I,J-1)
	       WW(2)=V(J)
    	       NN(2)=NLX(I,J)
	       GLC=GLC+GLCQ(WW,NN,I)
	    endif
	 endif
      RETURN
      END
!-------------------------------------------------------
 SUBROUTINE stefan1D(UU,NY,HY,TTT,TAUM,TMIN,STEP,I,NMS, &
 ITMAX,E1,UWK,NBOUND,flux,BOUND,GLC,GLKY,SUNWATER)
	USE THERMO
      IMPLICIT NONE
      INTEGER NMS,NY,ITMAX,NBOUND,I
      DIMENSION NMS(NY)
      real*8 HY(NY),UU(NY)
      REAL*8 TTT,TAUM,TINI,TMIN,STEP
      REAL*8 E1,UWK
      REAL*8 BOUND,flux,GLC,GLKY,SUNWATER

      INTEGER J,IT
      REAL*8 RAB1,RAB2,AKAPA2,AMU2,Q2
      REAL*8 A,B,C,D
      REAL*8 ALF(NY),BET(NY)
      real*8 U1(NY),U2(NY)
      REAL*8 T1,T,TLET,TAU
      REAL*8 EEY,EEY1,abs1,abs2
      REAL TAU1
      integer :: aaa

      T1=TTT
      TAU1=-1.0
      TAU=TAUM
      UU=U(i,:)
64    continue
      T=T1+TAU
      TLET=T
      U1=UU
      IT=1
      ALF(2)=0.D0
      BET(2)=BOUND(TLET,I)
22    continue
      IF(IT.GT.ITMAX) THEN
	TAU=TAU/2.D0
	TAU1=-1.0
	GOTO 64
      ENDIF
      DO J=2,NY-1
        D=GLC(U1,I,J)/TAU
        A=2.D0*GLKY(U1(J),I,J,TLET)/(HY(J)*(HY(J)+HY(J+1)))
        B=2.D0*GLKY(U1(J+1),I,J+1,TLET)/(HY(J+1)*(HY(J)+HY(J+1)))
        C=A+B+D
        ALF(J+1)=B/(C-A*ALF(J))
        BET(J+1)=(A*BET(J)+D*UU(J))/(C-A*ALF(J))
      ENDDO
      
      RAB1=GLKY(U1(NY),I,NY,TLET)
      RAB2=GLC(U1,I,NY)
      AKAPA2=2.D0*RAB1/(((RAB2*HY(NY)*HY(NY))/TAU+2.D0*RAB1))
      Q2=RAB1*flux
      AMU2=(UU(NY)*RAB2/TAU+2.D0*Q2/HY(NY))/(RAB2/TAU+2.D0*RAB1 &
	 /HY(NY)**2.D0)
      IF(DABS(AKAPA2)>1.D0) then 
        PRINT*,'YOU CAN NOT APPLY PROGONKA ON OY - CHANGE STEPS'
        print*,rab1,rab2,akapa2
        STOP
      endif
      IF (NBOUND.EQ.2)THEN
       U2(NY)=(AMU2+AKAPA2*BET(NY))/(1.D0-ALF(NY)*AKAPA2)
      ELSE
       U2(NY)=flux
      ENDIF
      do J=1,NY-1
        U2(NY-J)=ALF(NY-J+1)*U2(NY-J+1)+BET(NY-J+1)
      enddo

      IF(TAU>TMIN) then !GOTO 11
	do J=1,NY
	  EEY=SUNWATER(U2(J),NMS(J),I)
	  EEY1=SUNWATER(U1(J),NMS(J),I)
	  abs1=DABS(EEY-EEY1)
	  abs2=DABS(U1(J)-U2(J))
	  IF((abs1.GT.UWK).or.(abs2.GT.E1)) then 
	    U1=U2
            IT=IT+1
            GOTO 22
	  endif 
	enddo
      endif
      IF(T.LT.TTT+STEP-1.D-12)THEN
	  T1=T
	  UU=U2
	  IF(TAU1>0) then 
            if(TAU.LT.TAUM) then
        	  TAU=TAU*2.D0
        	  TAU1=-1.0
            end if
          else 
            TAU1=1.0
          endif
	  GOTO 64
      ELSEIF(T.GT.TTT+STEP+1.D-12)THEN
          TAU=(TTT+STEP-T1)
          goto 64
      ELSE
          UU=U2
      ENDIF
      
 END SUBROUTINE stefan1D


subroutine filexist(filename)
    character*64 filename
    logical chf
    inquire(file=filename,exist=chf)
    if (.not.chf) then 
	    write(*,'(/'' FILE '',a, '' DOESNT EXIST'')')trim(filename)
	    stop
    endif
end subroutine filexist!-----------------------------------------------

