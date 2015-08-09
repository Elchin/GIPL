! Geophysical Instatitue Permafrost Laboratory model version 2 GIPLv2
! version 2 is a numerical transient model that employs phase changes and the effect of the unfrozen volumetric  water content in the non-homogeniuos soil texture 
! Original version of the model developed by Romanovsky and Tipenko 2004 and described in Marchenko et al., (2008)
! Current version been significanlty modefied from its predicessor and using the IRF coding design
! This version is maintained by E. Jafarov at INSTAAR, CU Boulder
! Please site Jafarov et al., (2012) work when using it.
module bnd
    integer n_temp                                          ! number of upper boundary points for temperature (input)
    real*8,allocatable::  utemp_time(:), utemp(:,:)         ! upper boundary time and temperature (input)
    real*8,allocatable::  utemp_time_i(:), utemp_i(:,:)     ! time and upper boundary temprature (interpolated)
    integer n_snow                                          ! number of upper boundary points for snow (input)
    real*8 ,allocatable:: snd_time(:),snd(:,:)              ! upper boundary snow time and snow depth (input)
    integer n_stcon
    real*8 ,allocatable:: stcon_time(:),stcon(:,:)          ! snow thermal conductivity time and itself (input)
    real*8 ,allocatable:: snd_i (:,:), stcon_i (:,:)        ! snow depth and thermal conductivity (interpolated)
    real*8 TINIR,time_restart
    integer lbound                                          ! 1 const temp, 2 heat flux condition at the bottom boundary

! Parameter read from cmd file
    integer :: restart                                      ! 0/1 start from previous time step / start from the begining
    real*8 	:: STEP                                         ! step is the timestep in the example it is 1 yr
    real*8 	:: TAUM                                         ! taum is the convergence parameter used by the stefan subroutine
    real*8 	:: TMIN                                         ! tmin minimal timestep used in the Stefan subroutine
    real*8  :: time_beg,time_end                            ! inbegin time, end time
    integer :: itmax                                        ! maximum number of iterations in Stefan subroutine
    integer :: n_time                                       ! number of time steps that temp will be averaged over
    integer :: n_frz_max                                    ! maximum number of freezing fronts
    real*8  :: smooth_coef                                  ! smoothing factor
    real*8  :: unf_water_coef                               ! unfrozen water coefficient
    real*8  :: n_sec_day                                    ! number of second in a day
    real*8  :: frz_frn_max,frz_frn_min                      ! freezing front min and max depth [meters]
    real*8  :: sat_coef                                     ! saturation coefficient [dimensionless, fraction of 1]

end module bnd

module thermo
    real*8 L_fus                                            ! Latent heat of fusion [W/mK]
    real*8 sea_level                                        ! how many meter above the sea level the borehole is

! thermo physical parameters of soil for each soil layer
    real*8,allocatable:: vwc(:,:)                           ! volumetric water content
    real*8,allocatable:: a_coef(:,:),b_coef(:,:)            ! a and b unfrozen water curve coefficients
    real*8,allocatable:: temp_frz(:,:)                      ! temperature freezing depression
    real*8,allocatable:: EE(:,:)
    real*8,allocatable:: hcap_frz(:,:),hcap_thw(:,:)        ! soil layer heat capacity thawed/frozen
    real*8,allocatable:: tcon_frz(:,:),tcon_thw(:,:)        ! soil layer thermal conductivity thawed/frozen

    real*8 shcap                                            ! heat capacity of snow (constant)

    real*8, allocatable :: temp(:,:)                        ! soil temperature
    real, allocatable:: n_bnd_lay(:,:)                      ! number of boundaries between layer in soil
    integer k0
	

    integer, allocatable :: snow_code(:),veg_code(:)        ! (not necccessary) required for runing in parallel
    integer, allocatable :: geo_code(:),gt_zone_code(:)     ! (not necccessary) required for runing in parallel
    real*8, allocatable  :: temp_grd(:)                     ! temprature gradient at the lower boundary

    real*8 ,allocatable:: RES(:,:)                          ! unified variable for the writing results into the file

end module thermo

module grd
    integer, parameter :: n_lay=10                          ! total allowed number of soil layer
    integer :: n_site                                       ! number of sites
    integer :: n_grd                                        ! total number of grid points with depth (grid.txt)
    real*8,allocatable:: zdepth(:),dz(:)                    ! vertical grid and distance between grid point 'zdepth(n_grd)'
    integer,allocatable:: lay_id(:,:)                       ! layer index
    integer :: m_grd                                        ! number of grid points to store in res file
    integer,allocatable:: zdepth_id(:)                      ! index vector of stored grid points 'zdepth_id(m_grid)'
    integer,allocatable:: NUMSL(:)
    integer :: n_ini                                        ! number of vertical grid cells in init file
    real*8, allocatable :: zdepth_ini(:),ztemp_ini(:,:)     ! depth and correspoding initial temperature (time=0) 'zdepth_ini(n_ini)'


end module grd

module alt
    integer,allocatable::n_frz_frn(:,:)                     ! number of freezing front (e.g. when freezup is about to happened)
    integer,allocatable::i_time(:)                          ! internal time step with the the main loop
    real*8 ,allocatable::z_frz_frn(:,:,:)                   ! depth of the freezing front
end module alt


program permgis
use bnd
use thermo
use grd
use alt

implicit none
! all function names start with letter 'f'
    external futemp,ftcon,fsat_unf_water,fapp_hcap
! functions
    real*8 funf_water                                       ! unfrozen water
    real*8 fdunf_water                                      ! derivative of the unfrozen water
    real*8 fhcap                                            ! heat capacity
    real*8 ftcon                                            ! thermal conductivity
    real*8 fapp_hcap                                        ! apparent heat capacity
    real*8 futemp                                           ! temperature interpolation
    real*8 fsnow_level                                      ! snow level
    real*8 fsat_unf_water                                   ! saturated unfrozen water
! variables
    real*8 ,allocatable:: res_save(:,:)                     ! save results into 2D array
    real*8 ,allocatable:: dfrz_frn(:)                       ! depth of the freezing front
    real frz_up_time_cur                                    ! freezeup time current (within a year)
    real frz_up_time_tot                                    ! freezeup time global
! counters (time,steps)
    real*8 time_s,time_e                                    ! internal start and end times
    real*8, allocatable:: time_loop(:)                      ! main looping time
    real*8, allocatable:: time_cur(:)                       ! current time (e.g. current day)
    integer :: n_itime                                      ! total number of internal time steps
! other counters
    integer :: I,II,J,JJ
! output file names
    character(64) :: restart_file,result_file,aver_res_file
! putput files formats
    character(210) :: FMT1,FMT2                             ! results formating type

! MPI varriables----------------------------------------------
! This variable are for parallel code
    integer   my_rank  ! My process rank.
    integer   mype     ! The number of processes.
    integer   nx1,nx2  ! partitioning with dimensions nx2-nx1
    integer   subset ! n_site/mype
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

    subset = n_site/mype      ! subset length
    nx1 = my_rank*subset+1
    nx2 = (my_rank+1)*subset
    if ((my_rank+1).eq.mype) nx2=(my_rank+1-mype)*subset + n_site

    allocate(res_save(m_grd+3,n_site),STAT=IERR)
    allocate(dfrz_frn(n_time),STAT=IERR)
    allocate(time_loop(n_site),STAT=IERR)
    allocate(time_cur(n_site),STAT=IERR)
    allocate(RES(n_time,m_grd+3),STAT=IERR)
!________________INITIALIZATION________________________

!      n_itime=n_time*time_end+3
    n_itime=n_time+2
    time_s=STEP*DBLE(n_time*time_beg)
    time_e=STEP*DBLE(n_time*time_end)
    time_loop=0.0D0
    TINIR=0.0D0
!-----------------------------------------------initial data read
    call init_cond(restart,my_rank,nx1,nx2)
	
    write(aver_res_file,'(A10,I3.3,A4)') 'dump/mean_',my_rank+1,'.txt'
    write(result_file,'(A12,I3.3,A4)') 'dump/result_',my_rank+1,'.txt'
    write(restart_file,'(A11,I3.3,A4)') 'dump/start_',my_rank+1,'.txt'
    open(1,file=result_file,STATUS='unknown')
    open(2,file=aver_res_file,STATUS='unknown')
    open(3,file=restart_file,STATUS='unknown')
    write(FMT1,'(A30,I0,A12)')'(1x,I10,1x,F12.3,2(1x,F16.12),',m_grd,'(1x,F16.12))'
    write(FMT2,'(A28,I0,A40)')'(1x,I10,1x,F12.3,2(1x,F8.3),',m_grd,'(1x,F8.3),(1x,F8.3,1x,F12.3),(1x,F12.3))'

    allocate(utemp_time_i(n_itime),STAT=IERR)                  ! allocating interval varialbe after interation
    allocate(utemp_i(n_itime,n_site),STAT=IERR)
    allocate(snd_i(n_itime,n_site),STAT=IERR)
    allocate(stcon_i(n_itime,n_site),STAT=IERR)
	
!****************************************************************************	
    do I=1,n_time+2
        utemp_time_i(I)=time_restart+DBLE(I-1)*STEP
    enddo
    i_time=1

    do I=nx1,nx2
        if (lbound.EQ.2)temp_grd(i)=temp_grd(i)*zdepth(n_grd)
        do ii=1,NUMSL(i)
            temp_frz(ii,i)=-(vwc(ii,i)/a_coef(ii,i))**(1.d0/b_coef(ii,i))
        enddo
        call interpolate(utemp_time,utemp(:,I),n_temp,utemp_time_i,utemp_i(:,I),n_itime)
        call interpolate(snd_time,snd(:,I),n_snow,utemp_time_i,snd_i(:,I),n_itime)
        call snowfix(utemp_i(:,I),snd_i(:,I),n_itime)
        call interpolate(stcon_time,stcon(:,I),n_stcon,utemp_time_i,stcon_i(:,I),n_itime)
        call active_layer(I)
    enddo

! begining of the major loop
65 CONTINUE
    do I=nx1,nx2
        time_cur(I)=time_loop(I)+time_restart
        call save_results(I,time_cur(I),time_loop(I))
        6666  continue
        call stefan1D(temp(I,:),n_grd,dz,time_loop(I),TAUM,TMIN,STEP,I,lay_id(I,:) &
        ,ITMAX,smooth_coef,unf_water_coef,lbound,temp_grd(I),futemp,fapp_hcap,ftcon,fsat_unf_water)
!--------------------------------------------
!--------------------------------------------
!	    do j=1,n_grd			! WRITTING RESULTS
!	       write(1,'3(f12.3)') temp(I,j),zdepth(j),ftcon(temp(I,j),I,J,time_cur(I))
!	    enddo
        time_loop(I)=time_loop(I)+STEP
        time_cur(I)=time_loop(I)+time_restart
        if(i_time(I).LT.n_time)  then
            i_time(I)=i_time(I)+1
            call save_results(I,time_cur(I),time_loop(I))
            call active_layer(I)
            GOTO 6666
        endif
        if(time_s.LT.time_e.AND.time_loop(nx1).GT.time_s)then
            do II=1,n_time			! WRITTING RESULTS
                write(1,FMT1) I, (RES(II,JJ),JJ=1,m_grd+3)
            enddo
        endif
      	do jj=1,m_grd+3
            res_save(jj,i)=sum((RES(:,JJ)))
        enddo
     enddo

      i_time=1
      do I=nx1,nx2
           frz_up_time_cur=-7777.D0
	   frz_up_time_tot=frz_up_time_cur
           do J=2,n_time
              if((n_frz_frn(J,I)-n_frz_frn(J-1,I)).EQ.-2)then
        	if(z_frz_frn(J-1,n_frz_frn(J-1,I),I).GE.frz_frn_min) frz_up_time_cur=SNGL(RES(J,1))
              endif
           enddo
	   if(frz_up_time_cur.GT.0.0)then
	       frz_up_time_tot=AMOD(frz_up_time_cur,REAL(n_time))
	       if(frz_up_time_tot.EQ.0.0)frz_up_time_tot=REAL(n_time)
           endif
           dfrz_frn=z_frz_frn(:,1,I)

  	   call save_results(I,time_cur(I),time_loop(I))

	   call active_layer(I)

	   !____WRITTING MEAN
	   write(2,FMT2) I,(res_save(JJ,I)/DBLE(n_time),JJ=1,m_grd+3), &
 		   dfrz_frn(n_time),frz_up_time_cur,frz_up_time_tot	
	   do j=1,n_itime
	     utemp_time_i(j)=time_cur(nx1)+DBLE(j-1)*STEP
	   enddo
      	   call interpolate(utemp_time,utemp(:,I),n_temp,utemp_time_i,utemp_i(:,I),n_itime)
      	   call interpolate(snd_time,snd(:,I),n_snow,utemp_time_i,snd_i(:,I),n_itime)
	   call snowfix(utemp_i(:,I),snd_i(:,I),n_itime)
	   call interpolate(stcon_time,stcon(:,I),n_stcon,utemp_time_i,stcon_i(:,I),n_itime)
      enddo

      rewind(3) ! -------------start file writting begin     
      write(3, * ) time_restart
!      write(3, * ) time_cur(1)
      do J=1,n_grd
         write (3,* ) ( temp(II,J),II=nx1,nx2)
      enddo     ! -------------start file writting end     

      TINIR=time_loop(nx1)
!************************************************************
6444    if(time_loop(nx1).LT.time_e)GO TO 65
	close(1);close(2);close(3)
!   call MPI_FINALIZE(ierr)

end ! end of main program

subroutine initialize

use bnd
use thermo
use grd
use alt

implicit none

    integer IREAD,ierr
    integer :: i,j,k,z_num

    real*8, allocatable:: Y(:)
    real*8 glm,hcscale
	
    real*8 ,allocatable ::gtzone(:,:)
    character*64 stdummy
    character*64 fconfig

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

    fconfig='gipl_config.cfg'
    call filexist(fconfig)
    open(60,file=fconfig)
        read(60,'(A)')finput
        read(60,'(A)')boundf
        read(60,'(A)')snowf
        read(60,'(A)')rsnowf
        read(60,'(A)')inif
        read(60,'(A)')cmdf
        read(60,'(A)')cetkaf
        read(60,'(A)')fvegetation
        read(60,'(A)')fgeology
    close(60)

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
   read(60,*)n_site
        allocate(snow_code(n_site),STAT=IERR)
        allocate(veg_code(n_site),STAT=IERR)
        allocate(geo_code(n_site),STAT=IERR)
        allocate(gt_zone_code(n_site),STAT=IERR)
        allocate(temp_grd(n_site),STAT=IERR)
        do i=1,n_site
          read(60,*) IREAD,snow_code(i),veg_code(i),geo_code(i),&
                    gt_zone_code(i),temp_grd(i)
        enddo
    close(60)
!   print*, trim(finput),' has been read'
      
    open(60,file=BOUNDF)
    read(60,*)n_temp
    allocate(utemp_time(n_temp),STAT=IERR)
    allocate(utemp(n_temp,n_site),STAT=IERR)
    do i=1,n_temp
        read(60,*) utemp_time(I),(utemp(I,J),J=1,n_site)
    enddo
    close(60)
!   print*,trim(boundf),' has been read'

    open(60,file=RSNOWF)
    read(60,*)n_stcon
    allocate(stcon_time(n_stcon),STAT=IERR)
    allocate(stcon(n_stcon,n_site),STAT=IERR)
    do i=1,n_stcon
        read(60,*) stcon_time(i),(stcon(i,J),J=1,n_site)
    enddo
    close(60)
!   print*,trim(rsnowf),' has been read'

    open(60,file=SNOWF)
    read(60,*)n_snow
    allocate(snd_time(n_snow),STAT=IERR)
    allocate(snd(n_snow,n_site),STAT=IERR)
    do I=1,n_snow
       read(60,*) snd_time(i),(snd(i,J),J=1,n_site)
    enddo
      close(60)
!   print*,trim(snowf),' has been read'

    open(60,file=inif,action='read')
    read(60,*)z_num,n_ini!,time_restart
    allocate(zdepth_ini(n_ini),STAT=IERR)
    allocate(ztemp_ini(n_ini,n_site),STAT=IERR)
    allocate(gtzone(n_ini,z_num+1),STAT=IERR)
    read(60,*)stdummy
        do i=1,n_ini
            read(60,*) (gtzone(i,j),j=1,z_num+1)
        enddo
    close(60)
!   print*,trim(inif),'has been read'

    time_restart=utemp_time(1)
    zdepth_ini(:)=gtzone(:,1)
    do i=1,n_site
        k=gt_zone_code(i)
        ztemp_ini(:,I)=gtzone(:,k+1)
    enddo

    open(60,FILE=cmdf)
        read(60,*)restart
        read(60,*)STEP,TAUM,TMIN
        read(60,*) time_beg,time_end
        read(60,*) smooth_coef,unf_water_coef,itmax  !smoothing factor | unfrozen water parameter | max number of iterations
        read(60,*) n_sec_day,n_time ! number of second in a day [sec] | number of time steps (in the example number of days in a year )
        read(60,*) sea_level,n_frz_max
        read(60,*) frz_frn_min,frz_frn_max
        read(60,*) sat_coef
    close(60)
!   print*,trim(cmdf),' has been read'


    open(60,file=cetkaf)
    read(60,*)n_grd
    allocate(zdepth(n_grd),STAT=IERR)
        do i=1,n_grd
            read(60,*) zdepth(i)
        enddo
        read(60,*)m_grd
        allocate(zdepth_id(m_grd),STAT=IERR)
        do j=1,m_grd
            read(60,*)zdepth_id(j)
        enddo
    close(60)
!   print*,trim(cetkaf),' has been read'

! note: that all max NUMSL layers has to be read or it will a give segmantation error
!      n_lay=10!MAXVAL(NUMSL)
!----------------------------------------------------      
    open (60, file=fvegetation)
    read(60,*) vln ! reads numbers of  classes
    allocate(A1(n_lay,vln),STAT=IERR) ! vwc
    allocate(A2(n_lay,vln),STAT=IERR) ! a_coef
    allocate(A3(n_lay,vln),STAT=IERR) ! b_coef
    allocate(A4(n_lay,vln),STAT=IERR) ! hcap_frz
    allocate(A5(n_lay,vln),STAT=IERR)  !hcap_thw
    allocate(A6(n_lay,vln),STAT=IERR)  !tcon_frz
    allocate(A7(n_lay,vln),STAT=IERR)  !tcon_thw
    allocate(A8(vln,n_lay),STAT=IERR)  !bot_cond
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
!   print*,trim(fvegetation),' has been read'

    open (60, file='in/geo.txt')
    read(60,*) gln ! reads numbers of  classes
    allocate(B1(n_lay,gln),STAT=IERR) ! vwc
    allocate(B2(n_lay,gln),STAT=IERR) ! a_coef
    allocate(B3(n_lay,gln),STAT=IERR) ! b_coef
    allocate(B4(n_lay,gln),STAT=IERR) ! hcap_frz
    allocate(B5(n_lay,gln),STAT=IERR)  !hcap_thw
    allocate(B6(n_lay,gln),STAT=IERR)  !tcon_frz
    allocate(B7(n_lay,gln),STAT=IERR)  !tcon_thw
    allocate(B8(gln,n_lay),STAT=IERR)  !bot_cond
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


      allocate(vwc(n_lay,n_site),STAT=IERR)
      allocate(a_coef(n_lay,n_site),STAT=IERR)
      allocate(b_coef(n_lay,n_site),STAT=IERR)
      allocate(EE(n_lay,n_site),STAT=IERR)
      allocate(hcap_frz(n_lay,n_site),STAT=IERR)
      allocate(hcap_thw(n_lay,n_site),STAT=IERR)
      allocate(tcon_frz(n_lay,n_site),STAT=IERR)
      allocate(tcon_thw(n_lay,n_site),STAT=IERR)
      allocate(numsl(n_site),STAT=IERR)
      allocate(n_bnd_lay(n_site,n_lay+1),STAT=IERR)

      do i = 1,n_site
	layer_thick=0
	n_bnd_lay(i,1)=layer_thick
	  	layer_thick=0
	n_bnd_lay(i,1)=layer_thick
 	do j=1,num_vl(veg_code(i))
	   vwc(J,I)=A1(j,veg_code(i));
	   a_coef(J,I)=A2(j,veg_code(i));
	   b_coef(J,I)=A3(j,veg_code(i));
	   hcap_thw(J,I)=A4(j,veg_code(i));
	   hcap_frz(J,I)=A5(j,veg_code(i));
	   tcon_thw(J,I)=A6(j,veg_code(i));
	   tcon_frz(J,I)=A7(j,veg_code(i));
	   if (j.eq.1) then 
	      layer_thick=A8(veg_code(i),j)
	   else
	      layer_thick=layer_thick+A8(veg_code(i),j);
	   endif
	   n_bnd_lay(i,j+1)=layer_thick
	   EE(J,I)=0
!	     write(*,'(3(f8.3),2(f12.1),3(f8.3))') vwc(J,I),a_coef(J,I),b_coef(J,I), &
!		  hcap_thw(J,I),hcap_frz(J,I),tcon_thw(J,I),tcon_frz(J,I),n_bnd_lay(i,j+1)
        enddo
	k=1
	NUMSL(I)=num_vl(veg_code(i))+num_gl(geo_code(i))
 	do j=num_vl(veg_code(i))+1,NUMSL(I)
	   vwc(J,I)=B1(k,geo_code(i));
	   a_coef(J,I)=B2(k,geo_code(i));
	   b_coef(J,I)=B3(k,geo_code(i));
	   hcap_thw(J,I)=B4(k,geo_code(i));
	   hcap_frz(J,I)=B5(k,geo_code(i));
	   tcon_thw(J,I)=B6(k,geo_code(i));
	   tcon_frz(J,I)=B7(k,geo_code(i));
	   EE(J,I)=0
           layer_thick=layer_thick+B8(geo_code(i),k);
	   n_bnd_lay(i,j+1)=layer_thick!B8(geo_code(i),j)
	   k=k+1
        enddo
	   n_bnd_lay(i,NUMSL(I)+1)=zdepth(n_grd)
      enddo

	shcap=840000.0
	lbound=2 !1 Dirichlet, 2 heat flux condition at the bottom boundary

	allocate(Y(n_grd),STAT=IERR)
	allocate(dz(n_grd),STAT=IERR)
	glm=zdepth(n_grd) 
	Y=zdepth/GLM
	do i=2,n_grd
	 dz(i)=y(i)-y(i-1)
	enddo
	HCSCALE=glm*glm/n_sec_day
	hcap_frz=hcap_frz*HCSCALE
	hcap_thw=hcap_thw*HCSCALE
  	shcap=shcap*HCSCALE
	L_fus=HCSCALE*333.2*1.D+6

	allocate(temp(n_site,n_grd),STAT=IERR)
	allocate(lay_id(n_site,n_grd),STAT=IERR)
	allocate(i_time(n_site),STAT=IERR)
	allocate(z_frz_frn(n_time,n_frz_max,n_site),STAT=IERR)
      	allocate(n_frz_frn(n_time,n_site),STAT=IERR)
	allocate(temp_frz(n_lay,n_site),STAT=IERR)	
	
	call  NSLOJS(n_lay,NUMSL,n_site,n_grd,zdepth,n_bnd_lay,lay_id)
	

end subroutine initialize

subroutine init_cond(q,curr,first,last)

use bnd
use thermo
use grd

implicit none
   integer q,curr,first,last
   integer i,j
   character*64 inif

    if(q.EQ.1)then !restart=1 means reading initial data from
        do I=first,last
            call interpolate(zdepth_ini,ztemp_ini(:,I),n_ini,zdepth,temp(I,:),n_grd)
        enddo
    elseif(restart.EQ.0)then  			!restart=0 enbales spinup
        write(INIF,'(A11,I3.3,A4)') 'dump/start_',curr+1,'.txt'
!       write(INIF,'(A9,I3.3,A4)') 'in/start_',my_rank+1,'.txt'
        open(60,file=INIF,action='READ')
            read(60,*)time_restart              ! day number in restart file
            do J=1,n_grd
                read (60,* ) ( temp(i,j),i=first,last)
            enddo
        close(60)
    endif

end subroutine init_cond

subroutine active_layer(k)

use bnd
use thermo
use grd
use alt

implicit none

    integer :: k,j,jj
    real*8 GA,GB,YFRON,GX,GY
    real*8 fsat_unf_water

    z_frz_frn(i_time(k),:,k)=sea_level
    n_frz_frn(i_time(k),k)=0
    do 1329 JJ=1,n_grd-1
        J=n_grd-JJ
        if (zdepth(J).GE.sea_level.AND.zdepth(J+1).LE.frz_frn_max)then
            GA=fsat_unf_water(temp(k,J),lay_id(k,J),k)
            GB=fsat_unf_water(temp(k,J+1),lay_id(k,J+1),k)
            if((GA-sat_coef)*(GB-sat_coef).LE.0.D0) then
                GY=(GA-GB)/(zdepth(J)-zdepth(J+1))
                GX=(GA+GB-GY*(zdepth(J)+zdepth(J+1)))/2.D0
                if(GY.EQ.0.D0) then
                    YFRON=(zdepth(J)+zdepth(J+1))/2.D0
                else
                    YFRON=(sat_coef-GX)/GY
                endif
            else
            GOTO 1329
        endif
        if(n_frz_frn(i_time(k),k).LT.n_frz_max)then
            n_frz_frn(i_time(k),k)=n_frz_frn(i_time(k),k)+1
            z_frz_frn(i_time(k),n_frz_frn(i_time(k),k),k)=YFRON
            endif
        endif
    1329 CONTINUE

end subroutine active_layer

subroutine save_results(k,time1,time2)
USE thermo
USE grd
USE alt
implicit none
	integer :: k,j
     	real*8  :: time1,time2
      	real*8  :: futemp,fsnow_level

	  RES(i_time(k),1)=time1
          RES(i_time(k),2)=futemp(time2,k)
          RES(i_time(k),3)=fsnow_level(k,time2)
	  do  J=1,m_grd
            RES(i_time(k),J+3)=temp(k,zdepth_id(J))
	  enddo

end subroutine save_results

!________________________________________________
!__________________FUNCTIONS_____________________
!________________________________________________
      REAL*8 FUNCTION funf_water(V,NNN,I)
	USE thermo
        IMPLICIT REAL*8(A-H,O-Z)
     
	TFP=temp_frz(NNN,I) ! change I to k0 everywhere except TFP
	E=EE(NNN,I)
    	X=vwc(NNN,I)
      	ACL=a_coef(NNN,I)
      	BCL=b_coef(NNN,I)
	
        IF(V.LE.TFP-E)THEN
	       funf_water=ACL*((DABS(V))**BCL)
	ELSEIF(V.GT.TFP)THEN
	       funf_water=X
      	ELSE
	       funf_water=ACL*((DABS(TFP-E))**BCL)
	       funf_water=funf_water+(X-funf_water)*(V+E-TFP)/E
	ENDIF
      RETURN
      END
!-----------------------------------------------
      REAL*8 FUNCTION fsat_unf_water(V,NNN,I)!Saturated unforzen water
      USE thermo
      IMPLICIT REAL*8(A-H,O-Z)
      
	TFP=temp_frz(NNN,I)
	E=EE(NNN,I)
    	X=vwc(NNN,I)
      	ACL=a_coef(NNN,I)
      	BCL=b_coef(NNN,I)
        IF(V.LE.TFP-E)THEN
           fsat_unf_water=ACL*((DABS(V))**BCL)
        ELSEIF(V.GT.TFP)THEN
           fsat_unf_water=X
        ELSE
           fsat_unf_water=ACL*((DABS(TFP-E))**BCL)
           fsat_unf_water=fsat_unf_water+(X-fsat_unf_water)*(V+E-TFP)/E
        ENDIF
	fsat_unf_water=fsat_unf_water/X
      RETURN
      END
!-----------------------------------------------
      REAL*8 FUNCTION fdunf_water(V,NNN,I)
      USE thermo
      IMPLICIT REAL*8(A-H,O-Z)

	TFP=temp_frz(NNN,I)
	E=EE(NNN,I)
	X=vwc(NNN,I)
	ACL=a_coef(NNN,I)
	BCL=b_coef(NNN,I)

	IF(V.LE.TFP-E)THEN
		fdunf_water=-BCL*ACL*((DABS(V))**(BCL-1.0D0))
	ELSEIF(V.GT.TFP)THEN
      		fdunf_water=0.0D0
	ELSE
		fdunf_water=ACL*((DABS(TFP-E))**BCL)
		fdunf_water=(X-fdunf_water)/E
	ENDIF
      RETURN
      END
!----------------------------------------
      REAL*8 FUNCTION futemp(T,I)
      USE bnd
      REAL*8 T
      INTEGER I,II
	II=1+IDINT((T-TINIR)/STEP)
	futemp=utemp_i(II,I)+(T+time_restart-utemp_time_i(II)) &
	  *(utemp_i(II+1,I)-utemp_i(II,I))/(utemp_time_i(II+1)-utemp_time_i(II))
      RETURN
      END
!----------------------------------------
subroutine snowfix(air_temp,stcon,n)

real*8, intent (in)  :: air_temp(n)
real*8, intent (out) :: stcon(n)
integer :: n

   if(air_temp(1).gt.0.and.stcon(1).gt.0)stcon(1)=0 
   do i=2,n 
      if(air_temp(i).gt.0.and.stcon(i).gt.0)then 
	if (stcon(i-1).eq.0)stcon(i)=0 ! puts zeros only at the begining of the year
      endif
   enddo

return
end subroutine snowfix

!----------------------------------------
SUBROUTINE interpolate(XIN,YIN,NIN,XOUT,YOUT,n_itime)
	! Linear interpolation
      REAL*8 XIN(NIN),YIN(NIN)
      REAL*8 XOUT(n_itime),YOUT(n_itime)
      INTEGER NIN,n_itime
      DO  I=1,n_itime
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
      SUBROUTINE NSLOJS(n_lay,NUMSL,n_site,n_grd,zdepth,n_bnd_lay,lay_id)
	!assigns correspond layer to the point with with depth
	!starting from surface to the bottom
        INTEGER n_site,n_grd,n_lay,lay_id,NUMSL
	DIMENSION lay_id(n_site,n_grd),NUMSL(n_site)
	REAL*8 zdepth(n_grd)
	real n_bnd_lay(n_site,n_lay+1)

      do J=1,n_site
	do 6 I=1,n_grd
	lay_id(J,I)=NUMSL(J)
      	  do K=1,NUMSL(J)-1
             IF ( n_bnd_lay(J,K).LE.zdepth(I).AND.zdepth(I).LT.n_bnd_lay(J,K+1))THEN
	     lay_id(J,I)=K
	     GOTO 6
	  ENDIF
	  enddo
6       CONTINUE
      enddo
      RETURN
      END
!----------------------------------------
      REAL*8 FUNCTION fsnow_level(I,t)
	USE bnd
	REAL*8 T
	INTEGER I,II
      II=1+IDINT((T-TINIR)/STEP)
      fsnow_level=snd_i(II,I)+(T+time_restart-utemp_time_i(II))* &
	(snd_i(II+1,I)-snd_i(II,I))/(utemp_time_i(II+1)-utemp_time_i(II))
	RETURN
      END
!-----------------------------------------------

      REAL*8 FUNCTION ftcon(V,I,J,time_cur)
      use bnd
      USE grd
      USE thermo
      IMPLICIT REAL*8(A-H,O-Z)
      integer :: II
      
       XSNOW=sea_level
       dsnow=sea_level-fsnow_level(i,time_cur)
       NS=lay_id(I,J)
       IF(zdepth(j).le.dsnow)THEN  		!atmosphere
	      ftcon=1.d4
       ELSEIF (zdepth(j).Lt.XSNOW)THEN	!snow
              II=1+IDINT((time_cur-TINIR)/STEP)
	      ftcon=stcon_i(II,I)+(time_cur+time_restart-utemp_time_i(II))* &
 			(stcon_i(II+1,I)-stcon_i(II,I))/(utemp_time_i(II+1)-utemp_time_i(II))
       ELSE				!ground
              WC=funf_water(V,NS,I)/vwc(NS,I)
	      ftcon=(tcon_thw(NS,I)**WC)*(tcon_frz(NS,I)**(1.0-WC))
       ENDIF
   	RETURN
      END
!----------------------------------------
	REAL*8 FUNCTION fhcap(V,NNUS,I)
	USE thermo
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION NNUS(2),V(2),CL(2)
        
	H=1/(V(1)-V(2)) 
        IF(DABS(V(1)-V(2)).LT.1.D-6) THEN
           fhcap=0.5d0*(fdunf_water(V(1),NNUS(1),I)+fdunf_water(V(2),NNUS(2),I))
	else
	  if (nnus(1).ne.nnus(2))THEN
      	   fhcap=0.5D0*( H*(funf_water(V(1),NNUS(1),I)-funf_water(V(2),NNUS(1),I))+ &
	   	        H*(funf_water(V(1),NNUS(2),I)-funf_water(V(2),NNUS(2),I))    )
	  ELSE
	   fhcap=H*(funf_water(V(1),NNUS(1),I)-funf_water(V(2),NNUS(2),I))
          ENDIF
	endif
        fhcap=L_fus*DABS(fhcap)
      RETURN
      END FUNCTION fhcap
!----------------------------------------
!----------------------------------------
      REAL*8 FUNCTION fapp_hcap(V,I,J)       ! Apparent heat capacity
      USE thermo
      USE grd
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(n_grd),WW(2),NN(2)
	NS=lay_id(I,J)
        XSNOW=sea_level
	if(zdepth(J).lE.XSNOW)then
 	  fapp_hcap=shcap            ! heat capacity for snow
	else		    
	   WC=funf_water(V(J),NS,I)/vwc(NS,I)
	   fapp_hcap=hcap_thw(NS,I)*WC+hcap_frz(NS,I)*(1.0-WC)
!	   fapp_hcap=GLCUW(V(J),NS,I)	   
!             write(2,*)fapp_hcap
	   if(J.GT.(1).AND.J.LT.n_grd)then
	       WW(1)=(V(J-1)+V(J))/2.D0
	       NN(1)=lay_id(I,J-1)
	       WW(2)=V(J)
	       NN(2)=lay_id(I,J)
	       fapp_hcap=fapp_hcap+fhcap(WW,NN,I)*dz(J)/(dz(J+1)+dz(J))
	       WW(1)=V(J)
	       NN(1)=lay_id(I,J)
	       WW(2)=(V(J+1)+V(J))/2.D0
    	       NN(2)=lay_id(I,J+1)
	       fapp_hcap=fapp_hcap+fhcap(WW,NN,I)*dz(J+1)/(dz(J+1)+dz(J))
	    elseif(J.EQ.1)then
	       WW(1)=V(J)
	       NN(1)=lay_id(I,J)
	       WW(2)=(V(J+1)+V(J))/2.D0
    	       NN(2)=lay_id(I,J+1)
	       fapp_hcap=fapp_hcap+fhcap(WW,NN,I)
	    elseif(J.EQ.n_grd)then
	       WW(1)=(V(J-1)+V(J))/2.D0
	       NN(1)=lay_id(I,J-1)
	       WW(2)=V(J)
    	       NN(2)=lay_id(I,J)
	       fapp_hcap=fapp_hcap+fhcap(WW,NN,I)
	    endif
	 endif
      RETURN
      END
!-------------------------------------------------------
 SUBROUTINE stefan1D(UU,n_grd,dz,time_loop,TAUM,TMIN,STEP,I,NMS, &
 ITMAX,smooth_coef,unf_water_coef,NBOUND,flux,futemp,fapp_hcap,ftcon,fsat_unf_water)
	USE thermo
      IMPLICIT NONE
      INTEGER NMS,n_grd,ITMAX,NBOUND,I
      DIMENSION NMS(n_grd)
      real*8 dz(n_grd),UU(n_grd)
      REAL*8 time_loop,TAUM,time_restart,TMIN,STEP
      REAL*8 smooth_coef,unf_water_coef
      REAL*8 futemp,flux,fapp_hcap,ftcon,fsat_unf_water

      INTEGER J,IT
      REAL*8 RAB1,RAB2,AKAPA2,AMU2,Q2
      REAL*8 A,B,C,D
      REAL*8 ALF(n_grd),BET(n_grd)
      real*8 U1(n_grd),U2(n_grd)
      REAL*8 T1,T,time_cur,TAU
      REAL*8 EEY,EEY1,abs1,abs2
      REAL TAU1
      integer :: aaa

      T1=time_loop
      TAU1=-1.0
      TAU=TAUM
      UU=temp(i,:)
64    continue
      T=T1+TAU
      time_cur=T
      U1=UU
      IT=1
      ALF(2)=0.D0
      BET(2)=futemp(time_cur,I)
22    continue
      IF(IT.GT.ITMAX) THEN
	TAU=TAU/2.D0
	TAU1=-1.0
	GOTO 64
      ENDIF
      DO J=2,n_grd-1
        D=fapp_hcap(U1,I,J)/TAU
        A=2.D0*ftcon(U1(J),I,J,time_cur)/(dz(J)*(dz(J)+dz(J+1)))
        B=2.D0*ftcon(U1(J+1),I,J+1,time_cur)/(dz(J+1)*(dz(J)+dz(J+1)))
        C=A+B+D
        ALF(J+1)=B/(C-A*ALF(J))
        BET(J+1)=(A*BET(J)+D*UU(J))/(C-A*ALF(J))
      ENDDO
      
      RAB1=ftcon(U1(n_grd),I,n_grd,time_cur)
      RAB2=fapp_hcap(U1,I,n_grd)
      AKAPA2=2.D0*RAB1/(((RAB2*dz(n_grd)*dz(n_grd))/TAU+2.D0*RAB1))
      Q2=RAB1*flux
      AMU2=(UU(n_grd)*RAB2/TAU+2.D0*Q2/dz(n_grd))/(RAB2/TAU+2.D0*RAB1 &
	 /dz(n_grd)**2.D0)
      IF(DABS(AKAPA2)>1.D0) then 
        PRINT*,'YOU CAN NOT APPLY PROGONKA ON OY - CHANGE STEPS'
        print*,rab1,rab2,akapa2
        STOP
      endif
      IF (NBOUND.EQ.2)THEN
       U2(n_grd)=(AMU2+AKAPA2*BET(n_grd))/(1.D0-ALF(n_grd)*AKAPA2)
      ELSE
       U2(n_grd)=flux
      ENDIF
      do J=1,n_grd-1
        U2(n_grd-J)=ALF(n_grd-J+1)*U2(n_grd-J+1)+BET(n_grd-J+1)
      enddo

      IF(TAU>TMIN) then !GOTO 11
	do J=1,n_grd
	  EEY=fsat_unf_water(U2(J),NMS(J),I)
	  EEY1=fsat_unf_water(U1(J),NMS(J),I)
	  abs1=DABS(EEY-EEY1)
	  abs2=DABS(U1(J)-U2(J))
	  IF((abs1.GT.unf_water_coef).or.(abs2.GT.smooth_coef)) then 
	    U1=U2
            IT=IT+1
            GOTO 22
	  endif 
	enddo
      endif
      IF(T.LT.time_loop+STEP-1.D-12)THEN
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
      ELSEIF(T.GT.time_loop+STEP+1.D-12)THEN
          TAU=(time_loop+STEP-T1)
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

