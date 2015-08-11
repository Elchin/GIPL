! Geophysical Instatitue Permafrost Laboratory 2 model version 01 (GIPL2v01)
! GIPL2 is a numerical transient model that employs phase changes and the effect of the unfrozen volumetric  water content in the non-homogeniuos soil texture
! Original version of the model developed by Romanovsky and Tipenko 2004 and described in Marchenko et al., (2008)
! Current version been significanlty modefied from its predicessor and using the IRF coding design
! This version is maintained by E. Jafarov at INSTAAR, CU Boulder
! Please site Jafarov et al., (2012) work when using it.

program gipl

use bnd
use thermo
use grd
use alt

implicit none

call initialize
call run_model
call finilize

end

subroutine run_model
use bnd
use thermo
use grd
use alt

implicit none
! all function names start with letter 'f'
    external futemp,ftcon,fsat_unf_water,fapp_hcap
    
! functions
    real*8 ftcon                                            ! thermal conductivity
    real*8 fapp_hcap                                        ! apparent heat capacity
    real*8 futemp                                           ! temperature interpolation
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
    integer :: i_site,j_time,i_grd,i_lay

    real*8, allocatable :: z(:)                             ! vertical grid
    real*8 :: hcscale                                       ! nondimensional factor
    integer :: ierr


    !call initialize
    !call run_model

    allocate(z(n_grd),STAT=ierr)
    allocate(dz(n_grd),STAT=ierr)

    z=zdepth/zdepth(n_grd)                              ! normalizing soil grid
    do i_grd=2,n_grd
        dz(i_grd)=z(i_grd)-z(i_grd-1)
    enddo
    hcscale=zdepth(n_grd)*zdepth(n_grd)/n_sec_day
    hcap_frz=hcap_frz*hcscale
    hcap_thw=hcap_thw*hcscale
    hcap_s=hcap_snow*hcscale
    L_fus=hcscale*333.2*1.D+6

    call  assign_layer_id(n_lay,n_lay_cur,n_site,n_grd,zdepth,n_bnd_lay,lay_id)

    allocate(res_save(m_grd+3,n_site),STAT=ierr)
    allocate(dfrz_frn(n_time),STAT=ierr)
    allocate(time_loop(n_site),STAT=ierr)
    allocate(time_cur(n_site),STAT=ierr)
    allocate(res_vars(n_time,m_grd+3),STAT=ierr)
!________________INITIALIZATION________________________
!      n_itime=n_time*time_end+3
    n_itime=n_time+2
    time_s=time_step*DBLE(n_time*time_beg)
    time_e=time_step*DBLE(n_time*time_end)
    time_loop=0.0D0
    time_internal=0.0D0
!-----------------------------------------------initial data read

    call init_cond(restart,n_site)

    allocate(utemp_time_i(n_itime),STAT=ierr)                  ! allocating interval varialbe after interation
    allocate(utemp_i(n_itime,n_site),STAT=ierr)
    allocate(snd_i(n_itime,n_site),STAT=ierr)
    allocate(stcon_i(n_itime,n_site),STAT=ierr)
	
!****************************************************************************	
    do j_time=1,n_time+2
        utemp_time_i(j_time)=time_restart+DBLE(j_time-1)*time_step
    enddo
    i_time=1

    do i_site=1,n_site
        if (lbound.EQ.2)temp_grd(i_site)=temp_grd(i_site)*zdepth(n_grd)
        do i_lay=1,n_lay_cur(i_site)
            temp_frz(i_lay,i_site)=-(vwc(i_lay,i_site)/a_coef(i_lay,i_site))**(1.d0/b_coef(i_lay,i_site))
        enddo
        call interpolate(utemp_time,utemp(:,i_site),n_temp,utemp_time_i,utemp_i(:,i_site),n_itime)
        call interpolate(snd_time,snd(:,i_site),n_snow,utemp_time_i,snd_i(:,i_site),n_itime)
        call snowfix(utemp_i(:,i_site),snd_i(:,i_site),n_itime)
        call interpolate(stcon_time,stcon(:,i_site),n_stcon,utemp_time_i,stcon_i(:,i_site),n_itime)
        call active_layer(i_site)
    enddo

    open(1,file=result_file,STATUS='unknown')
    open(2,file=aver_res_file,STATUS='unknown')
    open(3,file=restart_file,STATUS='unknown')

! begining of the major loop
65 CONTINUE
    do i_site=1,n_site
        time_cur(i_site)=time_loop(i_site)+time_restart
        call save_results(i_site,time_cur(i_site),time_loop(i_site))
        6666  continue
        call stefan1D(temp(i_site,:),n_grd,dz,time_loop(i_site),i_site,lay_id(i_site,:), &
                    temp_grd(i_site),futemp,fapp_hcap,ftcon,fsat_unf_water)
        time_loop(i_site)=time_loop(i_site)+time_step
        time_cur(i_site)=time_loop(i_site)+time_restart
        if(i_time(i_site).LT.n_time)  then
            i_time(i_site)=i_time(i_site)+1
            call save_results(i_site,time_cur(i_site),time_loop(i_site))
            call active_layer(i_site)
            GOTO 6666
        endif
        if(time_s.LT.time_e.AND.time_loop(1).GT.time_s)then
            do j_time=1,n_time			! WRITTING RESULTS
                write(1,FMT1) i_site, (res_vars(j_time,i_grd),i_grd=1,m_grd+3)
            enddo
        endif
        do i_grd=1,m_grd+3
            res_save(i_grd,i_site)=sum((res_vars(:,i_grd)))
        enddo
     enddo

      i_time=1
      do i_site=1,n_site
           frz_up_time_cur=-7777.D0
           frz_up_time_tot=frz_up_time_cur
           do j_time=2,n_time
              if((n_frz_frn(j_time,i_site)-n_frz_frn(j_time-1,i_site)).EQ.-2)then
                if(z_frz_frn(j_time-1,n_frz_frn(j_time-1,i_site),i_site).GE.frz_frn_min) frz_up_time_cur=SNGL(res_vars(j_time,1))
              endif
      enddo

      if(frz_up_time_cur.GT.0.0)then
            frz_up_time_tot=AMOD(frz_up_time_cur,REAL(n_time))
            if(frz_up_time_tot.EQ.0.0)frz_up_time_tot=REAL(n_time)
      endif

      dfrz_frn=z_frz_frn(:,1,i_site)
      call save_results(i_site,time_cur(i_site),time_loop(i_site))
      call active_layer(i_site)

    !____WRITTING MEAN
      write(2,FMT2) i_site,(res_save(i_grd,i_site)/DBLE(n_time),i_grd=1,m_grd+3), &
                                dfrz_frn(n_time),frz_up_time_cur,frz_up_time_tot
      do j_time=1,n_itime
            utemp_time_i(j_time)=time_cur(1)+DBLE(j_time-1)*time_step
      enddo
      call interpolate(utemp_time,utemp(:,i_site),n_temp,utemp_time_i,utemp_i(:,i_site),n_itime)
      call interpolate(snd_time,snd(:,i_site),n_snow,utemp_time_i,snd_i(:,i_site),n_itime)
      call snowfix(utemp_i(:,i_site),snd_i(:,i_site),n_itime)
      call interpolate(stcon_time,stcon(:,i_site),n_stcon,utemp_time_i,stcon_i(:,i_site),n_itime)
    enddo

    rewind(3) ! -------------start file writting begin
    write(3, * ) time_restart
!   write(3, * ) time_cur(1)
    do i_grd=1,n_grd
        write (3,* ) ( temp(i_site,i_grd),i_site=1,n_site)
    enddo     ! -------------start file writting end

    time_internal=time_loop(1)

    if(time_loop(1).LT.time_e)GO TO 65
    close(1);close(2);close(3)

   ! call finilize

end subroutine run_model

subroutine initialize

use bnd
use thermo
use grd
use alt

implicit none

    integer IREAD,ierr
    integer :: i,j,k,z_num


    real*8 ,allocatable ::gtzone(:,:)
    character*64 stdummy
    character*64 fconfig

    character*64 finput,boundf,snowf,rsnowf,inif
    character*64 cmdf,cetkaf,fvegetation,fgeology

    real*8,allocatable:: A1(:,:),A2(:,:),A3(:,:),A4(:,:),A5(:,:)
    real*8,allocatable:: A6(:,:),A7(:,:),A8(:,:)
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
! input files
    read(60,'(A)')finput
    read(60,'(A)')boundf
    read(60,'(A)')snowf
    read(60,'(A)')rsnowf
    read(60,'(A)')inif
    read(60,'(A)')cmdf
    read(60,'(A)')cetkaf
    read(60,'(A)')fvegetation
    read(60,'(A)')fgeology
!output files
    read(60,'(A)')aver_res_file
    read(60,'(A)')result_file
    read(60,'(A)')restart_file
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
    allocate(snow_code(n_site),STAT=ierr)
    allocate(veg_code(n_site),STAT=ierr)
    allocate(geo_code(n_site),STAT=ierr)
    allocate(gt_zone_code(n_site),STAT=ierr)
    allocate(temp_grd(n_site),STAT=ierr)
    do i=1,n_site
        read(60,*) IREAD,snow_code(i),veg_code(i),geo_code(i),&
                gt_zone_code(i),temp_grd(i)
    enddo
close(60)
!   print*, trim(finput),' has been read'
      
open(60,file=BOUNDF)
    read(60,*)n_temp
    allocate(utemp_time(n_temp),STAT=ierr)
    allocate(utemp(n_temp,n_site),STAT=ierr)
    do i=1,n_temp
        read(60,*) utemp_time(I),(utemp(I,J),J=1,n_site)
    enddo
close(60)
!   print*,trim(boundf),' has been read'

open(60,file=RSNOWF)
    read(60,*)n_stcon
    allocate(stcon_time(n_stcon),STAT=ierr)
    allocate(stcon(n_stcon,n_site),STAT=ierr)
    do i=1,n_stcon
        read(60,*) stcon_time(i),(stcon(i,J),J=1,n_site)
    enddo
close(60)
!   print*,trim(rsnowf),' has been read'

open(60,file=SNOWF)
    read(60,*)n_snow
    allocate(snd_time(n_snow),STAT=ierr)
    allocate(snd(n_snow,n_site),STAT=ierr)
    do I=1,n_snow
       read(60,*) snd_time(i),(snd(i,J),J=1,n_site)
    enddo
close(60)
!   print*,trim(snowf),' has been read'

open(60,file=inif,action='read')
    read(60,*)z_num,n_ini!,time_restart
    allocate(zdepth_ini(n_ini),STAT=ierr)
    allocate(ztemp_ini(n_ini,n_site),STAT=ierr)
    allocate(gtzone(n_ini,z_num+1),STAT=ierr)
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
    read(60,*)time_step,TAUM,time_step_min
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
    allocate(zdepth(n_grd),STAT=ierr)
    do i=1,n_grd
        read(60,*) zdepth(i)
    enddo
    read(60,*)m_grd
    allocate(zdepth_id(m_grd),STAT=ierr)
    do j=1,m_grd
        read(60,*)zdepth_id(j)
    enddo
close(60)
!   print*,trim(cetkaf),' has been read'

! note: that all max n_lay_cur layers has to be read or it will a give segmantation error
!      n_lay=10!MAXVAL(n_lay_cur)
!----------------------------------------------------      
open (60, file=fvegetation)
    read(60,*) vln ! reads numbers of  classes
    allocate(A1(n_lay,vln),STAT=ierr) ! vwc
    allocate(A2(n_lay,vln),STAT=ierr) ! a_coef
    allocate(A3(n_lay,vln),STAT=ierr) ! b_coef
    allocate(A4(n_lay,vln),STAT=ierr) ! hcap_frz
    allocate(A5(n_lay,vln),STAT=ierr)  !hcap_thw
    allocate(A6(n_lay,vln),STAT=ierr)  !tcon_frz
    allocate(A7(n_lay,vln),STAT=ierr)  !tcon_thw
    allocate(A8(vln,n_lay),STAT=ierr)  !bot_cond
    allocate(veg_class(vln),STAT=ierr) !veg_class
    allocate(num_vl(vln),STAT=ierr)  !num_vl number of vegetation layers
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
    allocate(B1(n_lay,gln),STAT=ierr) ! vwc
    allocate(B2(n_lay,gln),STAT=ierr) ! a_coef
    allocate(B3(n_lay,gln),STAT=ierr) ! b_coef
    allocate(B4(n_lay,gln),STAT=ierr) ! hcap_frz
    allocate(B5(n_lay,gln),STAT=ierr)  !hcap_thw
    allocate(B6(n_lay,gln),STAT=ierr)  !tcon_frz
    allocate(B7(n_lay,gln),STAT=ierr)  !tcon_thw
    allocate(B8(gln,n_lay),STAT=ierr)  !bot_cond
    allocate(geo_class(gln),STAT=ierr) !geo_class
    allocate(num_gl(gln),STAT=ierr)  !num_vl number of lithologic layers
    do I = 1,gln
        read(60,*)geo_class(i),num_gl(i)
        do j=1,num_gl(i)
        read(60,*)B1(J,I),B2(J,I),B3(J,I), &
                B4(J,I),B5(J,I),B6(J,I),B7(J,I),B8(I,J)
        enddo
    enddo
close(60)
!      print*,trim(fgeology),' has been read'

allocate(vwc(n_lay,n_site),STAT=ierr)
allocate(a_coef(n_lay,n_site),STAT=ierr)
allocate(b_coef(n_lay,n_site),STAT=ierr)
allocate(EE(n_lay,n_site),STAT=ierr)
allocate(hcap_frz(n_lay,n_site),STAT=ierr)
allocate(hcap_thw(n_lay,n_site),STAT=ierr)
allocate(tcon_frz(n_lay,n_site),STAT=ierr)
allocate(tcon_thw(n_lay,n_site),STAT=ierr)
allocate(n_lay_cur(n_site),STAT=ierr)
allocate(n_bnd_lay(n_site,n_lay+1),STAT=ierr)

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
	n_lay_cur(I)=num_vl(veg_code(i))+num_gl(geo_code(i))! maximum number of soil layer = organic layers + mineral layers
 	do j=num_vl(veg_code(i))+1,n_lay_cur(I)
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
    n_bnd_lay(i,n_lay_cur(I)+1)=zdepth(n_grd)
enddo

!deallocate(utemp_time)

allocate(temp(n_site,n_grd),STAT=ierr)
allocate(lay_id(n_site,n_grd),STAT=ierr)
allocate(i_time(n_site),STAT=ierr)
allocate(z_frz_frn(n_time,n_frz_max,n_site),STAT=ierr)
allocate(n_frz_frn(n_time,n_site),STAT=ierr)
allocate(temp_frz(n_lay,n_site),STAT=ierr)
! setting formats for the results and mean files
write(FMT1,'(A30,I0,A12)')'(1x,I10,1x,F12.3,2(1x,F16.12),',m_grd,'(1x,F16.12))'
write(FMT2,'(A28,I0,A40)')'(1x,I10,1x,F12.3,2(1x,F8.3),',m_grd,'(1x,F8.3),(1x,F8.3,1x,F12.3),(1x,F12.3))'

end subroutine initialize

subroutine finilize
    use bnd
    use thermo
    use grd
    use alt
    implicit none

    deallocate(utemp_time)
    deallocate(utemp)
    deallocate(utemp_time_i)
    deallocate(utemp_i)
    deallocate(snd_time)
    deallocate(snd)
    deallocate(stcon_time)
    deallocate(stcon)
    deallocate(snd_i)
    deallocate(stcon_i)
    deallocate(vwc)
    deallocate(a_coef)
    deallocate(b_coef)
    deallocate(temp_frz)
    deallocate(EE)
    deallocate(hcap_frz)
    deallocate(hcap_thw)
    deallocate(tcon_frz)
    deallocate(tcon_thw)
    deallocate(temp)
    deallocate(n_bnd_lay)
    deallocate(snow_code)
    deallocate(veg_code)
    deallocate(geo_code)
    deallocate(gt_zone_code)
    deallocate(temp_grd)
    deallocate(res_vars)
    deallocate(n_lay_cur)
    deallocate(lay_id)
    deallocate(zdepth_id)
    deallocate(zdepth_ini)
    deallocate(ztemp_ini)
    deallocate(n_frz_frn)
    deallocate(i_time)
    deallocate(z_frz_frn)

end subroutine finilize

subroutine init_cond(rest_flag,nsite)

use bnd
use thermo
use grd

implicit none
   integer rest_flag,nsite
   integer isite,igrd
   character*64 inif

if(rest_flag.EQ.1)then                      !restart=1 means reading initial data from
    do isite=1,nsite
        call interpolate(zdepth_ini,ztemp_ini(:,isite),n_ini,zdepth,temp(isite,:),n_grd)
    enddo
elseif(rest_flag.EQ.0)then                  !restart=0 enbales spinup
    write(INIF,'(A11,I3.3,A4)') 'dump/start_',nsite+1,'.txt'
    open(60,file=INIF,action='READ')
        read(60,*)time_restart              ! day number in restart file
        do igrd=1,n_grd
            read (60,* ) ( temp(isite,igrd),isite=1,nsite)
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

subroutine save_results(isite,time1,time2)
use thermo
use grd
use alt

implicit none
    integer :: isite,igrd
    real*8 :: time1,time2
    real*8 :: futemp,fsnow_level

    res_vars(i_time(isite),1)=time1
    res_vars(i_time(isite),2)=futemp(time2,isite)
    res_vars(i_time(isite),3)=fsnow_level(isite,time2)
    do  igrd=1,m_grd
        res_vars(i_time(isite),igrd+3)=temp(isite,zdepth_id(igrd))
    enddo

end subroutine save_results

!________________________________________________
!__________________FUNCTIONS_____________________
!________________________________________________
real*8 function funf_water(T,NNN,I)
use thermo
implicit none
    real*8, intent(in) :: T ! temprature
    integer, intent(in) :: NNN, I
    real*8 :: temp_dep
    real*8 :: a,b,e
    real*8 :: theta

temp_dep=temp_frz(NNN,I) ! change I to k0 everywhere except temp_dep
e=EE(NNN,I)
theta=vwc(NNN,I)
a=a_coef(NNN,I)
b=b_coef(NNN,I)
	
IF(T.LE.temp_dep-e)THEN
    funf_water=a*((DABS(T))**b)
ELSEIF(T.GT.temp_dep)THEN
    funf_water=theta
ELSE
    funf_water=a*((DABS(temp_dep-e))**b)
    funf_water=funf_water+(theta-funf_water)*(T+e-temp_dep)/e
endif
return
end function funf_water
!-----------------------------------------------
real*8 function fsat_unf_water(T,NNN,I)!Saturated unforzen water
use thermo

!IMPLICIT REAL*8(A-H,O-Z)
implicit none
    real*8, intent(in) :: T
    integer, intent(in) :: NNN, I
    real*8 :: temp_dep
    real*8 :: a,b,e
    real*8 :: theta
      
temp_dep=temp_frz(NNN,I) ! freezing temprature depression
e=EE(NNN,I)
theta=vwc(NNN,I)
a=a_coef(NNN,I)
b=b_coef(NNN,I)
IF(T.LE.temp_dep-e)THEN
    fsat_unf_water=a*((DABS(T))**b)
ELSEIF(T.GT.temp_dep)THEN
    fsat_unf_water=theta
ELSE
    fsat_unf_water=a*((DABS(temp_dep-e))**b)
    fsat_unf_water=fsat_unf_water+(theta-fsat_unf_water)*(T+e-temp_dep)/e
ENDIF
fsat_unf_water=fsat_unf_water/theta
return
end function fsat_unf_water
!----------------------------------------
real*8 function fdunf_water(T,NNN,I)
use thermo
!IMPLICIT REAL*8(A-H,O-Z)
implicit none
    real*8, intent(in) :: T ! temprature
    integer, intent(in) :: NNN, I
    real*8 :: temp_dep
    real*8 :: a,b,e
    real*8 :: theta

temp_dep=temp_frz(NNN,I)
e=EE(NNN,I)
theta=vwc(NNN,I)
a=a_coef(NNN,I)
b=b_coef(NNN,I)

if(T.LE.temp_dep-e)THEN
    fdunf_water=-b*a*((DABS(T))**(b-1.0D0))
elseif(T.GT.temp_dep)THEN
    fdunf_water=0.0D0
else
    fdunf_water=a*((DABS(temp_dep-e))**b)
    fdunf_water=(b-fdunf_water)/e
endif
return
end function fdunf_water

!----------------------------------------
real*8 function futemp(T,I)
use bnd
implicit none
    real*8 T
    integer I,II

II=1+IDINT((T-time_internal)/time_step)
futemp=utemp_i(II,I)+(T+time_restart-utemp_time_i(II)) &
    *(utemp_i(II+1,I)-utemp_i(II,I))/(utemp_time_i(II+1)-utemp_time_i(II))
return
end function futemp
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
subroutine interpolate(XIN,YIN,NIN,XOUT,YOUT,n_itime)
! Linear interpolation
real*8, intent(in) :: XIN(NIN),YIN(NIN)
real*8, intent(out) :: XOUT(n_itime),YOUT(n_itime)
integer :: NIN,n_itime
do I=1,n_itime
    if(XOUT(I).LE.XIN(1))THEN
        YOUT(I)=YIN(1)
        GOTO 1
    elseif(XOUT(I).GT.XIN(NIN))THEN
        YOUT(I)=YIN(NIN)
        GOTO 1
    else
        do J=1,NIN-1
            if (XIN(J).LT.XOUT(I).AND.XOUT(I).LE.XIN(J+1))THEN
                YOUT(I)=YIN(J)+(XOUT(I)-XIN(J))*(YIN(J+1)-YIN(J))/(XIN(J+1)-XIN(J))
                GOTO 1
            endif
        enddo
    endif
    1 continue
enddo
return
end

!----------------------------------------
subroutine assign_layer_id(n_lay,n_lay_cur,n_site,n_grd,zdepth,n_bnd_lay,lay_id)
!assigns correspond layer id to the grid point
!starting from surface to the bottom
implicit none

    integer :: n_site,n_grd,n_lay
    integer :: lay_id(n_site,n_grd),n_lay_cur(n_site)
    real*8 :: zdepth(n_grd)
    real*8 :: n_bnd_lay(n_site,n_lay+1)
    integer :: isite,igrd,ilay

do isite=1,n_site
    do 6 igrd=1,n_grd
        lay_id(isite,igrd)=n_lay_cur(isite)
        do ilay=1,n_lay_cur(isite)-1
            if ( n_bnd_lay(isite,ilay).LE.zdepth(igrd).AND.zdepth(igrd).LT.n_bnd_lay(isite,ilay+1))then
                lay_id(isite,igrd)=ilay
                GOTO 6
            endif
        enddo
    6 continue
enddo
return
end

real*8 function fsnow_level(site_id,time)
use bnd
real*8 :: time
integer :: site_id,II

II=1+IDINT((time-time_internal)/time_step)
fsnow_level=snd_i(II,site_id)+(time+time_restart-utemp_time_i(II))* &
            (snd_i(II+1,site_id)-snd_i(II,site_id))/(utemp_time_i(II+1)-utemp_time_i(II))
return
end function fsnow_level

real*8 function ftcon(T,id,j,time_cur)
use bnd
use grd
use thermo

implicit real*8(A-H,O-Z)
integer :: II
      
gr_sur=sea_level
dsnow=sea_level-fsnow_level(id,time_cur)
NS=lay_id(id,j)
if(zdepth(j).le.dsnow)then                                  !atmosphere
    ftcon=1.d4
elseif (zdepth(j).Lt.gr_sur)then                            !snow
    II=1+IDINT((time_cur-time_internal)/time_step)
        ftcon=stcon_i(II,id)+(time_cur+time_restart-utemp_time_i(II))* &
        (stcon_i(II+1,id)-stcon_i(II,id))/(utemp_time_i(II+1)-utemp_time_i(II))
else                                                        !ground
    WC=funf_water(T,NS,id)/vwc(NS,id)
    ftcon=(tcon_thw(NS,id)**WC)*(tcon_frz(NS,id)**(1.0-WC))
endif
return
end function ftcon

!----------------------------------------
real*8 function fhcap(T,NNUS,I)
use thermo

IMPLICIT REAL*8(A-H,O-Z)
DIMENSION NNUS(2),T(2)
        
H=1/(T(1)-T(2))
if(DABS(T(1)-T(2)).LT.1.D-6) THEN
    fhcap=0.5d0*(fdunf_water(T(1),NNUS(1),I)+fdunf_water(T(2),NNUS(2),I))
else
    if (nnus(1).ne.nnus(2))THEN
        fhcap=0.5D0*( H*(funf_water(T(1),NNUS(1),I)-funf_water(T(2),NNUS(1),I))+ &
                H*(funf_water(T(1),NNUS(2),I)-funf_water(T(2),NNUS(2),I)) )
    else
        fhcap=H*(funf_water(T(1),NNUS(1),I)-funf_water(T(2),NNUS(2),I))
    endif
endif
fhcap=L_fus*DABS(fhcap)
return
end function fhcap
!----------------------------------------
!----------------------------------------
real*8 function fapp_hcap(T,I,J)       ! Apparent heat capacity
use thermo
use grd

implicit real*8(A-H,O-Z)
DIMENSION T(n_grd),WW(2),NN(2)

li=lay_id(I,J)                          ! layer index
gr_sur=sea_level                        ! ground surface
if(zdepth(J).lE.gr_sur)then
    fapp_hcap=hcap_s                  ! heat capacity for snow
else
    WC=funf_water(T(J),li,I)/vwc(li,I)
    fapp_hcap=hcap_thw(li,I)*WC+hcap_frz(li,I)*(1.0-WC)
    if(J.GT.(1).AND.J.LT.n_grd)then
        WW(1)=(T(J-1)+T(J))/2.D0
        NN(1)=lay_id(I,J-1)
        WW(2)=T(J)
        NN(2)=lay_id(I,J)
        fapp_hcap=fapp_hcap+fhcap(WW,NN,I)*dz(J)/(dz(J+1)+dz(J))
        WW(1)=T(J)
        NN(1)=lay_id(I,J)
        WW(2)=(T(J+1)+T(J))/2.D0
        NN(2)=lay_id(I,J+1)
        fapp_hcap=fapp_hcap+fhcap(WW,NN,I)*dz(J+1)/(dz(J+1)+dz(J))
    elseif(J.EQ.1)then
        WW(1)=T(J)
        NN(1)=lay_id(I,J)
        WW(2)=(T(J+1)+T(J))/2.D0
        NN(2)=lay_id(I,J+1)
        fapp_hcap=fapp_hcap+fhcap(WW,NN,I)
    elseif(J.EQ.n_grd)then
        WW(1)=(T(J-1)+T(J))/2.D0
        NN(1)=lay_id(I,J-1)
        WW(2)=T(J)
        NN(2)=lay_id(I,J)
        fapp_hcap=fapp_hcap+fhcap(WW,NN,I)
    endif
endif

return
end
!-------------------------------------------------------
subroutine stefan1D(temps,n_grd,dz,time_loop,isite,lay_idx, &
flux,futemp,fapp_hcap,ftcon,fsat_unf_water)

!call stefan1D(temp(i_site,:),n_grd,dz,time_loop(i_site),i_site,lay_id(i_site,:), &
!temp_grd(i_site),futemp,fapp_hcap,ftcon,fsat_unf_water)

    use thermo
    use bnd

    implicit none

    integer, intent(inout) :: n_grd
    real*8, intent(inout) :: dz(n_grd),temps(n_grd)
    integer, intent(inout) :: lay_idx(n_grd)
    real*8, intent(inout) :: time_loop
    real*8 :: futemp,flux,fapp_hcap,ftcon,fsat_unf_water

    integer :: isite,i_grd,IT

! tridiagonal variables
    real*8 :: RAB1,RAB2,AKAPA2,AMU2,Q2
    real*8 :: A,B,C,D
    real*8 :: ALF(n_grd),BET(n_grd)
    real*8 :: EEY,EEY1,abs1,abs2

    real*8 :: temp_o(n_grd)             ! old temperature before tridiagonal method
    real*8 :: temp_n(n_grd)             ! new temperature after tridiagonal method

! time counter internal to this subroutine
    real*8 :: time_l                    ! loop time in a subroutine
    real*8 :: time_p                    ! present time in a subroutine
    real*8 :: timei                     ! main subroutine timer
    real :: time_swith                  ! for timei

    time_l=time_loop
    time_swith=-1.0
    timei=TAUM
    temps=temp(isite,:)
64  continue
    time_p=time_l+timei
    temp_o=temps
    IT=1
    ALF(2)=0.D0
    BET(2)=futemp(time_p,isite)
22  continue
    if(IT.GT.ITMAX) then
        timei=timei/2.D0
        time_swith=-1.0
        GOTO 64
    endif

    do i_grd=2,n_grd-1
        D=fapp_hcap(temp_o,isite,i_grd)/timei
        A=2.D0*ftcon(temp_o(i_grd),isite,i_grd,time_p)/(dz(i_grd)*(dz(i_grd)+dz(i_grd+1)))
        B=2.D0*ftcon(temp_o(i_grd+1),isite,i_grd+1,time_p)/(dz(i_grd+1)*(dz(i_grd)+dz(i_grd+1)))
        C=A+B+D
        ALF(i_grd+1)=B/(C-A*ALF(i_grd))
        BET(i_grd+1)=(A*BET(i_grd)+D*temps(i_grd))/(C-A*ALF(i_grd))
    enddo

!    RAB1=ftcon(temp_o(n_grd),isite,n_grd,time_cur)
    RAB1=ftcon(temp_o(n_grd),isite,n_grd,time_p)
    RAB2=fapp_hcap(temp_o,isite,n_grd)
    AKAPA2=2.D0*RAB1/(((RAB2*dz(n_grd)*dz(n_grd))/timei+2.D0*RAB1))
    Q2=RAB1*flux
    AMU2=(temps(n_grd)*RAB2/timei+2.D0*Q2/dz(n_grd))/(RAB2/timei+2.D0*RAB1 &
                                                        /dz(n_grd)**2.D0)
    if(DABS(AKAPA2)>1.D0) then
        print*,'Tridiagonal method is failed - chang you time step tau'
        print*,rab1,rab2,akapa2
        STOP
    endif

! assigns boundary condition check
    if (lbound.EQ.2)then
        temp_n(n_grd)=(AMU2+AKAPA2*BET(n_grd))/(1.D0-ALF(n_grd)*AKAPA2)
    else
        temp_n(n_grd)=flux
    endif

! calculates new tempratures
    do i_grd=1,n_grd-1
        temp_n(n_grd-i_grd)=ALF(n_grd-i_grd+1)*temp_n(n_grd-i_grd+1)+BET(n_grd-i_grd+1)
    enddo

    if(timei>time_step_min) then
        do i_grd=1,n_grd
            EEY=fsat_unf_water(temp_n(i_grd),lay_idx(i_grd),isite)
            EEY1=fsat_unf_water(temp_o(i_grd),lay_idx(i_grd),isite)
            abs1=DABS(EEY-EEY1)
            abs2=DABS(temp_o(i_grd)-temp_n(i_grd))
            if((abs1.GT.unf_water_coef).or.(abs2.GT.smooth_coef)) then
                temp_o=temp_n
                IT=IT+1
                GOTO 22
            endif
        enddo
    endif

    if(time_p.LT.time_loop+time_step-1.D-12)then
        time_l=time_p
        temps=temp_n
        if(time_swith>0) then
            if(timei.LT.TAUM) then
                timei=timei*2.D0
                time_swith=-1.0
            endif
        else
            time_swith=1.0
        endif
        GOTO 64
        elseif(time_p.GT.time_loop+time_step+1.D-12)then
            timei=(time_loop+time_step-time_l)
            goto 64
        else
            temps=temp_n
    endif

end subroutine stefan1D


subroutine filexist(filename)
implicit none
    character*64 filename
    logical chf

inquire(file=filename,exist=chf)
if (.not.chf) then
    write(*,'(/'' FILE '',a, '' DOESNT EXIST'')')trim(filename)
stop
endif

end subroutine filexist!-----------------------------------------------

