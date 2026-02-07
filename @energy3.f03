!-----------------------------------------------------------------
!% mpif90 -mcmodel=medium -fpic -O2 -o a.out @energy3.f03 &> log
!% mpiexec -n 5 a.out &
!% pspdf sai105.231g -- make plots by using "sai231g.ps"
!% In Windows 11, plot sai231g.pdf by sai231g.pdf
!-----------------------------------------------------------------
      program spin_dynamics 
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'mpif.h'
      include    'param-spinRL7.h'
!
      integer(C_INT) rank,size,igrp,kstart,np1,np2,            &
                 i1(0:num_proc),i2(0:num_proc),i3(0:num_proc), &
                 i4(0:num_proc),                               &
                 cnt_recv(0:num_proc),disp_recv(0:num_proc),   &
                 cnt_recvC(0:num_proc),disp_recvC(0:num_proc), &
                 cnt_send,ierror,mpierror,i_MC,ifcmp,kk
!
      real(C_DOUBLE) e_sp0,e_c_r0,e_LJ0,e_sp1,e_c_r1,e_LJ1
      logical  if_LJ0
!
      real(C_DOUBLE) x,y,z,spx,spy,spz,ch,fx,fy,fz,fxC,fyC,fzC, &
                     Jint,r_ij,rintC,fc1,fc2,fcLJ
      common/partcl1/ &
                 x(np0),y(np0),z(np0),spx(np0),spy(np0),spz(np0), &
                 ch(np0)
      common/partcl2/ &
                 fx(np0),fy(np0),fz(np0),fxC(np0),fyC(np0), &
                 fzC(np0),Jint(nbx,np0),r_ij(nbx,np0),rintC(nbx,np0)
!
      real(C_DOUBLE) vx(np0),vy(np0),vz(np0),mass(np0),ag(np0),ep(np0),&
                 fkx(np0),fky(np0),fkz(np0),sp2(np0),spx0(np0),     &
                 spy0(np0),spz0(np0),spx1(np0),spy1(np0),spz1(np0), &
                 wsp1(3,np0),wsp2(3,np0),Jint_1,           &
                 u,u1,u2,spx00(np0),spy00(np0),spz00(np0), &
                 u_1,u_2,aspx(np0),aspy(np0),aspz(np0)
!
      integer(C_INT) nintS,lintS,nintC,lintC,if_LJ,spec,site
      common/partcl3/ & 
                 nintS(np0),lintS(nbx,np0),nintC(np0),lintC(nbx,np0), &
                 if_LJ(np0),spec(np0),site(np0)
!
      integer(C_INT) nlist(np0),lmax,k0,iac1,iac2,irej,modes,kwrite, &
                 n_of_j,np10,np100,ia,ja,ka,i1x,i2x,i1y,i2y,i1z,i2z, &
                 n_MCsteps,nt_p3m
!
      real(C_DOUBLE) t8,dt,dts,dt0,dth,Bex,Bey,Bez,spin2,spin3, &
                 Jaa,Jbb,Jab,J00,B00,Bapx,Bapy,Bapz,Bap,        &
                 tau_b,tau_diss,t_adv,fw,fw00,                  &
                 Temp,TCurie,g,mue_B,hbar,KJoule,Kcal,mol,kT,eV,&
                 omg_b,rnd,                                     &
                 f1,g1,b1,prb,qsx,qsy,qsz,rsx,rsy,rsz,hh1,hh2,  &
                 qq,rqq,dthg,dthqq,rdiv,spmin,                  &
                 alp,xx,yy,zz,rr,rr0,toler,deps,Usys8(nhs),     &
                 ss,smax,sav,t_unit,e_unit,a_unit,m_unit,m_Fe,m_O, &
                 pi,pi2,th,ph,tmax,tmax0,cptot,unif1(2),unif2(2),  &
                 xleng0,yleng0,zleng0,sss,rlist(np0),           &
                 buffer1(4),buffer2(4),del_en,wtime,            &
                 fc3,J_ki,wfdt,vth0,vth_O,vth_Fe,vth,           &
                 svx,svy,svz,sqrt2,vmax1,dgaus2,vsq1,vsq2,      &
                 vx1,vy1,vz1,mas,e_sp,e_c_r,e_LJ,e_Coulomb_p3m
!
      real(C_DOUBLE) x0(np0),y0(np0),z0(np0),vx0(np0),vy0(np0),vz0(np0)
      real(C_DOUBLE) x_0(np0),y_0(np0),z_0(np0)
      real(C_DOUBLE) dtwr,dtwr2,cputime
!
      real(C_DOUBLE) rad_Fe,rad_O,elj_Fe,elj_O,rcut,rcutC
      common/atoms/  rad_Fe,rad_O,elj_Fe,elj_O,rcut,rcutC
!
      integer(C_INT) it,is,iw,iwa,iwb,istop,                  &
                    iwrt1,iwrt2,iwrta,iwrtb,iter,itermax,     &
                    i00,i,j,k,l,nsite,notconv,nframe,mx,my,mz,&
                    ir0,lsite(num_proc,np0),ic,nstep_MC,nstep_MCt, &
                    ifdt,kmax,if_vres(np0),ifv
      common/parm1/ it,is
      common/parm2/ dtwr,dtwr2
      common/parm3/ t8,pi,dt,tmax,cptot
      common/parm4/ mx,my,mz
      common/parm5/ n_MCsteps,nt_p3m
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      common/imemo/ iwa,iwb
      common/ranfff/ ir0
      common/itera/  toler,itermax
!
      real(C_float) spinx,spinz,spin7,Bextx,Bextz,magx,magy,magz, &
                    Usys,conv,aitr,psdt,tfix,uss,usb,tsx,tsy,tsz,sum_mb,&
                    U_Fe,U_O,ds_Fe,ds_O,fdt4,vdt4,timeh,dtrhs
      real(C_float) sx1(3),sy1(3),sz1(3),sx2(3),sy2(3),sz2(3),sn1(3), &
                    csx,csy,csz,axis(100),freq(100),t4,Bex4,Bez4, &
                    spx4(np0),spy4(np0),spz4(np0),    &
                    ch4(np0),x4(np0),y4(np0),z4(np0), &
                    vx4(np0),vy4(np0),vz4(np0)
      real(C_DOUBLE) wx,wy,wz,wn,wsx,wsy,wsz,wp, &
                    wx1,wy1,wz1,wn1,uav,wt1,wx7,wy7,wz7,wn7,uav7,wt7,  &
                    ssx,ssy,ssz,tht,phi,psi,bbw,tsz0,atsz,av_tsz(nhs), &
                    ub,um,ub1,um1,uu1(3),uu2(3),ranff,fdt8,vdt8,ds1,ds2
      common/ehist/ spinx(nhs),spinz(nhs),spin7(nhs), &
                    Bextx(nhs),Bextz(nhs),magx(nhs),magy(nhs),magz(nhs),&
                    Usys(nhs),conv(nhs),aitr(nhs),psdt(nhs),tfix(nhs),  &
                    uss(nhs),usb(nhs),tsx(nhs,3),tsy(nhs,3),tsz(nhs,3), &
                    sum_mb(nhs),U_Fe(nhs),U_O(nhs),ds_Fe(nhs),ds_O(nhs),&
                    fdt4(nhs),vdt4(nhs),timeh(nhs)
!
      logical       MC_first/.true./,ft06_start/.true./
!
      real(C_DOUBLE) alpha,xleng,yleng,zleng
      integer(C_INT) PP
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
!
      character*8    label,cdate*10,ctime*8,plot_ch*8,char*2
      real(C_float)  time
      common/headr1/  label,cdate
      common/headr2/  time
!
!
      praefix= '107'
      suffix2= '1a' ! '1f' ! '1e'
!
      praefixs = '/home/mtanaka/MPI_spin/SAI'//praefix 
      praefixi = '/home/mtanaka/MPI_spin/sai'//praefix 
      praefixc = '/home/mtanaka/MPI_spin/sai'//praefix 
!
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_comm_world,rank,ierror)
      call mpi_comm_size (mpi_comm_world,size,ierror)
!
      if(rank.eq.0) then
        open (unit=12,file=praefixc//'.12'//suffix2, &
                                  status='unknown',form='unformatted')
!
        write(06,*) 'Energy plot'
        write(06,*) 'file=',praefixc//'.12'//suffix2
!
        read(12) t8,xleng,yleng,zleng,rcut,rcutC,Temp,TCurie, &
                 tmax,dt,cptot
        read(12) x,y,z,vx,vy,vz,ch,mass,ag,ep
        read(12) x_0,y_0,z_0,rintC
        read(12) spx,spy,spz,sp2,spx00,spy00,spz00,r_ij
        read(12) aspx,aspy,aspz
        read(12) Jint,u,spin2,spin3,n_MCsteps,tau_diss
        read(12) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,toler,itermax
        read(12) dtwr,dtwr2,fw00,atsz,tsz0,av_tsz
        read(12) it,is,iwa,iwb,ir0
        read(12) np1,np2,spec,site,nintS,lintS,nintC,lintC,if_LJ
        read(12) i1,i2,i3,i4,disp_recv,cnt_recv,disp_recvC,cnt_recvC
        read(12) spinx,spinz,spin7,Bextx,Bextz,magx,magy,magz, &
                 Usys,conv,aitr,psdt,tfix,uss,usb,tsx,tsy,tsz,sum_mb, &
                 U_Fe,U_O,fdt4,vdt4,timeh
        read(12) e_sp0,e_c_r0,e_LJ0,if_LJ0
!
!% mpif90 -mcmodel=medium -fpic -O2 -o a.out @energy.f03 &> log
!% mpiexec -n 5 a.out &
!!      open (unit=33,file='fort.33',form='formatted')
!         do i=1,nhs
!         write(33,333) i,Usys(i),fdt4(i),timeh(i)
! 333     format(i6,1p3d10.2)
!         end do
!!        close (33)
!!      end if
!
        write(06,*) 't8=',t8
        write(06,*) 'np1,np2=',np1,np2
        write(06,*) 'if_LJ0=',if_LJ0
!
        close(12)
!
!             ++++++++++++++++++++++++++++
        open (unit=79,file=praefixc//'.79'//suffix2//'.ps', &
                                                 form='formatted')
        write(06,*) 'File=',praefixc//'.79'//suffix2//'.ps'
        time= t8  ! in /headr2/
!
        nframe= 4
        call gopen (nframe)
!
        call lplots 
        call gclose
!
        close (79)
      end if
!
      call mpi_finalize  (ierror)
!
      stop
      end program spin_dynamics 
!
!
!------------------------------------------------------------------------
      subroutine date_and_time7 (date_now,time_now)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer, dimension(8) :: ipresent_time
      character(len=10) :: date_now
      character(len=8)  :: time_now

      call date_and_time (values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)') &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end subroutine date_and_time7
!
!
!------------------------------------------------------
      subroutine lplots 
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit      none
      include      'param-spinRL7.h'
!
      real(C_float) spinx,spinz,spin7,Bextx,Bextz,magx,magy,magz, &
                    Usys,conv,aitr,psdt,tfix,uss,usb,tsx,tsy,tsz, &
                    sum_mb,U_Fe,U_O,ds_Fe,ds_O,fdt4,vdt4,timeh,   &
                    cosd,sind,B00,Bmw,mue_B,hh,av_mz(nhs),ss,s2,  &
                    emax,emax1,emax2,emax3,emax7,emin,emin1,      &
                    emin2,emin3,emin7,Temp3(nhs),www,dEE(nhs)
      common/ehist/ spinx(nhs),spinz(nhs),spin7(nhs), &
                    Bextx(nhs),Bextz(nhs),magx(nhs),magy(nhs),magz(nhs),&
                    Usys(nhs),conv(nhs),aitr(nhs),psdt(nhs),tfix(nhs),  &
                    uss(nhs),usb(nhs),tsx(nhs,3),tsy(nhs,3),tsz(nhs,3), &
                    sum_mb(nhs),U_Fe(nhs),U_O(nhs),ds_Fe(nhs),ds_O(nhs),&
                    fdt4(nhs),vdt4(nhs),timeh(nhs)
!
      real(C_DOUBLE) spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      real(C_DOUBLE) t8,xleng,yleng,zleng,pi,dt,tmax,cptot
      common/parm3/ t8,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
!
      character*8    label,cdate*10,ctime*8
      integer(C_INT) i,j,k,it,is,is0,ns,ILN,ILG,nxtick,nytick
      common/headr1/ label,cdate
      common/parm1/  it,is
!
      ILN= 1
      ILG= 2
      hh= 0.7
      call symbol (1.0,18.5,hh,'Spin dynamics (MD)',0.,18)
!
!-------------------------------------------------
!* One-period averaged mz(t)= mz - <mz>
!   40 data points are generated per period
!-------------------------------------------------
!
      is0= av_start*is  ! Start averaging at this step
!
      do 103 k= 1,is
      av_mz(k)= magz(k)
  103 continue
!
      do 100 k= 40,is-40
      ss= 0
      ns= 0
!
      do 130 i= k-39,k+40
      ss= ss +magz(i)
      ns= ns +1
  130 continue
!
      av_mz(k)= ss/ns
  100 continue
!
!
      ss= 0
      s2= 0
      ns= 0
!
      do 200 k= is0,is-40
      ss= ss +(magz(k) -av_mz(k))*Bextz(k)   ! magz= -g*<sz>
      s2= s2 +(magz(k) -av_mz(k))**2         ! <dm**2>
      ns= ns +1
  200 continue
!
      mue_B= 9.27410e-21   ! e*hbar/2mc
      B00= 100             ! gauss
      Bmw= B00*sqrt(Bapx**2 +Bapy**2 +Bapz**2)
!
!
      do i=1,nhs
      Temp3(i)= 0
      end do
!
      do i= 5,is-3
      www= 0
!
      do j= 1,5
      www= www +(Usys(i+3-j) -Usys(5))
      end do
!
      Temp3(i)= www/5.
      end do
!
      do i= 8,is-3
      dEE(i)= Temp3(i) -Temp3(i-10)
      end do
!
      nxtick= 3
      nytick= 3
!
!% mpif90 -mcmodel=medium -fpic -O2 -o a.out @energy3.f03 &> log
!% mpiexec -n 5 a.out &
      call lplmax (Usys,emax7,emin7,is)
      call lplot1 (2,4,is,timeh,Usys,emax7,emin7,ILN,nxtick,nytick,&
                   'Usys    ',8,'        ',8,'        ',8,0)
!
      call lplmax (U_Fe,emax1,emin1,is)
      call lplmax (U_O, emax2,emin2,is)
      emax= max(emax1,emax2)
      emin= min(emin1,emin2)
!
      call lplot1 (2,5,is,timeh,U_Fe,emax,0.,ILN,nxtick,nytick,&
                   'U_Fe    ',8,'        ',8,'        ',8,0)
      call lplot1 (2,6,is,timeh,U_o, emax,0.,ILN,nxtick,nytick,&
                   'U_O     ',8,' time   ',8,'        ',8,0)
!
      call lplmax (spinx,emax1,emin1,is)
      call lplmax (spinz,emax2,emin2,is)
      emax = amax1(emax1,emax2,-emin1,-emin2)
      emin = -emax
      call lplot1 (3,4,is,timeh,spinx,emax,emin,ILN,nxtick,nytick,&
                   'spin-Z. ',8,'        ',8,'        ',8,0)
      call lplot1 (3,5,is,timeh,spinz,emax,emin,ILN,nxtick,nytick,&
                   'spin-z.7',8,'        ',8,'        ',8,0)
!
      call lplot1 (3,6,is,timeh,uss,emax7,emin7,ILN,nxtick,nytick,&
                   'Uss   . ',8,' time   ',8,'        ',8,0)
!     call lplmax (spin7,emax,emin,is)
!     emax= 1.2*emax
!     call lplot1 (3,6,is,timeh,spin7,emax,0.,ILN,nxtick,nytick,&
!                  'spin-7  ',8,' time   ',8,'        ',8,0)
!------------------------
      call chart
!------------------------
!
!% mpif90 -mcmodel=medium -fpic -O2 -o a.out @energy3.f03 &> log
!% mpiexec -n 5 a.out &
      call lplmax (Temp3,emax1,emin1,is)
      call lplot1 (2,4,is,timeh,Temp3,emax1,0.,ILN,nxtick,nytick,&
                   'Temp    ',8,'        ',8,'        ',8,0)
!
      call lplmax (dEE,emax1,emin1,is)
      call lplot1 (2,5,is,timeh,dEE,emax1,-emax1,ILN,nxtick,nytick,&
                   'dT/dt   ',8,'        ',8,'        ',8,0)
!------------------------
      call chart
!------------------------
!
      if(.true.) return
!
      call lplmax (uss,emax1,emin1,is)
      call lplmax (usb,emax2,emin2,is)
      emax = amax1(emax1,emax2,-emin1,-emin2)
      emin = -emax
!
      call lplot1 (2,4,is,timeh,uss,emax7,emin7,ILN,nxtick,nytick,&
                   'Uss   . ',8,'        ',8,'        ',8,0)
      call lplot1 (2,5,is,timeh,usb,emax7,emin7,ILN,nxtick,nytick,&
                   'Usb     ',8,' time   ',8,'        ',8,0)
!------------------------
      call chart
!------------------------
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine lplot1 (ix,iy,npt1,x,y,ymax,ymin,IL,nxtick, &
                         nytick,lab1,n1,lab2,n2,lab3,n3,iskip)
!-----------------------------------------------------------------------
!  <<warning>>  order and number of arguments /lplot/ have been changed.
!               also, x (time) is defined for all range.
!               date: 5/18/96 at mit.
!***********************************************************************
!   il=1................ linear plot of (x,y)
!   il=2................ log10 plot of (x,log y)
!***********************************************************************
!
      use, intrinsic :: iso_c_binding
!
      integer(C_INT) ix,iy,npt1,IL,nxtick,nytick,n1,n2,n3,iplot
      real(C_float)  x(npt1),y(npt1),u(npt1),v(npt1),   &
                     xmax,xmin,ymax,ymin,time
      dimension  xcm(6),ycm(7),pl(6),pr(6),ql(7),qr(7)
!
      character*8    lab1,lab2,lab3,label,cdate*10
      common/headr1/ label,cdate
      common/headr2/ time
      common/pplcom/ nfine,pl1(10),pr1(10),ql1(10),qr1(10), &
                     xmin1(10),xmax1(10),ymin1(10),ymax1(10)
!
      data  xcm/21.0, 2*10.00, 3*6.00/,       &
            ycm/15.0, 2*6.80, 4*4.00/,        &
            pl/2.0,  2.0,14.0, 2.0,9.0,16.0/, &
            ql/2.3, 10.5,2.3, 14.0,9.5,5.0,0.5/
!
      call set_width (1.0)
!
      iplot=1
      go to 1
!
!-----------------------------------------------------------------------
      entry hplot1 (ix,iy,npt1,x,y,ymax,ymin,IL,nxtick,nytick,  &
                    lab1,n1,lab2,n2,lab3,n3,iskip)
!-----------------------------------------------------------------------
      iplot=2
!
    1 npt= npt1
      isc= 1
!
      do 5 i=1,6
    5 pr(i)= pl(i) +xcm(i)
!
      do 6 j=1,7
    6 qr(j)= ql(j) +ycm(j)
!
!                 ******************************************************
!*                **  Make a copy before the top-left frame is drawn. **
!                 ******************************************************
      hh= 0.70
      i1= iabs(ix)
      j1= iabs(iy)
      if(i1.ge.3) go to 10
      if(j1.eq.3.or.j1.ge.5) go to 10
!                                              ************************
!                                              ** label of the page. **
!                                              ************************
      call symbol (20.3,0.1,hh,'t=',0.,2)
      call values (21.5,0.1,hh,time,0.,101)
!
   10 continue
!
      do 23 i=1,npt
   23 u(i)= x(i)
      xmax= u(npt)
      xmin= u(1)
!                             ************************************
!                             ** three-point average if il > 0  **
!                             ************************************
      if(il.gt.0) then
        v(1)=   y(1)
        v(npt)= y(npt)
        v(npt-1)= y(npt-1)
!
        do 37 i=2,npt-2
   37   v(i)= 0.33333*(y(i-1)+y(i)+y(i+1))
      else
        do 38 i=1,npt
   38   v(i)= y(i)
      end if
!                                                *****************
!                                                **  log. scale **
!                                                *****************
      if(iabs(il).eq.2) then
         do 40 i=1,npt
         if(v(i).gt.0.) then
            v(i)= alog10(v(i))
         else
            v(i)= -10.
         end if
   40    continue
      end if
!                                **************************************
!                                ** Set a new scale and draw a frame.**
!                                **************************************
      if(iplot.eq.2) then
         ymax= -1.e10
         ymin=  1.e10
!
         do 50 i= 1,npt
         ymax= amax1(ymax,v(i))
         ymin= amin1(ymin,v(i))
   50    continue
!
         if(ymin.ge.0.) then
           ymax= 1.1*ymax
           ymin= 0.
         else
           ymax= amax1(0.,ymax)
           ymin= 1.1*ymin
         end if
      end if
!
      if(ymax.le.ymin) ymax= ymin+1.0
      if(iabs(il).eq.2) then
         if(ymax.gt.0.0) ymax= ymax+1.0
      end if
!
      dx= (xmax-xmin)/xcm(i1)
      dy= (ymax-ymin)/ycm(j1)
      x0= xmin
      y0= ymin
!
      call scalex (pl(i1),ql(j1),x0,y0,dx,dy,isc)
!
      pl1(isc)= pl(i1)
      pr1(isc)= pr(i1)
      ql1(isc)= ql(j1)
      qr1(isc)= qr(j1)
      xmin1(isc)= xmin
      xmax1(isc)= xmax
      ymax1(isc)= ymax
      ymin1(isc)= ymin
!                                                      *************
!                                                      **  Frame. **
!                                                      *************
      call plot (pl(i1),ql(j1),3)
      call plot (pl(i1),qr(j1),2)
      call plot (pr(i1),qr(j1),2)
      call plot (pr(i1),ql(j1),2)
      call plot (pl(i1),ql(j1),2)
!                                                    ******************
!                                                    **  Tick marks. **
!                                                    ******************
      scx= xcm(i1)/(nxtick+1)
      scy= ycm(j1)/(nytick+1)
!
      x0= pl(i1)
      y1= ql(j1)
      y4= qr(j1)
      y2= y1 +0.25
      y3= y4 -0.25
!
      do 62 k=1,nxtick
      x0= x0 +scx
      call plot (x0,y1,3)
      call plot (x0,y2,2)
      call plot (x0,y3,3)
      call plot (x0,y4,2)
   62 continue
!
      y0= ql(j1)
      x1= pl(i1)
      x4= pr(i1)
      x2= x1 +0.25
      x3= x4 -0.25
!
      do 63 k=1,nytick
      y0= y0 +scy
      call plot (x1,y0,3)
      call plot (x2,y0,2)
      call plot (x3,y0,3)
      call plot (x4,y0,2)
   63 continue
!                                                     **************
!                                                     ** Numbers. **
!                                                     **************
!
      hhs= 0.6
      call number (pl(i1)-0.8,ql(j1)-0.45,hhs,xmin,0.,101)
      call number (pr(i1)-1.1,ql(j1)-0.45,hhs,xmax,0.,101)
!
      call number (pl(i1)-1.8,ql(j1)     ,hhs,ymin,0.,101)
      call number (pl(i1)-1.8,qr(j1)-0.30,hhs,ymax,0.,101)
!
!                                                     **************
!                                                     **  Labels. **
!                                                     **************
      xc= 0.5*(pl(i1)+pr(i1))
      xu= xc -1.60
      xd= xc -0.20*n2/2
!
      yr= qr(j1)+0.15
      yl= ql(j1)-0.70
!
      call symbol (xu,yr,hh,lab1,0.,n1)
      call symbol (xd,yl,hh,lab2,0.,n2)
!
      xl= pl(i1)-1.50
      yc= 0.5*(ql(j1)+qr(j1))
      call symbol (xl,yc,hh,lab3,0.,n3)
!                                     **********************************
!                                     **  No plot is made if npt1 < 0 **
!                                     **********************************
   70 if(npt1.lt.0) return
!
      isk= 0
      if(iskip.ne.0) then
        call set_width (2.0)
      end if
!
      if(iskip.eq.1) then
        itot= 5
        ibr = 3 ! 2
      end if
      if(iskip.eq.2) then
        itot= 13
        ibr = 9
      end if
      if(iskip.eq.7) then
        itot= 3
        ibr = 1
      end if
!
      call plotl (u(1),v(1),isc,3)
!**
      if(iplot.eq.1) then
         do 100 i=1,npt
         isk= isk +1
!
         if(iskip.eq.0) then
           call plotl (u(i),v(i),isc,2)
         else
           if(mod(isk,itot).lt.ibr) then
             call plotl (u(i),v(i),isc,2)
           else
             call plotl (u(i),v(i),isc,3)
           end if
         end if
  100    continue
      else
         do 120 i=1,npt-1
         call plotl (u(i+1),v(i)  ,isc,2)
         call plotl (u(i+1),v(i+1),isc,2)
  120    continue
      end if
!**
      call plotl (u(npt),v(npt),isc,3)
      call set_width (1.0)
!
      return
      end subroutine lplot1
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      include    'param-spinRL7.h'
!
      integer(C_INT) is
      real(C_float) f(is),fmax,fmin
!
      fmax= -1.e10
      fmin=  1.e10
!
      do 100 i= 1,is
      fmax= amax1(fmax,f(i))
      fmin= amin1(fmin,f(i))
  100 continue
!
      return
      end subroutine lplmax 
!
!
!------------------------------------
      subroutine set_width (awid)
!------------------------------------
      use, intrinsic :: iso_c_binding
      real(C_float) amid
!
      write(79,*) 'stroke'
      write(79,10) awid
   10 format(f4.1,' setlinewidth')
!
      return
      end subroutine set_width 
!
!
!-----------------------------------------------------------------------
      subroutine cplot3 (q,xmax8,ymax8,zmax8,char,nc)
!-----------------------------------------------------------------------
!***********************************************************************
!*   contour plots of scalar quantities.                               *
!***********************************************************************
      use, intrinsic :: iso_c_binding
!
      include    'param-spinRL7.h'
      parameter  (mx1=meshx+1,my1=meshy+1,mz1=meshz+1)
      parameter  (mx=meshx,my=meshy,mz=meshz)
!
      character*8  char,label,cdate*10
      common/headr1/ label,cdate
      common/headr2/ time
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(mx1),pxc(mx1),pxl(mx1),pyr(my1),pyc(my1), &
                     pyl(my1),pzr(mz1),pzc(mz1),pzl(mz1)
      real(C_DOUBLE) q(mx,my,mz),xmax8,ymax8,zmax8
      real(C_float)  a(2048),b(2048),ww(2048),cut(200,4)
!
      j0= my/2 +1
      k0= mz/2 +1
      xmax= xmax8
      ymax= ymax8
      zmax= zmax8
!
!* 1. Plot at k= k0: subscript j first.
!
      npx= 0
      ij= 0
      qc = 1./16.
!***
      do 10 i= 1,mx
      npx= npx +1
      ir= pxr(i)
      il= pxl(i)
!
      npy= 0
      do 10 j= 1,my
      npy= npy +1
      jr= pyr(j)
      jl= pyl(j)
!
      ij= ij+1
      a(ij)= qc*(   q(ir,jr,k0) +2.*q(ir,j,k0)    +q(ir,jl,k0)  &
                +2.*q(i ,jr,k0) +4.*q(i ,j,k0) +2.*q(i ,jl,k0)  &
                +   q(il,jr,k0) +2.*q(il,j,k0)    +q(il,jl,k0) )
   10 continue
!
!* 2. Plot at j= j0: subscript k first.
!
      npx= 0
      ij= 0
      qc = 1./16.
!***
      do 20 i= 1,mx
      npx= npx +1
      ir= pxr(i)
      il= pxl(i)
!
      npz= 0
      do 20 k= 1,mz
      npz= npz +1
      kr= pzr(k)
      kl= pzl(k)
!
      ij= ij+1
      b(ij)= qc*(   q(ir,j0,kr) +2.*q(ir,j0,k)    +q(ir,j0,kl)  &
                +2.*q(i ,j0,kr) +4.*q(i ,j0,k) +2.*q(i ,j0,kl)  &
                +   q(il,j0,kr) +2.*q(il,j0,k)    +q(il,j0,kl) )
   20 continue
!
!
      hh = 0.70
      call symbol (0.1,18.2,hh,label,0.,8)
      call symbol (13.0,0.7,hh,cdate,0.,10)
      call symbol (13.0,0.1,hh,'t =',0.,3)
      call values (999.0,999.0,hh,time,0.,101)
!
      xl1=  1.8
      xr1=  9.3
      xl2= 10.0
      xr2= 17.5
!
      zl=  1.0
      zr=  zl +(xr1 -xl1)*zmax/xmax
      if(zr.gt.10.) zr= 10.
!                          <--- limit elongated y-length.
!
      yl=  zr +1.
      yr=  yl +(xr1 -xl1)*ymax/xmax
      if(yr.gt.25.) yr= 25.
!                          <--- limit elongated y-length.
!
      xc1= 0.5*(xr1+xl1)
      xc2= 0.5*(xr2+xl2)
      yc=  0.5*(yr+yl)
      zc=  0.5*(zr+zl)
!
!---------------------------------------------
!*  **Maximum of the vectors**
!---------------------------------------------
!
      am2= 0.
      am4= 0.
!
      do 100 ij= 1,npx*npy
      am2= amax1(am2,abs(a(ij)))
  100 continue
!
      do 200 ij= 1,npx*npz
      am4= amax1(am4,abs(b(ij)))
  200 continue
!
      ams= amax1(am2,am4)
      if(ams.lt.1.e-10) ams=999.0
!
      call symbol (zl,0.10,hh,'scalar.max=',0.,11)
      call values (999.0,999.0,hh,ams,0.,101)
!
!---------------------------------------------
!*  (1): Contours in (x,z) plane.
!---------------------------------------------
!
      call setscl (0.,0.,zmax,xmax,zl,xl2,zr,xr2,gdz,gdx, &
                   nc,char,6,' (x-z)',0.4, &
                   1,'z',0.4,1,'x',0.4,1)
!
      call values (zl-0.45,xl2-0.5,hh,0.0,0.,101)
      call values (zl-1.3,xr2-0.3, hh,xmax,0.,101)
      call values (zr-1.3,xl2-0.5, hh,zmax,0.,101)
!
      nxz= npx*npz
      call daisho (b,nxz,wamin,wamax)
!
      ncontr= 11
      call eqcntr (b,ww,npz,npx,zl,xl2,zr,xr2,wamin,0.0,wamax, &
                   ncontr,1)
!
!---------------------------------------------
!*  (2): Contours in (x,y) plane.
!---------------------------------------------
!
      call setscl (0.,0.,ymax,xmax,yl,xl2,yr,xr2,gdy,gdx, &
                   1,' ',1,' ',0.4, &  
                   1,'y',0.4,1,'x',0.4,1)
!
      call values (yl-0.45,xl2-0.5,hh,0.0,0.,101)
      call values (yl-1.3,xr2-0.3, hh,xmax,0.,101)
      call values (yr-0.3,xl2-0.5, hh,ymax,0.,101)
!
      nxy= npx*npy
      call daisho (a,nxy,wamin,wamax)
!
      ncontr= 11
      call eqcntr (a,ww,npy,npx,yl,xl2,yr,xr2,wamin,0.0,wamax, &
                   ncontr,1)
!
!---------------------------------------------
!*  (3): Cut plots.
!---------------------------------------------
!
      do 300 jj= 1,4
      j= (my/4)*(jj-1) +1
!
      do 300 i= 1,mx
      cut(i,jj)= q(i,j,k0)
  300 continue
!
!
      amax7= -1.e+10
      amin7=  1.e+10
!
      do 320 jj= 1,4
      do 320 i= 1,mx
      amax7= amax1(cut(i,jj),amax7)
      amin7= amin1(cut(i,jj),amin7)
  320 continue
!
      if(amax7.lt.0.) amax7= 0.
      if(amin7.gt.0.) amin7= 0.
!
!
      dd= amax7 -amin7
      dx= (yr -yl)/6.
      dy= xr1 -xl1
!
      do 340 jj= 1,4
      xo= yl +1.5*(jj-1)*dx
      xu= xo +dx
      call plot (xo,xl1,3)
      call plot (xu,xl1,2)
      call plot (xu,xr1,2)
      call plot (xo,xr1,2)
      call plot (xo,xl1,2)
!
!* zero line.
      x1= xo +dx*(0. -amin7)/dd
      call plot (x1,xl1,3)
      call plot (x1,xr1,2)
!
      x1= xo +dx*(cut(1,jj) -amin7)/dd
      y1= xl1 
      call plot (x1,y1,3)
!
      do 340 i= 1,mx
      x1= xo  +dx*(cut(i,jj) -amin7)/dd
      y1= xl1 +dy*(i-1)/float(mx-1)
!
      call plot (x1,y1,2)
  340 continue
!---------------------
      call chart
!---------------------
      return
      end subroutine cplot3 
!
!
!-----------------------------------------------------------------------
      subroutine eqcntr(u,w,nx,ny,xl,yl,xr,yr,umin,ubund,umax, &
                        lank,iwaku)
!-----------------------------------------------------------------------
!  << eqcntr >>
!         presented by kunihiko.watanabe 14.nov.1989
!         reviced   by hisanori.takamaru 16.mar.1990
!-----------------------------------------------------------------------
!     1. function
!        (1) to draw tokosen
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) u       nx,ny     (i)       world value
!        (2) w       nx,ny     (i)       work array (real*4)
!        (3) xl,xr,yl,yr       (i)       absolute coordinate value
!        (4) umin,umax         (i)       height of max & min
!                                        umin>umax : automatic control
!        (5) ubund             (i)       draw dash line (u < ubund)
!        (6) lank              (i)       number of draw lines
!        (7) iwaku             (i)       =1 : draw frame
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      dimension u(1),w(1)
!
      if (nx.lt.2) return
      if (ny.lt.2) return
      if (xr.lt.xl) return
      if (yr.lt.yl) return
!
      nxy = nx*ny
      nxm1 = nx - 1
      nym1 = ny - 1
!
      dx = (xr-xl)/ float(nxm1)
      dy = (yr-yl)/ float(nym1)
!
      umax1 = umax
      umin1 = umin
!
      if(umax1.gt.(1.000001*umin1)) then
!
        do 10 i = 1 , nxy
          w(i) = u(i) - umin1
          if(u(i).gt.umax1) w(i) = umax1 - umin1
          if(u(i).lt.umin1) w(i) = 0.
   10   continue
!
      else
!
        umax1=-1.e+30
        umin1= 1.e+30
        do 20 i = 1 , nxy
          umax1=amax1(umax1,u(i))
          umin1=amin1(umin1,u(i))
   20   continue
        do 25 i = 1 , nxy
          w(i) = u(i) - umin1
   25   continue
!
      endif
!
!------------------------------------------------
      if(umax1.le.(1.000001*umin1))  return
!------------------------------------------------
!
      if(iwaku.eq.1) then
        call plot(xl,yl,3)
        call plot(xr,yl,2)
        call plot(xr,yr,2)
        call plot(xl,yr,2)
        call plot(xl,yl,2)
        call plot(xl,yl,3)
      endif
!
      uld = float(lank+1) / (umax1-umin1)
      eps = 1.0e-8
!
      nxym1 = nxm1*nym1
      do 9000  ijnxy1 = 1,nxym1
        j = (ijnxy1-1)/nxm1 + 1
        i = ijnxy1 - (j-1)*nxm1
!
          i1 = i + nx * (j - 1)
          i2 = i1 + 1
          i3 = i1 + 1 + nx
          i4 = i1 + nx
!
          u1 =  w(i1) * uld
          u2 =  w(i2) * uld
          u3 =  w(i3) * uld
          u4 =  w(i4) * uld
!
          k1 = ifix(u1)
          k2 = ifix(u2)
          k3 = ifix(u3)
          k4 = ifix(u4)
!
          j1 = iabs(k2-k1)
          j2 = iabs(k3-k2)
          j3 = iabs(k4-k3)
!
          if(j1.ne.0) then
            do 1000 ll = 1 , j1
              u0 = float(ll) + float(min0(k1,k2))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
!
              if(abs(u2-u1).lt.eps)                 go to 1000
!
              x1 = xl + dx * ( (u0-u1)/(u2-u1) + float(i-1) )
              y1 = yl + dy * float(j-1)
!
              if( ((u3-u0)*(u2-u0)).gt.0. )         go to 1100
              if( ( (u0-u2).gt.0. ).and.( (u0-u4).gt.0. ) ) go to 1100
              if( abs(u3-u2).lt.eps )               go to 1100
!
                x2 = xl + dx * float(i)
                y2 = yl + dy * ( (u0-u2)/(u3-u2) + float(j-1) )
!
                call wdash(x1,y1,x2,y2,jouge)
!
 1100         continue
              if( ((u4-u0)*(u3-u0)).gt.0. )         go to 1200
              if( ((u1-u0)*(u3-u0)).gt.0. )         go to 1200
              if( ((u2-u0)*(u4-u0)).gt.0. )         go to 1200
              if( abs(u4-u3).lt.eps )               go to 1200
!
                x2 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
                y2 = yl + dy * float(j)
                call wdash(x1,y1,x2,y2,jouge)
!
 1200         continue
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 1300
              if( ( (u0-u1).gt.0. ).and.( (u0-u3).gt.0. ) ) go to 1300
              if( abs(u1-u4).lt.eps )               go to 1300
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy*((u0-u1)/(u4-u1)+float(j-1))
                call wdash(x1,y1,x2,y2,jouge)
 1300         continue
 1000       continue
!
          endif
!
          if(j2.ne.0) then
!
            do 2000 ll = 1 , j2
              u0 = float(ll) + float(min0(k2,k3))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
              if( abs(u3-u2).lt.eps )               go to 2000
!
              x1 = xl + dx * float(i)
              y1 = yl + dy * ( (u0-u2)/(u3-u2) + float(j-1) )
!
              if( ((u4-u0)*(u3-u0)).gt.0. )         go to 2100
              if( ( (u0-u1).gt.0. ).and.( (u0-u3).gt.0. ) ) go to 2100
              if( abs(u4-u3).lt.eps )               go to 2100
!
                x2 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
                y2 = yl + dy * float(j)
!
                call wdash(x1,y1,x2,y2,jouge)
!
 2100         continue
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 2200
              if( ((u1-u0)*(u3-u0)).gt.0. )         go to 2200
              if( ((u2-u0)*(u4-u0)).gt.0. )         go to 2200
              if( abs(u1-u4).lt.eps )               go to 2200
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy * ( (u0-u1)/(u4-u1)+float(j-1) )
                call wdash(x1,y1,x2,y2,jouge)
 2200         continue
 2000       continue
!
          endif
!
          if(j3.ne.0) then
!
            do 3000 ll = 1 , j3
              u0 = float(ll) + float(min0(k3,k4))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
              if( abs(u4-u3).lt.eps )               go to 3000
!
              x1 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
              y1 = yl + dy * float(j)
!
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 3100
              if( ( (u0-u2).gt.0. ).and.( (u0-u4).gt.0. ) ) go to 3100
              if( abs(u1-u4).lt.eps )               go to 3100
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy * ( (u0-u1)/(u4-u1) + float(j-1) )
                call wdash(x1,y1,x2,y2,jouge)
 3100         continue
 3000       continue
          endif
 9000 continue
!
      return
      end subroutine eqcntr
!
!
!-----------------------------------------------------------------------
      subroutine setscl (wminx,wminy,wmaxx,wmaxy, xl,yl,xr,yr,gdx,gdy, &
                 n1,char1,n2,char2, hight1,                            &
                 nnx,charx,hightx, nny,chary,highty, iwaku)
!-----------------------------------------------------------------------
!  << setscl >>                   /char1/
!                          wmaxy  +--------------------+  (xl,yl)
!                               y |             (xr,yr)|  (xr,yr) on 0
!                               r |                    |
!                               a |                    |
!    (wminx,wminy)              h |                    |
!    (wmaxx,wmaxy) on is        c |                    |
!                                 |(xl,yl)             |
!                          wminy  +--------+--+--------+
!                                 wminx  /charx/       wmaxx
!-----------------------------------------------------------------------
!
!     setscl
!
!     1. function
!        (1) to scale the graphics by calcomp specifications
!     2. arguments            (i/o)     (meaning)
!        (1) wminx,wmaxx,
!            wminy,wmaxy       (i)       world coordinate value
!        (2) xl,xr,yl,yr       (i)       absolute coordinate value
!        (3) gdx,gdy           (o)       scaling factor of coordinate
!                                        from world to absolute
!        (4) char1,charx,cahry (i)       title on graph,x-axis,y-axis
!        (5) iwaku             (i)       draw frame (0:off ; 1:on)
!                                         999 : write out only title,
!                                                    not draw otherwise
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!             (** symbol **)
!             (** number **)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      character*1  char1(1),char2(1),charx(1),chary(1)
!
      if (wmaxx.le.wminy) goto 9999
      if (wmaxx.le.wminy) goto 9999
      if (xr.le.xl)       goto 9999
      if (yr.le.yl)       goto 9999
!
      gdx= (xr-xl)/(wmaxx-wminx)
      gdy= (yr-yl)/(wmaxy-wminy)
!
      xc = 0.5*( xr + xl )
      yc = 0.5*( yr + yl )
!
      if (n1 .gt.0) then
        if (hight1.gt.0) then
          xs1= xc -0.5*n1*hight1
          xs2= xs1 +(n1+1)*hight1
          call symbol(xs1,yr+0.1,hight1,char1(1),0.,n1)
          call symbol(xs2,yr+0.1,hight1,char2(1),0.,n2)
        end if
      end if
!-----------------------------------------------------------------------
      if (iwaku.eq.999) return
!-----------------------------------------------------------------------
!
      if (iwaku.eq.1) then
        call plot (xl,yl,3)
        call plot (xl,yr,2)
        call plot (xr,yr,2)
        call plot (xr,yl,2)
        call plot (xl,yl,2)
        call plot (999.,999.0,3)
      end if
!
      if (nnx.gt.0) then
        if (hightx.gt.0) then
          call symbol(xc-0.5*hightx*nnx,yl-0.5,hightx,charx(1),0.,1)
          do 200 nnx1=2,nnx
  200     call symbol(999.0,999.0,hightx,charx(nnx1),0.,1)
        end if
      end if
      if (nny.gt.0) then
        if (highty.gt.0) then
          call symbol(xl-0.5,yc-0.5*highty*nny,highty,chary(1),0.,1)
          do 300 nny1=2,nny
  300     call symbol(999.0,999.0,highty,chary(nny1),0.,1)
        end if
      else if(nny.lt.0) then
        if (highty.gt.0) then
          call symbol(xc-0.5*highty*nny,yc,highty,chary(1),0.,1)
          do 400 nny1=2,nny
  400     call symbol(999.0,999.0,highty,chary(nny1),0.,1)
        end if
      end if
!
      return
!
!-----------------------------------------------------------------------
!
 9999 continue
      write(6,*) '**********  abnormal world coordinate ********'
      write(6,*) '      '
      write(6,*) '    wmaxx =',wmaxx,' wminx = ',wminx
      write(6,*) '    wmaxy =',wmaxy,' wminy = ',wminy
      write(6,*) '    xl,yl,xr,yr =',xl,yl,xr,yr
      write(6,*) '    fctr  =',fctr
      write(6,*) '      '
      call chart
      call symbol(1.0,10.0,0.2,' abnormal world coordinate call',0.,31)
      call symbol(1.0,09.0,0.2,' wmaxx =',0.,8)
      call number(999.0,999.0,0.2,wmaxx,0.,2)
      call symbol(1.0,08.5,0.2,' wminx =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,08.0,0.2,' wmaxy =',0.,8)
      call number(999.0,999.0,0.2,wmaxy,0.,2)
      call symbol(1.0,07.5,0.2,' wminy =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,07.0,0.2,' fctr  =',0.,8)
      call number(999.0,999.0,0.2,fctr,0.,2)
      call symbol(1.0,06.5,0.2,' xleft =',0.,8)
      call number(999.0,999.0,0.2,xl,0.,2)
      call symbol(1.0,06.0,0.2,' yleft =',0.,8)
      call number(999.0,999.0,0.2,yl,0.,2)
      call symbol(1.0,05.5,0.2,' xright=',0.,8)
      call number(999.0,999.0,0.2,xr,0.,2)
      call symbol(1.0,05.0,0.2,' yright=',0.,8)
      call number(999.0,999.0,0.2,yr,0.,2)
      call exit
      return
      end subroutine setscl 
!
!
!-----------------------------------------------------------------------
      subroutine scalex (xcm,ycm,x00,y00,dx,dy,isc)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      x0(isc)= x00
      y0(isc)= y00
      dxi(isc)= 1./dx
      dyi(isc)= 1./dy
!
      xl(isc)= xcm
      yl(isc)= ycm
!
      return
      end subroutine scalex 
!
!
!-----------------------------------------------------------------------
      subroutine plotl (x,y,isc,ipl)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      xcm= xl(isc) +dxi(isc)*(x -x0(isc))
      ycm= yl(isc) +dyi(isc)*(y -y0(isc))
!
      call plot (xcm,ycm,ipl)
!
      return
      end subroutine plotl 
!
!
!-----------------------------------------------------------------------
      subroutine values (x,y,height,val,theta,ifmat)
!-----------------------------------------------------------------------
!  << values >>
!     1. function
!        (1) to draw variable
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) x,y               (i)       absolute coordinate value
!        (2) height            (i)       draw out size on paper
!        (3) val               (i)       variable
!        (4) theta             (i)       angle
!        (5) ifmat             (i)       format type
!     3. called by
!             (** nothing **)
!     4. calls
!             (** number **)
!             (** symbol **)
!-----------------------------------------------------------------------
!        ifmat = (n100)*100 + keta
!        n100 = 0 : integer format
!        n100 = 1 : f format ::  number(x,y,height,val,theta,keta)
!        n100 = 2 : e format ::
!        n100 = 3 : power of ten format
!        n100 = othewise : not write out
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      real*4 val
      character chr13*13,chr12*12,chr3*3
      character*1 minus,zero,blank
      parameter(ratio = 6./7. )
      data minus/'-'/,zero/'0'/,blank/' '/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ifmat.lt.0) return
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n100 = ifmat/100
      keta = ifmat - n100*100
!
      if (n100.eq.0) then
        call number(x,y,height,val,theta,-1)
      else if (n100.eq.1) then
        call number(x,y,height,val,theta,keta)
      else if (n100.eq.2) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:4) = chr13(1:3)//'e'
          numsym = 4
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else if (val.eq.0) then
            chrval = val
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+2) = chr13(2:keta+2)//'e'
            numsym = keta + 2
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
        call symbol(999.,999.,height,chr3,theta,numsy1)
      else if (n100.eq.3) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:6) = chr13(1:3)//'x10'
          numsym = 6
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+5) = chr13(1:keta+2)//'x10'
            numsym = keta + 5
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+4) = chr13(2:keta+2)//'x10'
            numsym = keta + 4
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        sint = sin(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
!
!                                             *******************
!                                             ** exponent part **
!                                             *******************
!
        h2 = height * 5./7.
        x1 = (numsym+1)* height * ratio
        y1 = height * 4./7.
        if (abs(theta).lt.1e-04) then
          x1 = x + x1
          y1 = y + y1
        else
          x2 =     x1 * cost - y1 * sint
          y1 = y + x1 * sint + y1 * cost + h2*cost
          x1 = x + x2                    - h2*sint
        end if
        call symbol(x1,y1,h2,chr3,theta,numsy1)
      end if
      return
      end subroutine values 
!
!
!-----------------------------------------------------------------------
      subroutine wdash (x1,y1,x2,y2,ipen )
!-----------------------------------------------------------------------
!  << wdash  >>                      ver 2.00   16.mar.1990
!
!     1. function
!        (1) to draw line from (x1,y1) to (x2,y2) by wdash
!                            in absolute coordinate
!     2. arguments            (i/o)     (meaning)
!        (1) x1,x2,y1,y2       (i)       absolute coordinate value
!        (2) ipen              (i)       pen type of 'wdash'
!     3. called by
!             (** eqcntr  **)
!             (** wdashl  **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
!       ipen : meaning           - : 0.05 (cm)
!        1   :       line     -------------------
!        2   :  dash line     --- --- --- --- ---
!        3   :  dash line     -- -- -- -- -- -- --
!        4   :  dash line     - - - - - - - - - -
!        5   :  1 point dash  ---- - ---- - ---- -
!        6   :  2 point dash  --2.0-- - - --2.0--
!   otherwise:  line          ---------------------
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      h1  =  0.05
      h2  =  2.0 * h1
      h3  =  3.0 * h1
      h4  =  4.0 * h1
      h20 = 20.0 * h1
      call plot ( x1 , y1 , 3 )
      k = - 1
      if(ipen.lt.2) then
        go to 999
      else if(ipen.eq.2) then
        hh1 = h3
        hh2 = h1
      else if (ipen.eq.3) then
        hh1 = h2
        hh2 = h1
      else if (ipen.eq.4) then
        hh1 = h1
        hh2 = h1
      else if (ipen.eq.5) then
        hh1 = h4
        hh2 = h1
        hh3 = h1
        hh4 = h1
      else if (ipen.eq.6) then
        hh1 = h20
        hh2 = h1
        hh3 = h1
        hh4 = h1
        hh5 = h1
        hh6 = h1
      end if
      if(ipen.lt.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hhh
  200   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hhh
          d=d+hh1
          goto 200
        end if
      else if (ipen.eq.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hhh
  500   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hhh
          d=d+hh1
          goto 500
        end if
      else if (ipen.eq.6) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hh5
        hh5 = hh6
        hh6 = hhh
  600   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hh5
          hh5 = hh6
          hh6 = hhh
          d=d+hh1
          goto 600
        end if
      end if
  999 call plot ( x2 , y2 , ( 5 + k ) / 2 )
      call plot ( x2 , y2 , 3)
      return
      end subroutine wdash 
!
!
!-----------------------------------------------------------------------
      subroutine daisho(x,nx,xmin1,xmax1)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      dimension x(1)
!
      xmax1= x(1)
      xmin1= x(1)
      do 100 i=2,nx
      xmax1= amax1(xmax1,x(i) )
      xmin1= amin1(xmin1,x(i) )
  100 continue
!
      return
      end subroutine daisho
!
!
!***************************************************************
!*   This program package generates a unix postscript          *
!*   graphic file when called by calcomp-compatible /plot23.f/ *.  
!***************************************************************
!----------------------------------------------------------
!    Postscript header by fortran
!        t. ogino (nagoya university) February 27, 1992
!      modified to conform gsipp commands
!        Motohiko Tanaka (nifs)       November 23, 1993
!
!----------------------------------------------- 5/27/96 -------
!   This ps-adobe-2.0 header allows us full paging features in
!   the ghostview.  to scroll up the page (backward), click the 
!   page number and press two buttons of mouse simultaneously.
!
!   consult: A.Saitou (kyoto u.)  the definition of /@eop  
!   needs stroke for line drawings (not in the tex header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  this is an adobe-2.0 postscript file.
!
       write(79,10)
   10  format('%!ps-adobe-2.0',/      &
              '%%pages: (atend)',/    &
              '%%pageorder: ascend',/ &
              '%%endcomments',/       &
              '%%begindocument')
!
!%%%%%%%%%%%%%%%%%%% procedure defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(79,11) 
!  11 format('%%boundingbox: 150. 400. 550. 600.')
!
      write(79,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(79,23) 
   23 format('/tr {/times-roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/  &
             '/se {setfont} bind def',/    &
             '/ro {rotate}  bind def',/    &
             '/tl {translate} bind def',/  &
             '/sc {scale} bind def')
!
      write(79,24) 
   24 format('/@bop          % @bop -- begin the a new page',/ &
             '{erasepage newpath initgraphics',/ &
             '/saveimage save def',/ &
             '} bind def')
!
      write(79,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/   &
             ' saveimage restore',/ &
             '} bind def')
!
      write(79,26) 
   26 format('/@end          % @end -- done the whole shebang',/ &
             ' /end load def')
!
      write(79,27) 
   27 format('/dir 0 def')
!
      write(79,29) 
   29 format('/s             % string s -- show the string',/ &
             '{dir 1 eq',/                                    &
             ' {gsave currentpoint translate 90 rotate 0 0 moveto',/ &
             ' show grestore}',/ &
             ' {show} ifelse',/  &
             '} bind def')
!
      write(79,31)
   31 format('%%enddocument',/ &
             '%%endprolog',/   &
             '%%beginsetup',/  &
             '/resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%endsetup')
!
!%%%%%%%%%%%%%%%%%%% end of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(79,12) ipage,ipage
   12  format('%%page:',1x,i2,1x,i2)
!
       write(79,30) 
   30  format('%%beginpagesetup',/ &
              '%%endpagesetup',/   &
              '@bop')
!
!
!*  Set magnifying factor (gsipp to sun coordinate).
!   rotate and translate to output on a4-l paper.
!      left corner ...... (  0.,  0.)
!      right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(79,*) '90.0 ro'
       write(79,*) '50.0 -550.0 tl'
!
!*  if nfrm=4, four frames in a page (top-left frame).
!
       if(nfrm.eq.1) then
          write(79,*) '1.00 1.00 sc'
       else
          write(79,*) '0.50 0.50 sc'
          write(79,*) '0.0 550.0 tl'
       end if
!
       return
       end subroutine gopen 
!
!
!-----------------------------
       subroutine gclose
!-----------------------------
       call plote
       return
       end
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(79,10) 
   10  format('@eop')
       return
       end
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*     four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(79,10) 
   10     format('%%pagetrailer    % need for the page count')
!
          write(79,20) lpage,lpage
   20     format('%%page:',1x,i2,1x,i2)
!
          write(79,30) 
   30     format('%%beginpagesetup',/ &
                 '%%endpagesetup',/   &
                 '@bop')
!
          write(79,*) '90.0 ro'
          write(79,*) '50.0 -550.0 tl'
!
          if(nfrm.eq.1) then
             write(79,*) '1.00 1.00 sc'
          else
             write(79,*) '0.50 0.50 sc'
             write(79,*) '0.0  550.0 tl'
          end if
!
          return
       end if
!
!-----------------------------------------------------
!   First cancel the previous translation, then
!   make a new translation (scale factor alive).
!-----------------------------------------------------
!*   frames 2-4:
!
       if(loc.eq.1) then
          write(79,*) '  0.0 -550.0 tl'
          write(79,*) '700.0  550.0 tl'
       end if
!
       if(loc.eq.2) then
          write(79,*) '-700.0 -550.0 tl'
          write(79,*) '   0.0    0.0 tl'
       end if
!
       if(loc.eq.3) then
          write(79,*) '  0.0 0.0 tl'
          write(79,*) '700.0 0.0 tl'
       end if
!
       return
       end subroutine chart
!
!
!------------------------------------
       subroutine factor(fct)
!------------------------------------
       write(79,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end
!
!
!---------------------------------------
       subroutine newcolor (ic,r,g,b)
!---------------------------------------
!  ic= 3 tri-color
!  ic= 0 gray scale, r= 0. for black
!
       write(79,*) 'stroke'
!
       if(ic.eq.0) then
         write(79,10) 1.-r  ! 0. for black
   10    format(f4.1,' setgray')
       end if
!
       if(ic.eq.3) then
         write(79,30) r,g,b
   30    format(3f4.1,' setrgbcolor')
       end if
!
       return
       end subroutine newcolor 
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
       write(79,*) 'sn'
       pi1=0.40*float(i1-1)
       write(79,30) pi1
   30  format(f3.1,' sl')
       if(i2.ne.1) then
       write(79,*) '[2 2] 0 sd'
       endif
!
       return
       end subroutine newpen 
!
!
!-----------------------------
       subroutine linee
!-----------------------------
       write(79,*) 'st'
       return
       end
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
!
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
!
       if(ip.eq.3)  write(79,10) x,y
       if(ip.eq.2)  write(79,20) x,y
       if(ip.eq.-3) write(79,30) x,y
       if(ip.eq.-2) write(79,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
!       write(79,*) 'st'
!
       return
       end subroutine plot 
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       character isymb*80,ica*80,ich(80)*1
       equivalence (ica,ich(1))
!
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
!
       write(79,*) 'tr'
       write(79,10) h
   10  format(f5.1,' sf')
       write(79,*) 'se'
       write(79,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(79,30) ang
   30  format(f5.1,' ro')
!*
       ica= isymb
       write(79,*) '(',(ich(i),i=1,n),') s'
!
       return
       end subroutine symbol 
!
!
!-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       character  isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(79,*) 'tr'
       write(79,10) h
   10  format(f5.1,' sf')
       write(79,*) 'se'
!
       write(79,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(79,30) ang
   30  format(f5.1,' ro')
!
       write(isymb,40) anu
   40  format(1pe9.2)
       write(79,*) '(',isymb,') s'
!
       return 
       end subroutine number 
!
!
!---------------------------------------------------
       subroutine sunscl (x,y,h,n)
!---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
!
       if(x.eq.999.) then
         x= x0 +iabs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
!
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
!
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
!
       return
       end subroutine sunscl 
