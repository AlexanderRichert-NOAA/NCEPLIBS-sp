c Test sptruns and sptrunsv
c Interpolate heights and winds from a latlon grid to two antipodal polar stereographic grids.
c Subprograms GETGB and PUTGB from w3emc are referenced.
c  unit number 11 is the input latlon grib file
c  unit number 31 is the input latlon grib index file
c  unit number 51 is the output northern polar stereographic grib file
c  unit number 52 is the output southern polar stereographic grib file
c  nominal spectral truncation is r40
c  maximum input gridsize is 360x181
c  maximum number of levels wanted is 12
      parameter(lug=11,lui=0,lun=51,lus=52)
      parameter(iromb=1,maxwv=40,jf=360*181,kx=12)
      integer kp5(kx),kp6(kx),kp7(kx)
      integer kpo(kx)
      data kpo/1000,850,700,500,400,300,250,200,150,100,70,50/

      CALL BAOPENR(lug,'ref_gdaswave.t00z.wcoast.0p16.f000.grib1',IRET)
c      open(lui,file='tmidaily20112012.grb.idx')
c height
      km=12
      kp5=7
      kp6=100
      kp7=kpo
      call gs65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c winds
      km=12
      kp5=33
      kp6=100
      kp7=kpo
      call gv65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c
      stop 0
      end
c
      subroutine gs65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c  interpolates a scalar field using spectral transforms.
      integer kp5(km),kp6(km),kp7(km)
c  output grids are 65x65 (381 km true at latitide 60).
c  nh grid oriented at 280E; sh grid oriented at 100E.
      parameter(nph=32,nps=2*nph+1,npq=nps*nps)
      parameter(true=60.,xmesh=381.e3,orient=280.)
      parameter(rerth=6.3712e6)
      parameter(pi=3.14159265358979,dpr=180./pi) 
      real gn(npq,km),gs(npq,km)
      integer jpds(25),jgds(22),kpds(25,km),kgds(22,km)
      logical lb(jf)
      real f(jf,km)
c
      g2=((1.+sin(abs(true)/dpr))*rerth/xmesh)**2
      r2=2*nph**2
      rlatn1=dpr*asin((g2-r2)/(g2+r2))
      rlonn1=mod(orient+315,360.)
      rlats1=-rlatn1
      rlons1=mod(rlonn1+270,360.)
      jpds=-1
      jgds=-1
      do k=1,km
!        jpds(5)=kp5(k)
!        jpds(6)=kp6(k)
!        jpds(7)=kp7(k)
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,f(1,k),iret)
        print*, iret
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
      enddo
      idrt=kgds(1,1)
      imax=kgds(2,1)
      jmax=kgds(3,1)
c
      call sptruns(iromb,maxwv,idrt,imax,jmax,km,nps,
     &             0,0,0,jf,0,0,0,0,true,xmesh,orient,f,gn,gs)
c
      do k=1,km
        kpds(3,k)=27
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlatn1*1.e3)
        kgds(5,k)=nint(rlonn1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(orient*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=0
        kgds(11,k)=64
        call putgb(lun,npq,kpds(1,k),kgds(1,k),lb,gn(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=28
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlats1*1.e3)
        kgds(5,k)=nint(rlons1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(mod(orient+180,360.)*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=128
        kgds(11,k)=64
        call putgb(lus,npq,kpds(1,k),kgds(1,k),lb,gs(1,k),iret)
      enddo
c
      end
c
      subroutine gv65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c  interpolates a vector field using spectral transforms.
      integer kp5(km),kp6(km),kp7(km)
c  output grids are 65x65 (381 km true at latitide 60).
c  nh grid oriented at 280E; sh grid oriented at 100E.
c  winds are rotated to be relative to grid coordinates.
      parameter(nph=32,nps=2*nph+1,npq=nps*nps)
      parameter(true=60.,xmesh=381.e3,orient=280.)
      parameter(rerth=6.3712e6)
      parameter(pi=3.14159265358979,dpr=180./pi) 
      real un(npq,km),vn(npq,km),us(npq,km),vs(npq,km)
      integer jpds(25),jgds(22),kpds(25,km),kgds(22,km)
      logical lb(jf)
      real u(jf,km),v(jf,km)
c
      g2=((1.+sin(abs(true)/dpr))*rerth/xmesh)**2
      r2=2*nph**2
      rlatn1=dpr*asin((g2-r2)/(g2+r2))
      rlonn1=mod(orient+315,360.)
      rlats1=-rlatn1
      rlons1=mod(rlonn1+270,360.)
      jpds=-1
      do k=1,km
        jpds(5)=kp5(k)
        jpds(6)=kp6(k)
        jpds(7)=kp7(k)
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,u(1,k),iret)
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
        jpds=kpds(:,k)
        jgds=kgds(:,k)
        jpds(5)=jpds(5)+1
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,v(1,k),iret)
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
      enddo
      idrt=kgds(1,1)
      imax=kgds(2,1)
      jmax=kgds(3,1)
c
      call sptrunsv(iromb,maxwv,idrt,imax,jmax,km,nps,
     &              0,0,0,jf,0,0,0,0,true,xmesh,orient,u,v,
     &              .true.,un,vn,us,vs,.false.,dum,dum,dum,dum,
     &              .false.,dum,dum,dum,dum)
c
      do k=1,km
        kpds(3,k)=27
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlatn1*1.e3)
        kgds(5,k)=nint(rlonn1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(orient*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=0
        kgds(11,k)=64
        kpds(5,k)=kp5(k)
        call putgb(lun,npq,kpds(1,k),kgds(1,k),lb,un(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=27
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlatn1*1.e3)
        kgds(5,k)=nint(rlonn1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(orient*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=0
        kgds(11,k)=64
        kpds(5,k)=kp5(k)+1
        call putgb(lun,npq,kpds(1,k),kgds(1,k),lb,vn(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=28
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlats1*1.e3)
        kgds(5,k)=nint(rlons1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(mod(orient+180,360.)*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=128
        kgds(11,k)=64
        kpds(5,k)=kp5(k)
        call putgb(lus,npq,kpds(1,k),kgds(1,k),lb,us(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=28
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlats1*1.e3)
        kgds(5,k)=nint(rlons1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(mod(orient+180,360.)*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=128
        kgds(11,k)=64
        kpds(5,k)=kp5(k)+1
        call putgb(lus,npq,kpds(1,k),kgds(1,k),lb,vs(1,k),iret)
      enddo
c
      end
