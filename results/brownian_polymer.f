        implicit none
        integer *4 itime,ntime,idum,npart,i,j,k,itimeqb,g
        parameter(ntime=4000000,itimeqb=1000000,g=100)
        real *8 T,D,beta,alpha,T_eta,k_spr,rij,l0,ran2,gasdev,
     &  fx_vis,fy_vis,x(g),y(g),vx(g),vy(g),
     &  x_sep,y_sep,fx_int1(g),fy_int1(g),fx_spr,fy_spr,
     &  r,rnpart,rn,rn_avg,rn2_avg,xee,yee,xsumee,ysumee,x0,
     &  x2sumee,y2sumee,dt,vx2avg,vy2avg,etax,etay,sigma,
     &  fx_int2(g),fy_int2(g),V0,fint,fintx,finty,rinv,
     &  rinv2,rinv3,rinv6,rinv12
        parameter(alpha=1.0,k_spr=1000.0,dt=0.01,l0=1.0,sigma=1.0,
     &  V0=1.0)

        idum=13452879
        T=1.0
        beta=1.0/T
        D=alpha/beta
        T_eta=dsqrt(2.0*D*dt)
        rn_avg=0.0
        rn2_avg=0.0


        open(unit=10,file="InitalConfig.d",status="unknown")
        open(unit=20,file="FinalConfig.d",status="unknown")
        open(unit=30,file="Rn2.d",status="unknown",access='append')



        x0=0.0
        do npart=10,g,10
        do i=1,npart
          x(i)=x0+float(i)*l0
          y(i)=0.0
          vx(i)=dsqrt(2.0*D)*gasdev(idum)
          vy(i)=dsqrt(2.0*D)*gasdev(idum)
          fx_int1(i)=0.0
          fy_int1(i)=0.0
          fx_int2(i)=0.0
          fy_int2(i)=0.0
          write(10,*) x(i),y(i),vx(i),vy(i)
        enddo

        xsumee=0.0
        ysumee=0.0
        x2sumee=0.0
        y2sumee=0.0

        vx2avg=0.0
        vy2avg=0.0

        do itime=1,ntime

          do i=1,npart-1
              j=i+1
              x_sep=x(j)-x(i)
              y_sep=y(j)-y(i)
              rij= x_sep**2 + y_sep**2
              rij= dsqrt(rij)
              fx_spr= k_spr*(rij-l0)*(x_sep/rij)
              fy_spr= k_spr*(rij-l0)*(y_sep/rij)
              fx_int1(i)= fx_int1(i) + fx_spr
              fy_int1(i)= fy_int1(i) + fy_spr
              fx_int1(j)= fx_int1(j) - fx_spr
              fy_int1(j)= fy_int1(j) - fy_spr

            do k=i+1,npart
              if((k-i).gt.1) then
                x_sep=x(k)-x(i)
                y_sep=y(k)-y(i)
                r= x_sep**2 + y_sep**2
                r= dsqrt(r)
                rinv= 1.0/r
                rinv2= rinv*rinv
                rinv3=rinv2*rinv
                rinv6= rinv3*rinv3
                rinv12= rinv6*rinv6
                fint=V0*(12.0*(rinv12*rinv2)-6.0*(rinv6*rinv2))
                fintx= fint*x_sep
                finty= fint*y_sep
                fx_int2(i)= fx_int2(i) -fintx
                fy_int2(i)= fy_int2(i) -finty
                fx_int2(k)= fx_int2(k) +fintx
                fy_int2(k)= fy_int2(k) +finty
              endif
            enddo
          enddo

          do i=1,npart
            fx_vis= - alpha*vx(i)
            fy_vis= - alpha*vy(i)
            etay = T_eta*gasdev(idum)
            etax = T_eta*gasdev(idum)
            vx(i)=vx(i)+(fx_vis+fx_int1(i)+fx_int2(i))*dt+etax
            vy(i)=vy(i)+(fy_vis+fy_int1(i)+fy_int2(i))*dt+etay
            x(i)= x(i) + vx(i)*dt
            y(i)= y(i) + vy(i)*dt

            if(itime.gt.itimeqb) then
              vx2avg= vx2avg + vx(i)*vx(i)
              vy2avg= vy2avg + vy(i)*vy(i)
            endif
            fx_int1(i)=0.0
            fy_int1(i)=0.0
            fx_int2(i)=0.0
            fy_int2(i)=0.0
          enddo



          if(itime.gt.itimeqb) then
            xee= x(npart)- x(1)
            yee= y(npart)- y(1)
            xsumee= xsumee + xee
            ysumee= ysumee + yee
            x2sumee= x2sumee + xee*xee
            y2sumee= y2sumee + yee*yee
          endif


        enddo    !Time Loop ends here

        do i=1,npart
        write(20,*) x(i),y(i),vx(i),vy(i)
        enddo



        rn_avg=(dsqrt(xsumee**2 + ysumee**2))/float(ntime-itimeqb)
        rn2_avg=(x2sumee+y2sumee)/float(ntime-itimeqb)

        vx2avg= vx2avg/float((ntime-itimeqb)*npart)
        vy2avg= vy2avg/float((ntime-itimeqb)*npart)

        write(*,*) vx2avg,vy2avg


        write(30,*) npart,rn_avg,rn2_avg
        enddo

        stop
        end


        FUNCTION ran2(idum)
        INTEGER *4 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        real *8 ran2,AM,EPS,RNMX
        PARAMETER(IM1=2147483563,IM2=2147483399,AM=1.0/IM1,IMM1=IM1-1)
        PARAMETER(IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211)
        PARAMETER(IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7)
        PARAMETER(RNMX=1.-EPS)
        INTEGER *4 idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1)iy=iy+IMM1
        ran2=min(AM*iy,RNMX)
        return
        END function

        FUNCTION gasdev(idum)
        INTEGER *4 idum
        real *8 gasdev
        INTEGER *4 iset
        real *8 fac,gset,rsq,v1,v2,ran2
        SAVE iset,gset
        DATA iset/0/
        if (idum.lt.0) iset=0
        if (iset.eq.0) then
   1    v1=2.*ran2(idum)-1
        v2=2.*ran2(idum)-1
        rsq = v1**2 + v2**2
        if(rsq.ge.1..or.rsq.eq.0.) goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
        else
        gasdev=gset
        iset=0
        endif
        return
        END FUNCTION
