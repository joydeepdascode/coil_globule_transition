        integer *4 i,n,idum,bn
        real *8 alpha,eta_x,dt,t
        real *8 vx1,x1,vx2,x2,vx3,x3,vx4,x4,vx5,x5     
        real *8 gasdev,D,beta,tempeta,te,l,k
        parameter (n=300000,alpha=1.0,bn=5)
        vx1=0.0
        x1=0.0
        vx2=0.0
        x2=0.0
        vx3=0.0
        x3=0.0
        vx4=0.0
        x4=0.0
        vx5=0.0
        x5=0.0
        
        dt=0.01
        time=0.0
        idum=13425879
        T=1.0
        beta=1./T
        D=alpha/beta
        tempeta=dsqrt(2.0*D*dt)
        l=1.0
        k=1.0
        open(unit=10,file="vx.d",status="unknown")
        open(unit=20,file="x.d",status="unknown")
        
        do itime=1,ntime
        !write(10,*)vx1,vx2,vx3,vx4,vx5
        eta_x=gasdev(idum)
        vx1= vx1 - alpha*vx1*dt -k*(x2-x1-l) + eta_x*tempeta
        vx2= vx2 - alpha*vx2*dt -k*(x3-x1-2*x2) + eta_x*tempeta
        vx3= vx3 - alpha*vx3*dt -k*(x4-x2-2*x3) + eta_x*tempeta
        vx4= vx4 - alpha*vx4*dt -k*(x5-x3-2*x4) + eta_x*tempeta
        vx5= vx5 - alpha*vx5*dt -k*(x4-x5-l) + eta_x*tempeta
        !write(20,*)x1,x2,x3,x4,x5
        x1=x1+vx1*dt
        x2=x2+vx2*dt
        x3=x3+vx3*dt
        x4=x4+vx4*dt
        x5=x5+vx5*dt
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
1                       v1=2.*ran2(idum)-1
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

