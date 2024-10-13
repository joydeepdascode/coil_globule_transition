        integer *4 i,n,idum
        real *8 vx,vy,x,y,alpha,eta_x,eta_y,dt,t
	real *8 gasdev,D,beta,tempeta,te
        parameter (n=300000,alpha=1.0)
        vx=0.0
        x=0.0
	vy=0.0
	y=0.0 
        dt=0.01
        t=0.0
	idum=13425879
	te=1.0
        beta=1./te
	D=alpha/beta
	tempeta=dsqrt(2.0*D*dt)

        open(unit=10,file="vx.d",status="unknown")
        open(unit=20,file="vy.d",status="unknown")
        open(unit=30,file="x.d",status="unknown")
        open(unit=40,file="y.d",status="unknown")
        open(unit=50,file="singlep.d",status="unknown")

        do i=0,n
        write(10,*)vx,t
        write(20,*)vy,t
        write(30,*)x,t
        write(40,*)y,t
        write(50,*)x,y,vx,vy
	eta_x=gasdev(idum)
	eta_y=gasdev(idum)
	vx= vx - alpha*vx*dt + eta_x*tempeta
        x= x + vx*dt
        vy= vy - alpha*vy*dt + eta_y*tempeta
        y= y + vy*dt
        t=t+dt
        
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

