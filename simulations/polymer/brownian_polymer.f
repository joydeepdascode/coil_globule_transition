        integer *4 itime,ntime,idum,nparti,i,k
        parameter(ntime=300000,nparti=10)
        real *8 fx_vis(nparti),fy_vis(nparti),eta,dt,
     &  vx(nparti),x(nparti),vy(nparti),y(nparti),f_int(nparti),
     &  alpha,ran2,gasdev,temp_eta,D,T,beta,k_spr,
     &  fx_int(nparti),fy_int(nparti),l0,sum(ntime)
        
        parameter(alpha=1.0,k_spr=1.0,dt=0.01,l0=1.0)
        
        idum=13425879
        T=1.0
        beta=1./T
        D=alpha/beta
        temp_eta=dsqrt(2.0*D*dt)
        
        x(1)=0.0
        y(1)=0.0
        vx(1)=0.5
        vy(1)=0.5
        fx_int(1)=0.0
        fy_int(1)=0.0


        do i=2,nparti
                x(i)=float(i-1)*l0
                y(i)=0.0
                vx(i)=0.5*k_spr
                vy(i)=0.5*k_spr
                fx_int(i)=0.0
                fy_int(i)=0.0
                
        enddo

        open(unit=10,file="brownian_polymer.d",status="unknown")
        open(unit=20,file="qrownian_polymer.d",status="unknown")

        

        do itime=1,ntime
         do i=1,nparti-1
            write(10,*) x(i),y(i),vx(i),vy(i)
            fx_vis(i)=-alpha*vx(i)
            fy_vis(i)=-alpha*vy(i)

            j=i+1
            rr= (x(j)-x(i))**2 + (y(j)-y(i))**2
            rr= sqrt(rr)
            Fintsp= - k_spr*(rr-l0)
            Fintspx= Fintsp*(x(j)-x(i))/rr
            Fintspy= Fintsp*(y(j)-y(i))/rr
            fx_int(i)= fx_int(i) + Fintspx
            fy_int(i)= fy_int(i) + Fintspy
            fx_int(j)= fx_int(j) - Fintspx
            fy_int(j)= fy_int(j) - Fintspy

            eta = gasdev(idum)*temp_eta
            vx(i)=vx(i)+(fx_vis(i)+fx_int(i))*dt+eta
            vy(i)=vy(i)+(fy_vis(i)+fy_int(i))*dt+eta
            x(i)= x(i) + vx(i)*dt
            y(i)= y(i) + vy(i)*dt
         enddo
         write(10,*) x(nparti),y(nparti),vx(nparti),vy(nparti)
         vx(nparti)=vx(nparti)+(fx_vis(nparti)+fx_int(nparti))*dt+eta
         vy(nparti)=vy(nparti)+(fy_vis(nparti)+fy_int(nparti))*dt+eta
         x(nparti)= x(nparti) + vx(nparti)*dt
         y(nparti)= y(nparti) + vy(nparti)*dt

        enddo

        do i=1,nparti,10
        read(10,*)sum(i)

                


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


        



 
