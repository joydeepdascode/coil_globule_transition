        implicit none
        integer *4 itime,ntime,idum,npart,i,j,itimeqb,k,l
        parameter(ntime=4000000,k=100,npart=10,itimeqb=1000000)
        real *8 T,D,beta,alpha,T_eta,k_spr,rij,l0,ran2,gasdev,
     &  fx_vis(k),fy_vis(k),x(k),y(k),vx(k),vy(k),
     &  x_sep,y_sep,fx_int(k),fy_int(k),fx_spr,fy_spr,
     &  r1,rnpart,Rn,Rn_avg,Rn2_avg,xee,yee,xsumee,ysumee,x0,
     &  x2sumee,y2sumee,dt,vx2avg,vy2avg,x2avg(k),etax,etay,
     &  flen,fx_len(k),fy_len(k),ril,xsep,ysep,a,b
        parameter(alpha=1.0,k_spr=1.0,dt=0.01,l0=1.0,a=1.0,b=1.0)

        idum=13452879
        T=1.0
        beta=1.0/T
        D=alpha/beta
        T_eta=dsqrt(2.0*D*dt)
        Rn_avg=0.0
        Rn2_avg=0.0
        

        open(unit=10,file="Rn2.d",status="unknown")
        open(unit=20,file="eq.d",status="unknown")

        
        x0=0.0
        do i=1,npart
                x(i)=x0+float(i)*l0
                y(i)=0.0
                vx(i)=dsqrt(2.0*D)*gasdev(idum)
                vy(i)=dsqrt(2.0*D)*gasdev(idum)
                fx_int(i)=0.0
                fy_int(i)=0.0
                fx_vis(i)=0.0
                fy_vis(i)=0.0
                fx_len(i)=0.0 
                fy_len(i)=0.0 
                

        enddo

        xsumee=0.0
        ysumee=0.0
        x2sumee=0.0
        y2sumee=0.0

        vx2avg=0.0
        vy2avg=0.0

        do itime=1,ntime
        do i=1,npart
           if(i.eq.npart) goto 2
                   l=i+1
                        x_sep=x(l)-x(i)
                        y_sep=y(l)-y(i)
                        ril= x_sep**2 + y_sep**2
                        ril= dsqrt(rij)
                        fx_spr= k_spr*(ril-l0)*(x_sep/ril)
                        fy_spr= k_spr*(ril-l0)*(y_sep/ril)
                        fx_int(i)= fx_int(i) + fx_spr
                        fy_int(i)= fy_int(i) + fy_spr
                        fx_int(l)= fx_int(l) - fx_spr
                        fy_int(l)= fy_int(l) - fy_spr
2                       fx_vis(i)= - alpha*vx(i)
                        fy_vis(i)= - alpha*vy(i)
                                                
                        
                        etay = T_eta*gasdev(idum)
                        etax = T_eta*gasdev(idum)
             do j=1,npart                                  
                 if(j.ne.i-1 .and. j .ne. i+1.and.j.ne.i) then
                 xsep=x(j)-x(i)
                 ysep=y(j)-y(i)
                 rij=xsep**2+ysep**2
                 rij=dsqrt(rij)
                 flen=((-12*a)/rij**13)+((6*b)/rij**7)
                 fx_len(i)=fx_len(i)+flen*(xsep/rij)
                 fy_len(i)=fy_len(i)+flen*(ysep/rij)
                 endif
             enddo
             if(itime.gt.itimeqb) then
                vx2avg = vx2avg + vx(i)*vx(i)
                vy2avg = vy2avg + vy(i)*vy(i)
             endif
        enddo
        
        do i=1,npart
           vx(i)=vx(i)+(fx_vis(i)+fx_int(i)+fx_len(i))*dt+etax
           vy(i)=vy(i)+(fy_vis(i)+fy_int(i)+fy_len(i))*dt+etay
           x(i)= x(i) + vx(i)*dt
           y(i)= y(i) + vy(i)*dt
        enddo
               
                if(itime.gt.itimeqb) then
                xee= x(npart)- x(1)
                yee= y(npart)- y(1)
                xsumee= xsumee + xee
                ysumee= ysumee + yee
                x2sumee= x2sumee + xee*xee
                y2sumee= y2sumee + yee*yee
                

                endif

        do i=1,npart
         fx_int(i)=0.0
         fy_int(i)=0.0
         fx_len(i)=0.0
         fy_len(i)=0.0
        enddo


        enddo
       
        do i=1,npart
        write(30,*) x(i),y(i),vx(i),vy(i)
        enddo 



        Rn_avg=(dsqrt(xsumee**2 + ysumee**2))/float(ntime-itimeqb)
        Rn2_avg=(x2sumee+y2sumee)/float(ntime-itimeqb)
        
        vx2avg= vx2avg/float((ntime-itimeqb)*npart)
        vy2avg= vy2avg/float((ntime-itimeqb)*npart)

        write(*,*) vx2avg,vy2avg
 
        
        write(10,*) npart,Rn_avg,Rn2_avg
        
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
   1                    v1=2.*ran2(idum)-1
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

