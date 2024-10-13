c       USED TEMPERATURE  T= 0.0
c	asum1 represents average of vx**2
c	asum2 represents average of x**2
        integer *4 n,k,i,ibin,nbin,nx
        parameter (n=300000)
        real *8 binsize,x(n),vx(n),a_ct(n),vxmax,vxmin,xmax,xmin
        real *8 asum1,u,q,asum2
        k=1
        vxmax=5.0
        vxmin=-5.0
        xmax=5.0
        xmin=-5.0
        nbin = 200
        nx=0
        asum1=0
        asum2=0
        
        open(unit=10,file="x_vx_data.d",status="old")
        open(unit=20,file="pvx.d",status="unknown")
        open(unit=50,file="px.d",status="unknown")
        


        do i=1,n
        read(10,*)x(i),vx(i)
        enddo

        binsize=(vxmax-vxmin)/float(nbin)

        do ibin=1,nbin
        a_ct(ibin)=0
        enddo

        do i=1,n
        if(vx(i) .ge. vxmin .and. vx(i) .le. vxmax) then
        ibin=int((vx(i)-vxmin)/binsize)+1
        endif
        if(ibin .ge. nbin) ibin=nbin
        a_ct(ibin)=a_ct(ibin)+1
        nx=nx+1
        enddo

        
        do ibin=1,nbin
        a_ct(ibin)=a_ct(ibin)/(float(nx)*binsize)
        u=((ibin-1)*binsize + vxmin + binsize/2.0)
        write(20,*) u,a_ct(ibin)
        asum1=asum1+(u**2)*a_ct(ibin)*binsize
        enddo
        write(*,*) asum1

        do ibin=1,nbin
        a_ct(ibin)=0
        enddo
        
        nx=0
        
        do i=1,n
        if(x(i) .ge. xmin .and. x(i) .le. xmax) then
        ibin=int((x(i)-xmin)/binsize)+1
        endif
        if(ibin .ge. nbin) ibin=nbin
        a_ct(ibin)=a_ct(ibin)+1
        nx=nx+1
        enddo

        
        do ibin=1,nbin
        a_ct(ibin)=a_ct(ibin)/(float(nx)*binsize)
        q=((ibin-1)*binsize + xmin + binsize/2.0)
        write(50,*) q,a_ct(ibin)
        asum2=asum2+(q**2)*a_ct(ibin)*binsize
        enddo
        write(*,*) asum2


        stop
        end
