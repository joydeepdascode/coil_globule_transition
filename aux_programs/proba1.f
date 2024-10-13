        integer *4 n,k,i,bn,ibin,nbin
        parameter (n=300000)
        real *8 binsize,ll,ul,x(n),y(n),vx(n),vy(n),v(n),p(n),a_ct(n)
        real *8 asum1
        k=1
        vxmax=5.0
        vxmin=-5.0
        nbin = 200
        nx=0
        asum1=0

        
        open(unit=10,file="singlep.d",status="old")
        open(unit=20,file="pvx.d",status="unknown")
        open(unit=30,file="pvy.d",status="unknown")
        open(unit=40,file="pv.d",status="unknown")
        open(unit=50,file="px.d",status="unknown")
        open(unit=60,file="py.d",status="unknown")


        do i=1,n
        read(10,*)x(i),y(i),vx(i),vy(i)
        v(i)=sqrt(vx(i)**2+vy(i)**2)
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
        write(20,*) ((ibin-1)*binsize + vxmin + binsize/2.0),a_ct(ibin)
        asum1=asum1+a_ct(ibin)*binsize
        enddo
        write(*,*) asum1

        stop
        end
        
        

