        integer *4 n,k,i,bn,nbin,ct
        parameter (n=300000)
        real *8 binsize,ll,ul,x(n),y(n),vx(n),vy(n),v(n),p(n),mid
        k=1
        
        open(unit=10,file="singlep.d",status="old")
        open(unit=20,file="Pvx.d",status="unknown")

        do i=1,n
        read(10,*)x(i),y(i),vx(i),vy(i)
        v(i)=sqrt(vx(i)**2+vy(i)**2)
        enddo
        
        binsize=(vx(n)-vx(1))/float(n)
        do bn=1,n
        ll=vx(k)
        ul=vx(k)+binsize
        mid=(ll+ul)/2
        ct=0
        do i=1,n
        if(vx(i) .gt. ll .AND. vx(i) .le. ul) then
        ct=ct+1
        k=i
        endif
        enddo
        p(bn)=ct/binsize
        p(1)=ct/binsize + 1
        write(20,*) mid,p(bn)
        enddo
        
        stop
        end

