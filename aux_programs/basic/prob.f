        integer *4 ibin,nbin,i,n
        parameter (n=1000)
        real *8 x(n),xrange,binsize,acount(50)
        do i=1,n
        call random_seed()
        call random_number(x(i))
        enddo

        xrange=1.0
        nbin=50
        binsize=xrange/float(nbin)
        
        do ibin=1,nbin
        acount(ibin)=0.0
        enddo
        
        do i=1,n
        ibin=int(x(i)/binsize)+1
        if(ibin.gt.nbin) ibin=nbin
        acount(ibin)=acount(ibin)+1.0
        enddo

       

        do ibin=1,nbin
        acount(ibin)=acount(ibin)/float(n)
        write(10,*) (ibin-1)*binsize+binsize/2.0,acount(ibin)
        enddo

        stop
        end



