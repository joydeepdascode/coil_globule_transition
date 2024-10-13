        integer *4 i,n
        real *8 x,a
        
        n=10
        a=0.0

        open(unit=10,file="sum.d",status="unknown")
        open(unit=20,file="sum1.d",status="unknown")

        
        do i=1,n
        x=float(i)
        a=a+x
        if(mod(i,2).eq.0) then
        write(10,*) x,a
        endif
        if(mod(i,2).ne.0) then
        write(20,*) x,a
        endif
        
        enddo

        write(*,*) a

        stop
        end
        
