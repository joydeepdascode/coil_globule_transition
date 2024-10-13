        integer *4 i,j,n
        parameter(n=5)
        real *8 x(n),y(n),a(n),b(n),c(n),z(n)

        open(unit=10,file="sum.d",status="old")
        open(unit=20,file="sum1.d",status="old")
        open(unit=30,file="add-data.d",status="unknown")

        do i=1,n
        read(10,*) x(i),a(i)
        read(20,*) y(i),b(i)
        c(i)=a(i)+b(i)
        z(i)=x(i)+y(i)
        write(30,*) z(i),c(i)
        enddo
        

        stop
        end



