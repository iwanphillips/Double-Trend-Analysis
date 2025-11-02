program cp

    ! tries to find the most plausible change point
    ! in linear interpolation.
    ! inspired by Bradley's way of identification of fractal dimension
    
    ! test for all possible break points, fit two slopes on each subset,
    ! and calculate the rms errors.
    
    ! Here+ the two linear fits are requied to merge continuously
    ! at the break point. So the minimization is a constraint minimization
    ! and has to be done simultaneously for the two segments. 
    
    ! modification of previous versions: the time index in the second segment is
    ! set to start again at 0, for numerical stability this is essential
    
    implicit real*8 (a-h,o-z)
    parameter (n=50000)
    
    real*8 x(n)
    iseed1=4235243
    iseed2=64424321
    
    !l=50000
    !open(9, file='tseries.bin', access='stream', form='unformatted')
    !do i=1,l
    !    call gauss(y1,y2,iseed1,iseed2)
    !    x(i)=y1
    !    if (i.gt.25000) then
    !        x(i)=x(i)+(i-25000)*.000125
    	    !write(9,*)real(x(i))
    !    end if
    !    write(9) real(x(i))
    !enddo
    !close(9)

    !print*,l,' gaussian data with trend generated'
    
!    open(8, file='Potsdam_T2m_av_anom.dat')
!        do i=1,n
!            read(8, *, end=3) x(i)
!            write(9,*)real(xx)
!        enddo

    !integer :: io
    !open(8, file='Potsdam_T2m_av_anom.dat', status="old", action="read")
    !open(8, file='pot_.dat', status="old", action="read")
    open(8, file='potsd.dat', status="old", action="read") ! , status="old", action="read"
    open(9, file='pot_check.bin', status='replace', access='stream', form='unformatted')
    !read(*,*)
    do i=1,n
        !read(8, *, end=3) x(i)
        read(8, *, end=3) x(i)
        !write(9,*)real(xx)
        !write(9) real(xx)
        write(9) x(i)
        !print*, x(i)
    enddo
    close(8)
    close(9)

    3 continue !!! In Holger's code there was no continue !!!
    l=i-1
    print*,l,' data read from Potsdam_T2m_av_anom.dat'
    
    errmin=l*50 ! Maybe this should be made bigger l*50 -> l*500 
                ! due to having bigger variance from not detracting seasonal cycle
    
    open(10, file='error.bin', status='replace', access='stream', form='unformatted')
    do ic=365,l-365,365
    ! do ic=1000,l-1000,1000
        call linfit_c(x,ic,l,a1,b1,a2,b2,err)
        ! write(10,*)ic,real(err),real(a1),real(a2)
        write(10) real(err)
        write(12,*)real(a1),real(b1),real(a2),real(b2)
        write(13,*)0,b1
        write(13,*)ic,ic*a1+b1
        write(13,*)ic,b2
        write(13,*)l,b2+(l-ic)*a2
        write(13,*)
        if (err.lt.errmin) then
            errmin=err
            icopt = ic
            a1opt=a1
            a2opt=a2
            b1opt=b1
            b2opt=b2
        endif
    enddo
    close(10)

    print*,'optimal values found:'
    print*,'errmin =',errmin
    print*,'changepoint at',icopt
    print*,'slope and offset of first segment',a1opt,b1opt
    print*,'slope and offset of second segment',a2opt,b2opt
    !print*,'and in K/century:',a1opt*36525,a2opt*36525

    open(11, file='fit_results.bin', status='replace', access='stream', form='unformatted')
    write(11) dble( errmin )
    write(11) dble( icopt ) 
    write(11) dble( a1opt )
    write(11) dble( b1opt )
    write(11) dble( a2opt )
    write(11) dble( b2opt )
    close(11)
        
    write(14,*)0,b1opt
    write(14,*)icopt,icopt*a1opt+b1opt
    write(14,*)icopt,b2opt
    write(14,*)l,b2opt+(l-icopt)*a2opt
    write(14,*)
    print*,'optimal fit in fort.14'
    stop

end program cp 

!>-------------------------------------------------------

subroutine gauss(r1,r2,iseed1,iseed2)
    
    real*8 r1,r2,p,phi,r
    pii=8.d0*atan(1.d0)
    
    call RANDOM1(p,iseed1)
    call RANDOM1(phi,iseed2)
           
    phi=phi*pii
    r=sqrt(-log(1.d0-p)*2.d0)
    
    r1=r*sin(phi)
    r2=r*cos(phi)
    return
end subroutine gauss

subroutine RANDOM1(r,iseed)
    
    ! random number generator of Park & Miller
    ! random numbers in [0,1] !!!
    real*8 r
    integer*8 ia,im,ix
    ia=7**5
    im=2147483647
    ix=iseed
    ix=mod(ia*ix,im)
    r=dfloat(ix)/dfloat(im)
    iseed=ix
    return
end subroutine RANDOM1

subroutine linfit_c(x,k,l,a1,b1,a2,b2,err)
    implicit real*8 (a-h,o-z)
    real*8 x(50000)
    
    xisum1=0.d0
    xsum1=0.d0
    do i=1,k
        xsum1=xsum1+x(i)
        xisum1=xisum1+dfloat(i)*x(i)
    enddo

    xisum1=xisum1/dfloat(k)
    xsum1=xsum1/dfloat(k)
          
    xisum2=0.d0
    xsum2=0.d0

    do i=k+1,l
        xsum2=xsum2+x(i)
        xisum2=xisum2+dfloat(i-k)*x(i)
    enddo

    xisum2=xisum2/dfloat(l-k)
    xsum2=xsum2/dfloat(l-k)
    ! print*,'the sums'
    ! print*,xisum1,xsum1,xisum2,xsum2
    
    g1=1.d0/dfloat(l-k)+dfloat(l-k-1)/dfloat(4*k)/(dfloat(l-k)+.5d0)
    g2=dfloat(k+1)/3.d0
    g3=dfloat(l-k-1)/4.d0/(dfloat(l-k)+.5d0)
    g4=dfloat(k+1)/dfloat(l-k)/2.d0
    g5=3.d0*xisum2/2.d0/dfloat(k)/(dfloat(l-k)+.5d0)
    g6=2.d0/dfloat(k-1)*(dfloat(k)*xsum1-xisum1)
    
    a1=(xsum1/dfloat(l-k)+xsum2/dfloat(k)-g5-g6*g1)/(g4+g3-g2*g1)
    b1=2.d0/dfloat(k-1)*(dfloat(k)*xsum1-xisum1-a1*dfloat(k-1)*float(k+1)/6.d0)
    b2=b1+a1*dfloat(k)
    a2=3.d0*xisum2/dfloat(l-k+1)/(dfloat(l-k)+.5d0)-b2*3.d0/2.d0/(dfloat(l-k)+.5d0)
    ! print*,'the fitted parameters at cp=',k
    ! print*,a1,b1,a2,b2
    
    ! test successful: the lambdas come out identically????
    
    ! and check the original equations whether they are satified
    ! slambda1=(xsum1-a1*(k+1)*.5d0-b1)*dfloat(k)
    ! slambda2=(xsum2-a2*dfloat(l-k+1)*.5d0-b2)*dfloat(l-k)
    ! slambda1a=xisum1-a1*(k+1)*(k+.5d0)/3.d0-b1*(k+1)*.5d0
    ! zero=xisum2-a2*(l-k+1)*(l-k+.5d0)/3.d0-b2*(l-k+1)/2.d0
    ! print*,'the 3 lambdas'
    ! print*,slambda1,slambda1a,slambda2
    ! print*,'check the zero of the third equation',zero
    ! print*,'the parameters',a1,b1,a2,b2
    
    err=0.d0
    do i=1,k
        err=err+(x(i)-dfloat(i)*a1-b1)**2
    enddo

    do i=k+1,l
        err=err+(x(i)-dfloat(i-k)*a2-b2)**2
    enddo
    ! print*,err,'error when splitting at',k
    
    return
end subroutine linfit_c
















            ! read(8,*,end=3)xx
            !read(8, *, end=3) xx
            !3 continue
            !x(i)=xx