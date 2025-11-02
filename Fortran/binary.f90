program test
    integer :: i
    open(1, file='test.dat', access='stream', form='unformatted')
    do i = 1, 10
        write(1) i
        write(1) dble(i)**2
    end do
    close(1)
end program test 