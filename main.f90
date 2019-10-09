program main
    implicit none
    integer :: ierr

    ierr=0
    open(unit=19, file='GCOEFF.txt', status='old', action='read', iostat=ierr)
    close(unit=19)
    if ( ierr.eq.0 ) then
        write(*,*) "'GCOEFF.txt' file exits, now calc parity."
        call parity
        write(*,*) "**************End of parity calculation!*********************"
    else 
        write(*,*) "'GCOEFF.txt' file dosn't exit, now read WAVECAR."
        call wavetrans
        call parity
        write(*,*) "**************End of parity calculation!*********************"
    end if
end program main
