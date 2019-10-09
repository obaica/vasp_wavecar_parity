subroutine parity
    implicit none
    complex, parameter :: zi=(0.000, 0.000)
    integer :: i, j, l, n, m
    real:: V 
    complex, allocatable :: coeff(:, :)
    integer :: nspin, nkpt, nband, npmax
    integer, allocatable :: Gvector(:, :, :)
    real :: a1(3),a2(3),a3(3), b1(3),b2(3),b3(3)
    real :: kcurrent(3)
    real :: ktarget(3)
    integer :: bandindex, nplane
    complex :: energy
    real :: occ
    integer :: inot,rpt
    real :: r0(3), is(3)
    complex :: psi1, psi2
    real :: pindex
    complex, allocatable :: parityresult(:,:,:)
    real :: diff
    complex :: mp
    character :: pos

    
    a1=0.00;a2=0.00;a3=0.00
    b1=0.00;b2=0.00;b3=0.00
    rpt=5
    is=[0.00, 0.00, 0.00]

    open(unit=10, file='GCOEFF.txt', status='old')
    write(*,*) "Specified the target kpoint (kx, ky, kz):"
    read(*, *) ktarget(1), ktarget(2), ktarget(3)
    write(*,*) "Specified the Inversion Centre (x, y, z):"
    read(*,*) is(1), is(2), is(3)
    read(10, *) nspin
    read(10, *) nkpt
    read(10, *) nband
    read(10, *) V
    read(10, *) a1(1), a1(2), a1(3)
    read(10, *) a2(1), a2(2), a2(3)
    read(10, *) a3(1), a3(2), a3(3)
    read(10, *) b1(1), b1(2), b1(3)
    read(10, *) b2(1), b2(2), b2(3)
    read(10, *) b3(1), b3(2), b3(3)
    read(10, *) npmax
    allocate(coeff(nband, npmax))
    allocate(Gvector(nband, npmax, 3))
    allocate(parityresult(5,nband,2))
    Gvector=0.00
    coeff=0.00
    pindex=0
    
    do i=1, nkpt 
        inot=0 ! assume current k point not in target k list        
        read(10, *) kcurrent(1), kcurrent(2), kcurrent(3)
        do j=1, nband
            read(10, *) bandindex, nplane
            read(10, *) energy, occ
            do l=1, nplane
                read(10, *) Gvector(bandindex,l,1),Gvector(bandindex,l,2), Gvector(bandindex, l, 3), coeff(j, l)
        !        write(*, *) Gvector(bandindex, l, :)
            end do ! l=1, nplane
        end do ! j=1, nband

        !call inlist(kcurrent, ktarget, nktarget, inot)
        diff = sqrt((kcurrent(1)-ktarget(1))**2+(kcurrent(2)-ktarget(2))**2&
                     +(kcurrent(3)-ktarget(3))**2)
        if ( diff<=1d-5 ) then ! if the current kpt is target kpt
            write(*, *) 'current k point:', kcurrent
            write(*, *) 'target k point:', ktarget
            write(*,*) "#band index, parity, wf(+r)/wf(-r)"
            do l=1, nband
                do m=1, rpt
                    call rvector(r0)
                    call wavefun(Gvector(l,1:nplane,:),coeff(l,1:nplane),nplane,kcurrent, r0, is, V, psi1, psi2)
                    parityresult(m, l, 1)= psi1
                    parityresult(m, l, 2)= psi2
                end do
            end do ! l=1, nband   
            do l=1, nband
                mp=0.00
                do m=1,rpt
                    !if (mod(l,2).eq.0) cycle
                    mp=mp+parityresult(m,l,1)/parityresult(m,l, 2)
                end do
                mp = mp/dble(rpt) !take 5 points to avoid accident fault
                if (real(mp) < -0.99 .and. real(mp) > -1.01) then
                    pos="-"
                else if ( real(mp)>0.99 .and. real(mp)<1.01) then
                    pos="+"
                end if
                write(*,770)l, pos, mp
770 format (I10, A4, '   (', ES15.6, ', i', ES15.6, ')')
            end do
            exit
        end if ! if kcurrent in ktarget
    end do ! i=1, npkt
    close(10)
    write(*,*)"!!!!!! IF the parity of bands are not +/- or not showing "
    write(*,*)"!!!!!! or wf(+)/wf(-) not +/- (1., i0.)"
    write(*,*)"!!!!!! you need reconsider the INVERSION CENTRE!"
end subroutine parity

subroutine wavefun(G, C, npl, k ,r0, is, v, psi1, psi2)
! this subroutine aims at construct wavefuncion giving coefficients

    implicit none
    integer :: i
    complex, parameter :: pzi=2.0*3.141592654*(0.00,1.0000)
    complex, parameter :: zi=(0.00, 1.00)
    integer, intent(in) :: npl
    integer, intent(in) :: G(npl, 3)
    complex, intent(in) :: C(npl)
    real, intent(in) :: k(3), r0(3), v, is(3)
    real :: r1(3)
    complex, intent(out) :: psi1, psi2
    real :: ck(3)
    real :: dotpd0, dotpd1
    
    psi1=(0.000,0.000)
    psi2=(0.000,0.000)
    r1=2.0*is-r0
    do i=1,npl
        ck = G(i, :) + k
        dotpd0 = ck(1)*r0(1) + ck(2)*r0(2) + ck(3)*r0(3)
        dotpd1 = ck(1)*r1(1) + ck(2)*r1(2) + ck(3)*r1(3)
        psi1 = psi1 + C(i)*exp(pzi*dotpd0)/sqrt(v)
        psi2 = psi2 + C(i)*exp(pzi*dotpd1)/sqrt(v)
    end do
end subroutine wavefun

subroutine inlist(e, list, l, r)
!
!PURPUSE: check if element E in LIST or not
!PARAMETERS:
! e:    element
! list: target list
! l:    lenth of target list
! r:    result, r=1 for e IN list 
!
    implicit none
    integer, intent(in) :: l
    real,intent(in):: e(3), list(l, 3)
    integer, intent(inout) :: r
    integer :: n
    real:: tol, diff

    tol=1d-6
    diff=1d2
    
    do n=1, l
        diff=sqrt((e(1)-list(n,1))**2+(e(2)-list(n,2))**2+(e(3)-list(n,3))**2)
        if ( diff.lt.tol ) then
            r=1
            exit
        end if ! diff
    end do ! n=1, l

end subroutine inlist

subroutine rvector(r)
implicit none 
real, intent(out) :: r(3)
integer :: i
real :: x 
call random_seed()
do i=1,3
    call random_number(x)
    r(i)=2.0*x-1.0
end do
end subroutine rvector
