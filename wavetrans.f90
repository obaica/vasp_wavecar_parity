!!$************************* WaveTrans
!************************************
!!$*
!!$*   input the WAVECAR file in binary format from VASP, and output 
!!$*   in text format a file GCOEFF of G values and corresponding plane
!!$*   wave coefficients, together with energy eigenvalues and
!occupations
!!$*
!!$*   Compile with gfortran or ifort. Flag "-assume byterecl" is
!required
!!$*   for ifort.
!!$*
!!$*   version 2.0 - July 5, 2012 - R. M. Feenstra and M. Widom
!!$*   version 2.1 - Sept 30, 2012 - changed estimator for max. no. of
!!$*                                 plane waves
!!$*   version 1.3 - Sept 17, 2014 - updated 'c' value
!!$*
!!$*   format of GCOEFF file:
!!$*     no. wavevector values (integer)
!!$*     no. bands (integer)
!!$*     real space lattice vector (a1) x,y,z coefficients (real)
!!$*     real space lattice vector (a2) x,y,z coefficients (real)
!!$*     real space lattice vector (a3) x,y,z coefficients (real)
!!$*     recip. lattice vector 1 (b1) x,y,z coefficients (real)
!!$*     recip. lattice vector 2 (b2) x,y,z coefficients (real)
!!$*     recip. lattice vector 3 (b3) x,y,z coefficients (real)
!!$*     loop over spins
!!$*        loop through wavevector values:
!!$*           wavevector x,y,z components (real)
!!$*           energy (complex), occupation (real)
!!$*           loop through bands:
!!$*              index of band (integer), no. plane waves (integer)
!!$*              loop through plane waves:
!!$*                 ig1,ig2,ig3 values (integer) and coefficient
!(complex)
!!$*              end loop through plane waves
!!$*           end loop through bands
!!$*         end loop through wavevector values
!!$*     end loop over spins
!!$*
!!$*   the x,y,z components of each G value are given in terms of the
!!$*   ig1,ig2,ig3 values and the components of the recip. lattice
!vectors
!!$*   according to:
!!$*   ig1*b1_x + ig2*b2_x + ig3*b3_x,
!!$*   ig1*b1_y + ig2*b2_y + ig3*b3_y, and
!!$*   ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively
!!$*
!!$*   note that the energy eigenvalues are complex, as provided in the
!!$*   WAVECAR file, but the imaginary part is zero (at least for all
!cases
!!$*   investigated thus far)
!!$*     
subroutine wavetrans       
implicit real*8 (a-h, o-z)                                        
character :: soc
complex*8, allocatable :: coeff(:,:)
complex*16, allocatable :: cener(:)
real*8, allocatable :: occ(:)
integer, allocatable :: igall(:,:)
dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vtmp(3),sumkg(3),wk(3)
     
!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value;
!!$*   checks for discrepancy of any results between this and VASP
!values)

data c/0.262465831d0/ 
!!$*   data c/0.26246582250210965422d0/ 
pi=4.*atan(1.)
      
!!$*   input
soc='n'     
nrecl=24
do
    write(*,*)"Calculation with SoC ? 'y/n'"
    read(*,*) soc
    if ( soc.eq.'y' ) then
        write(*,*)"Calculation with SoC!"
        exit
    else if ( soc.eq.'n' ) then
        write(*,*)"Calculation without SoC!"
        exit
    else
        write(*,*)"Input error, y or n"
    end if
end do
    write(*,*)
open(unit=10,file='WAVECAR',access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost            
read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)
nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)
if(nprec.eq.45210) then
   write(6,*) '*** error - WAVECAR_double requires complex*16'
   stop
endif
write(6,*) 
write(6,*) 'record length  =',nrecl,' spins =',nspin, &
      'prec flag ',nprec
open(unit=10,file='WAVECAR',access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost            
open(unit=11,file='GCOEFF.txt')
read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)
nwk=nint(xnwk)
if (nwk.eq.0) write(*,*) 'no kpoints in WAVECAR or reading error'
nband=nint(xnband)
allocate(occ(nband))
allocate(cener(nband))
write(6,*) 'no. k points =',nwk
write(6,*) 'no. bands =',nband
write(6,*) 'max. energy =',sngl(ecut)
write(6,*) 'real space lattice vectors:'
write(6,*) 'a1 =',(sngl(a1(j)),j=1,3)
write(6,*) 'a2 =',(sngl(a2(j)),j=1,3)
write(6,*) 'a3 =',(sngl(a3(j)),j=1,3)
write(6,*) ' '
write(11,*) nspin
write(11,*) nwk
write(11,*) nband

!!$*   compute reciprocal properties

write(6,*) ' '
call vcross(vtmp,a2,a3)
Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
write(11,*) Vcell
write(6,*) 'volume unit cell =',sngl(Vcell)
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
do j=1,3
   b1(j)=2.*pi*b1(j)/Vcell
   b2(j)=2.*pi*b2(j)/Vcell
   b3(j)=2.*pi*b3(j)/Vcell
enddo
b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
write(6,*) 'reciprocal lattice vectors:'
write(6,*) 'b1 =',(sngl(b1(j)),j=1,3)
write(6,*) 'b2 =',(sngl(b2(j)),j=1,3)
write(6,*) 'b3 =',(sngl(b3(j)),j=1,3)
write(6,*) 'reciprocal lattice vector magnitudes:'
write(6,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)
write(6,*) ' '
write(11,*) (sngl(a1(j)),j=1,3)
write(11,*) (sngl(a2(j)),j=1,3)
write(11,*) (sngl(a3(j)),j=1,3)
write(11,*) (sngl(b1(j)),j=1,3)
write(11,*) (sngl(b2(j)),j=1,3)
write(11,*) (sngl(b3(j)),j=1,3)

phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
call vcross(vtmp,b1,b2)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
      
phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      
phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
if ( soc.eq.'y') then
    npmax=2*min0(npmaxA,npmaxB,npmaxC)
else if ( soc.eq.'n' ) then
    npmax=min0(npmaxA,npmaxB,npmaxC)
endif

write(6,*) 'max. no. G values; 1,2,3 =',nb1max,nb2max,nb3max
write(6,*) ' '

allocate (igall(3,npmax))
write(6,*) 'estimated max. no. plane waves =',npmax
write(11, *) npmax
allocate (coeff(npmax,nband))

!!$*   Begin loops over spin, k-points and bands
irec=2
do isp=1,nspin
   write(*,*) ' '
   write(*,*) '******'
   write(*,*) 'reading spin ',isp
   do iwk=1,nwk
      irec=irec+1
      read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
           (cener(iband),occ(iband),iband=1,nband)
      nplane=nint(xnplane)
      write(6,*) 'k point #',iwk,'  input no. of plane waves =', &
           nplane
      write(6,*) 'k value =',(sngl(wk(j)),j=1,3)
      write(11,*) (sngl(wk(j)),j=1,3)
      
!!$*   Calculate plane waves
      ncnt=0
      do ig3=0,2*nb3max
         ig3p=ig3
         if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
         do ig2=0,2*nb2max
            ig2p=ig2
            if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
            do ig1=0,2*nb1max
               ig1p=ig1
               if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
               do j=1,3
                  sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                       (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
               enddo
               gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
               etot=gtot**2/c
               if (etot.lt.ecut) then
                  ncnt=ncnt+1
                  igall(1,ncnt)=ig1p
                  igall(2,ncnt)=ig2p
                  igall(3,ncnt)=ig3p
               end if
            enddo
         enddo
      enddo
      if ( soc.eq.'y') then
         if (2*ncnt.ne.nplane) then
            write(6,*) 'ncnt=',ncnt, '// nplane=',nplane, 'npmax=', npmax
            write(6,*) '*** error - computed no. != input no.'
            stop
         endif
      else if ( soc.eq.'n' ) then
         if (ncnt.ne.nplane) then
            write(6,*) 'ncnt=',ncnt, '// nplane=',nplane, 'npmax=', npmax
            write(6,*) '*** error - computed no. != input no.'
            stop
         endif
      endif
      if (ncnt.gt.npmax) then
         write(6,*) '*** error - plane wave count exceeds estimate'
         stop
      endif
      
      do iband=1,nband
         irec=irec+1
         read(unit=10,rec=irec) (coeff(iplane,iband), &
              iplane=1,nplane)
      enddo
!!$*   output G values and coefficients
560      format('( ',g14.6,' , ',g14.6,' ) ',g14.6)            
570         format(3i6,'  ( ',g14.6,' , ',g14.6,' )')     
if ( soc.eq.'y' ) then
      do iband=1,nband
         write(11,*) iband,nplane
         write(11,560) cener(iband),occ(iband)
         do iplane=1,nplane
            write(11,570) (igall(j,mod(iplane-1, ncnt)+1),j=1,3), &
                 coeff(iplane,iband)
         enddo
      enddo
else
      do iband=1,nband
         write(11,*) iband,nplane
         write(11,560) cener(iband),occ(iband)
         do iplane=1,nplane
            write(11,570) (igall(j,iplane),j=1,3), &
                 coeff(iplane,iband)
         enddo
      enddo
end if
   enddo
enddo
end subroutine wavetrans

!!$*
!!$*   routine for computing vector cross-product
!!$*
subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross
      
      
