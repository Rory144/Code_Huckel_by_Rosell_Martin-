!------------------------------------------------------------------------
! Copyright (c) [2024] [Rosell Martín Gómez]
!
! Permission is granted to use, modify and distribute this
! software, with or without modification, under the terms of the
! MIT License. See the LICENSE file for more details.
! more details.
!-----------------------------------------------------------------------

Module calculations

!----------------------------------------------------------------------------------------!   
! Module that contains two subroutines: 
!----------------------------------------------------------------------------------------

 implicit none 

 contains
 
!---------------------------------------------------------------------------------------
!           1. Subroutine to diagonalize the hamiltonian matrix and calculates 
!                          the eigenvalues and the eigenvectors 
!--------------------------------------------------------------------------------------
subroutine calculation1(N, hamil_matrix, E)

implicit none 

!.......................................................................................
!                           Definition of the variables 
!......................................................................................
! 1. Input variables
 
 integer :: N, M, N_atom  
 real(8), intent(inout) :: hamil_matrix(N, N)       
 real(8), intent(out) :: E(N)             
 
 ! hamil_matrix -> matrix that contains the eigenvectors 
 ! E -> Final vector that contains the eigenvalues. 

! 2. Local variables

 integer :: lwork, info, i
 real(8), allocatable :: work (:)

 ! work = temporary array. lwork is the size of such matrix
!........................................................................................

!Allocation
 allocate(work(3*N))
 lwork = size(work)

!Call the function dsyev 
call dsyev('V', 'U', N, hamil_matrix, N, E, work, lwork, info)

 if (info /= 0) then
      write(*,*) "Problem in diagonalize_matrix"
         stop
 end if

! dsyev provides the eigenvalues and eigenvectors of a symmetric matrix. 
! This function is part of LAPACK, a library used for numerical calculations in linear algebra.
! 'U' ->  the triangular upper part of the matrix is stored
! info = 0 means that there were NO errors 

!...........................................................................................
!        Write the eigenvalues and the eigenvectors in a external file 
!............................................................................................

open(78 , file = "linear_10000.txt", status = "old", action = 'write', position = 'append')
open(95, file = "data.txt", action = 'read')

read(95, *) 
read(95, *) 
read(95, *) M 
read(95, *) N_atom 

N = N_atom * M 

! 1. Write the eigenvalues 

 write(78, '(a)') adjustl(repeat('=', 75))
 write(78, '(9x, a)') 'The number of atoms      The normalization     Eigenvalues'
 write(78, '(a)') adjustl(repeat('-', 75))
     
        do i= 1, N
                write (78,'(*(F20.12))') real(i), real(i)/real(N), E(i)
        end do

! 2. Write the diagonalized matrix or eigenvectors
  write(78, '(a)') repeat('=', 75)
  write(78, '(23x, a)') 'The diagonalized matrix is or eigenvectors are: '
  write(78, '(a)') repeat('-', 75)

        do i = 1, N
                write(78,'(5x, 1000f12.6)') hamil_matrix(i,:)
        end do

 end subroutine calculation1

!===============================================================================
!                     2. Subroutine to compute the GAP  
!===============================================================================

subroutine calculation2(N, E) 

 implicit none 

!....................................................................................
!                          Definition of the variables 
!..................................................................................
 integer :: N, N_atom, M
 integer :: ne, q, i, HOMO, LUMO
 real(8) :: gap, Fermi_level                         
 real(8), allocatable, intent(in) :: E(:)          
 real(8), allocatable :: occup(:)                  

 ! occup -> the highest energy level that one electron can occupy at the absolute zero temperature.
 ! HOMO -> Highest Occupied Molecular Orbital.
 ! LUMO -> Lower Unoccupied Molecular Orbital.
 ! GAP -> The difference the energy between the HOMO and LUMO. 
 ! Fermi_level -> the energy at which there is a 50% probability of finding an electron occupying a state. 
!----------------------------------------------------------------------------------

open(78, file = "linear_10000.txt", status = "old", action = 'write', position = 'append')
open(90, file = "data.txt", action = 'read')

read(90, *) q
read(90, *) 
read(90, *) M 
read(90, *) N_atom 

 write(78, '(a)') repeat('=', 75)
 write(78, '(15x, a, i10)') ' The number of identical cells is : ', M

 write(78, '(a)') repeat('=', 75)
 write(78, '(15x, a, i10)') ' The number of different atoms is : ', N_atom

      N = N_atom * M

 write(78, '(a)') repeat('=', 75)
 write(78,'(18x, a, i10)') 'The number of total atoms is: ', N

       ne = N - q

 write(78, '(a)') adjustl(repeat('=', 75))
 write(78, '(20x, a, i10)') 'The number of electrons is: ', ne


!Allocation
allocate(occup(N))
occup = 0 

!In function of the electron filling, we will have:

!One electron per atom pair, even number of electrons:
  if (mod(ne, 2).eq.0) then
        do i = 1, ne/2
            occup(i) = 2
        end do
 end if

!One electron per atom, odd number of electrons:
if (mod(ne, 2).ne.0) then
  do i = 1, (ne - 1)/2
    occup(i) = 2
  end do 
  occup((ne + 1)/2) = 1 
end if 


!...........................................................................................
!                     Write the following in the external file 
!............................................................................................
! 1. Write the orbital, the occupation vector and the energies 

 write(78, '(a)') adjustl(repeat('=', 75))
 write(78, '(9x, a)') 'The Orbital      Occupation vector      Energy(eV) ' 
 write(78, '(a)') adjustl(repeat('-', 75))
 
         do i = 1, N
                write (78, '(6x, i10, 11x, f12.3, 9x, f12.6)') i, occup(i),  E(i)
         end do

 write(78, '(a)') adjustl(repeat('=', 75))

! 2. Compute the system gap and write the frontier orbital and the gap  
LUMO = 0

do i = 1, size(occup)
 if (occup(i).eq.0) then
        LUMO = i 
      exit
 end if
end do

HOMO = LUMO - 1

write(78,'(32x, a)') 'HOMO'
write(78, '(a)') adjustl(repeat('-', 75))
 write(78, '(25x, a)') 'Orbital      Energy(eV) '
 write(78, '(a)') adjustl(repeat('-', 75))
write(78, '(25x, i8, 10x, f16.5)') HOMO, E(HOMO)
write(78, '(a)') repeat('=', 75)

write(78,'(32x, a)') 'LUMO'
write(78, '(a)') repeat('-', 75)
write(78, '(25x, a)') 'Orbital     Energy(eV) '
write(78, '(a)') repeat('-', 75)
write(78, '(25x, i8, 10x, f16.5)')  LUMO, E(LUMO)

GAP = E(LUMO) - E(HOMO)

write(78, '(a)') repeat('=', 75)
write(78,'(20x, a, 1x, f9.5, 1x, a)') 'HOMO-LUMO GAP is:', GAP, 'eV'
write(78, '(a)') repeat('=', 75)

! 3. Compute the system gap at the Fermi Level and write it.

Fermi_level = E(HOMO) + GAP/2

write(78, '(15x, a, 1x, f12.5, 1x, a)') 'The GAP at the Fermi Level is:', Fermi_level, 'eV'
write(78, '(a)') adjustl(repeat('=', 75))

end subroutine calculation2 

end module calculations
