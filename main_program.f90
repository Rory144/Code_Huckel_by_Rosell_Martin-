!------------------------------------------------------------------------
! Copyright (c) [2024] [Rosell Martín Gómez] 
!
! Permission is granted to use, modify and distribute this
! software, with or without modification, under the terms of the
! MIT License. See the LICENSE file for more details.
! more details.
!-----------------------------------------------------------------------

!....................................................................
! Program to compute and diagonalize the hamiltonian matrix corresponding 
!                 to a 1D polyene(linear or cyclic)
!...................................................................

program hamiltonian_matrix
      
! This program use the module called calculations: 
  use calculations 

  implicit none

!.................................................................
!                     Definition of variables
!.................................................................
integer :: i, j, k, g, p, un, N, N_atom, M, ne
integer :: ios 
character(len = 10)  :: chain
real(8) :: q 
real(8), allocatable :: hamil_matrix(:,:), E(:), alpha(:), beta(:)
  
! N -> is the dimension of the hamiltonian matrix 
! alpha -> diagonal elements of the hamiltonian matrix 
! beta(+) -> Adjacent diagonal elements at the upper part
! beta(-) -> Adjacent diagonal elements at the lower part
! chain ->   variable to define the shape of the chain of polyene
! hamil_matrix -> The hamiltonian matrix .
! E -> The matrix of eigenvalues.
! N_atom -> Number of different atoms in the polyene.
! M -> Number of identical cells.
! N -> dimension of the hamiltonian matrix.
! q ->  charge of the polyene.
! ne -> number of electrons.
!.................................................................

!..........................................................................................
!   First, read the variables from a external file and write the results and another file
!...........................................................................................

 open(78, file = "linear_10000.txt", status = "new", action = 'write')

 open(88, file = 'data.txt', action = 'read', status = 'old', iostat = ios)
 if (ios/= 0) then
         print *, "Error opening file 'data.txt'"
    stop
  end if

 read(88, *) 
 read(88, *) chain
 read(88, *) M
 read(88, *) N_atom

 N = N_atom * M 
!.............................................................................
!In this code, two types of dimerization are considered: 
!.............................................................................
!                         1. Bond alternation
!..............................................................................
!Two diffentes values of beta are considered: betha(+) and beta(-)
 allocate(beta(2))

 read(88, *) beta(1) 
 read(88, *) beta(2) 
!.........................................................................
!                         2. Atom alternation 
!........................................................................
!Diffent types of atoms or different values of diagonal parameter -> alpha
 allocate(alpha(N_atom))
 alpha = 0.d0 

 do i = 1, N_atom
        read(88, *) alpha(i)
 end do  

!..............................................................................
!   Inizialize and allocation of the hamiltonian matrix and eigenvalues vector
!................................................................................
 allocate(E(N))
 allocate(hamil_matrix(N, N))
 hamil_matrix = 0.d0 
!.................................................................................
!                   Definition of the alpha and beta vectors 
!................................................................................
k = 1
g = 1
p = 1

!When iterate for each row, i, and each column, j of the hamiltonian matrix, do the following: 

    do i = 1, N
       do j = 1, N
           
          ! Definition of alpha: 
           if (i.eq.j) then 
       
             hamil_matrix(i,j) = alpha(k)
                
               k = k + 1
                  
               if (k.gt.N_atom) then 
          
                   k = 1
            
               end if  
           
           end if 
     
           ! Definition for the two types of beta: 
           if (abs(i-j).eq.1) then
              
              !This is for those beta that is find in the upper part
               if (i.lt.j) then 
               
                    hamil_matrix(i,j) = beta(g)
               
                    g = g + 1
              
                        if (g.gt.2) then 
                            g = 1 
                        end if 
               end if 
                        
               !This is for beta that are in the lower part: 
               if (i.gt.j) then  
       
                   hamil_matrix(i,j) = beta(p)
                   
                   p = p + 1
                     
                     if (p.gt.2) then 
                     
                        p = 1

                     end if 
               end if 
           end if 
                   
        end do 
     end do 
          

!Consider another shape of the polyene: if the polyene is cyclic, write the word 'cyclic' in the input:  
  
     if (trim(chain).eq.'cyclic') then 
        
              hamil_matrix(1, N) = beta(2)
         
              hamil_matrix(N, 1) = beta(2)
     end if 

! Write the hamiltonian  matrix in a external file 
 
   write(78, '(a)') repeat('=', 75)
     write(78,'(21x, a)') 'The hamiltonian matrix is: ' 
   write(78, '(a)') repeat('-', 75)
     
   do i = 1, N
       write(78,'(1000f20.8)')  hamil_matrix(i,:)
   end do 

!Call the subroutines

call calculation1(N, hamil_matrix, E) 
call calculation2(N, E) 

close(88) 
close(78) 

end program hamiltonian_matrix
