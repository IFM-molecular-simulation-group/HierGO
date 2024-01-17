program read
implicit none

real,  dimension(:), allocatable :: x, y, z, x1, y1, z1
integer, dimension(:), allocatable :: nom, num, binsum
character(len=4), dimension(:), allocatable :: Resname
character(len=3), dimension(:), allocatable :: Atomtype, Atomtype1, atom1	
integer :: i, bintot, count,j, nat,counter, count1,index,p
character(len=50) :: format1 = "(I5,A4,3X,A3,I5,2X,F6.3,2X,F6.3,2X,F6.3)"
real :: alat,blat,clat
character(len=5), dimension(:), allocatable::residue,residue1, opls 
real, dimension(:), allocatable:: charge1
integer, dimension(:), allocatable:: oplsnum
integer:: k, natms, totatms,oplstypes,q



open(40, &
file='./output.gro',&
     status='old',action= 'read')                                                !!!incomplete gro file
open(50,&
file='./Oplsatom2', &
     status='old',action= 'read')                                               !!! file with OPLS atoms
open(60,& 
file='./result.gro', &
     status='new', action='write')


11 format (f10.5,f10.5,f10.5)
!12 format (A3,3X, f10.3, 3X, f10.3, 3X, f8.3 )
2 format (I6) 


oplstypes=17      !!!number of OPLS atoms types

!3 format (2X,A3)  

     
 !!!!#########read gro file    
   read(40,*) 
     read(40,*) nat
   
     allocate(nom(nat))
     allocate(num(nat))
     allocate(x1(nat))
     allocate(y1(nat))
     allocate(z1(nat))
     allocate(Atomtype1(nat))
     allocate(Resname(nat))  
     
     do i = 1, nat
       read(40,format1)nom(i), Resname(i), Atomtype1(i), num(i),x1(i), y1(i), z1(i)
     end do
       
    read(40,11) alat, blat, clat

!!!!!!####!!!!!!!!!!!!!!!!!read OPLS atoms and charge    Edit atomtypes and charges in input file 

4 format (A5,I3,2X,A4,1X,A3,2x,F7.4) 

     allocate(opls(oplstypes))
     allocate(oplsnum(oplstypes))
     allocate(residue1(oplstypes))
     allocate(atom1(oplstypes))
     allocate(charge1(oplstypes))

do q = 1, oplstypes
  read(50,4) opls(q), oplsnum(q), residue1(q), atom1(q), charge1(q) 
end do 


!!!!!!!!!!!writing new gro file
write(60,*) 'Graphene Oxide'
Write(60,2) nat  
 

  do i= 1,nat
    do q = 1, oplstypes
      if (Atomtype1(i).eq. atom1(q))   then
        write(60,format1)nom(i), adjustl(residue1(q)), Atomtype1(i), num(i), x1(i), y1(i), z1(i)
       end if
    end do 
  end do      


write(60,11) alat, blat, clat
   
close(40)
close(50)
close(60)

end program read
