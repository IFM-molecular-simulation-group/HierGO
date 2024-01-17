program gromacs_top
implicit none

real, dimension(:), allocatable:: charge, mass, charge1
integer, dimension(:), allocatable:: nr, resnr, cgnr,oplsnum
character(len=16)::header, titlebond, titleangles  
character(len=18)::comment,titledihedrals,titleimpropers
character(len=15)::molecule 
character(len=11)::atoms, bonds,angles
character(len=50)::titles 
character(len=13)::dihedrals,impropers


character(len=8),dimension(:), allocatable:: atmtype
character(len=3), dimension(:), allocatable:: atom, atom1, atom2
character(len=5), dimension(:), allocatable::residue,residue1, opls  
 
integer,dimension(:), allocatable:: bondi, bondj, funcb,bond2i, bond2j, func2b
integer, dimension(:), allocatable:: angi, angj, angl, funca
integer, dimension(:), allocatable:: dihedi, dihedj, dihedk, dihedl, funcd
integer,dimension(:), allocatable:: impropi, impropj, impropk, impropl, funci  
integer:: n, i,j,k,l,p,q,natoms, nbonds, nangles, count, nbonds1, ndihedrals,nimpropers,oplstypes
integer:: func,  funcp, nrexcl 
real:: sum1
character(len=80):: format1= "(I6,9X,A3,6X,I1,5X,A4,4X,A3,1X,I6,4X,F7.4,4X,F7.4)" 
character(len=100):: format2= "(I6,10X,A5,I3,6X,I1,5X,A4,4X,A3,1X,I6,4X,F7.4,4X,F7.4)"   



natoms=4477        !!!number of atoms
nbonds=6201      !!!number of bonds 
nangles=12448      !!!number of angles
ndihedrals=24528   !!!number of dihedrals
nimpropers=2978    !!!number of impropers
oplstypes=17      !!!number of OPLS atoms types




open(20,&
file='./Oplsatom2', &   !!!!!file with OPLS atoms
     status='old',action= 'read')
open(30,& 
file='./filename.top', &   !!!!half-baked tolpology
     status='old',action= 'read')
open(40, &
file='./GO_sheet_stitch.itp', & 
     status='new',action= 'write')




 
!!!!!!!!!!!!!!!!!read OPLS atoms and charge    Edit atomtypes and charges in input file 
4 format (A5,I3,2X,A4,1X,A3,2x,F7.4) 

     allocate(opls(oplstypes))
     allocate(oplsnum(oplstypes))
     allocate(residue1(oplstypes))
     allocate(atom1(oplstypes))
     allocate(charge1(oplstypes))
     
do q = 1, oplstypes
  read(20,4) opls(q), oplsnum(q), residue1(q), atom1(q), charge1(q) 
end do 



11 format (A16)  
12 format (A18)  
13 format (A9,5X,I1)
14 format (A10)   
15 format (A48) 
16 format (A9) 
17 format (A12)  
18 format (I6,1X,I6,1X,I1)   
19 format (I6,1X,I6,1X,I6,1X,I1) 
2 format (I6,X,I6,X,I6,X,I6)
10 format (A13) 
3 format (I6,1X,I6,1X,I6,1X,I6,1X,I1)  


read(30,11) header
read(30,12) comment
read(30,*) molecule, nrexcl
 
read(30,*)
read(30,14) atoms
read(30,15) titles
      
     allocate(nr(natoms))
     allocate(atmtype(natoms))
     allocate(resnr(natoms))
     allocate(residue(natoms))
     allocate(atom(natoms))
     allocate(cgnr(natoms))
     allocate(charge(natoms))
     allocate(mass(natoms))

do p = 1, natoms
  read(30,format1)nr(p),atmtype(p),resnr(p),residue(p),atom(p),cgnr(p),charge(p),mass(p)
end do 

read(30,*)
read(30,16)bonds
read(30,17)titlebond
      
     allocate(bondi(nbonds))
     allocate(bondj(nbonds))
     allocate(funcb(nbonds))

do i = 1, nbonds
  read(30,*) bondi(i), bondj(i), funcb(i)
end do

read(30,*)
read(30,14)angles
read(30,13)titleangles

     allocate(angi(nangles))
     allocate(angj(nangles))
     allocate(angl(nangles))
     allocate(funca(nangles))

do j = 1, nangles
  read(30,*) angi(j), angj(j), angl(j), funca(j)
end do

read(30,*)
read(30,10)dihedrals
read(30,12)titledihedrals

     allocate(dihedi(ndihedrals))
     allocate(dihedj(ndihedrals))
     allocate(dihedk(ndihedrals))
     allocate(dihedl(ndihedrals))
     allocate(funcd(ndihedrals))

do j = 1, ndihedrals
  read(30,*) dihedi(j), dihedj(j), dihedk(j), dihedl(j), funcd(j)
end do

read(30,*)
read(30,10)impropers
read(30,12)titleimpropers

     allocate(impropi(nimpropers))
     allocate(impropj(nimpropers))
     allocate(impropk(nimpropers))
     allocate(impropl(nimpropers))
     allocate(funci(nimpropers))

do j = 1, nimpropers
  read(30,*)  impropi(j), impropj(j), impropk(j), impropl(j), funci(j)  
end do

!!!!!!!Writing output!!!!!!!!!!!######################

write(40,11) header
write(40,12) comment
write(40,13) adjustl("GO       "), nrexcl
 
write(40,*)
!!!!###################atoms
write(40,14) atoms
write(40,15) titles


count=0
sum1=0.0d0
!!!!!!!!!!!!!!!!!!!change atomtypes to OPLS and charge
do p = 1, natoms 
  do q =1, oplstypes
    if (atom(p) .eq. atom1(q)) then     !!!!!atom names in OPLS file and .top file must match
       count=count+1
     write(40,format2)nr(p), opls(q),oplsnum(q),resnr(p),residue1(q),atom(p),cgnr(p),charge1(q),mass(p)
       !write(*,format1)nr(p),atmtype(p),resnr(p),residue(p),atom(p),cgnr(p),charge(p),mass(p)
       sum1=sum1+charge1(q)
    end if
   end do  
end do 


!!!!!!!!!!#############bonds
write(40,*)
write(40,16)bonds
write(40,17)titlebond
!
do i = 1, nbonds
write(40,18) bondi(i), bondj(i), funcb(i)
end do
!
!!!!################angles
write(40,*)
write(40,14)angles
write(40,13)titleangles
!
do j = 1, nangles
write(40,19) angi(j), angj(j), angl(j), funca(j)
end do
!count=0
!
!################proper dihedrals
func=3   !!!function for dihedrals 

write(40,*)
write(40,10)dihedrals
write(40,12)titledihedrals
!
do j = 1, ndihedrals
  write(40,3) dihedi(j), dihedj(j), dihedk(j), dihedl(j), func
end do
!################Impropers
write(40,*)
write(40,10)dihedrals
write(40,12)titleimpropers
!
do j = 1, nimpropers
  write(40,3)  impropi(j), impropj(j), impropk(j), impropl(j), func  
end do

!!!!###############1-4 Pairs
write(40,*)
write(40,*)"[ pairs ]"
write(40,17)titlebond

funcp=1  !!!!function for 1-4 pairs

do j = 1, ndihedrals
  write(40,18) dihedi(j), dihedl(j), funcp
end do

1 format (F10.6)
Write(*,1) sum1     !##### print total sum of charges 


close(30)

end program gromacs_top
