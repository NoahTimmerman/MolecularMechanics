!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  This is the molecule module. In this module two subroutines are used.
!!  The first subroutine reads the input text file and stores it in the
!!  type 'Atom' the first column is the element and the other 3 columns are
!!  the x, y and z coordinates. The second subroutine assigns the parameters 
!!  to the type 'parameters'. 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module MoleculeMod

    implicit none

    PRIVATE
    PUBLIC ::  Atom, Parameters, AssignParameters, AtomReading 

    type Atom
        character(1) :: Element
        real*8       :: x,y,z
    end type Atom

    type Parameters
        real*8 :: CCBondLength, CHBondLength, CCBondStretch, CHBondStretch, ForceConstantAngle, EquiAngle, V1, n, Y, ChargeH, &
        ChargeC, Coulombs, Avogadro, SigmaH, SigmaC, WellDepthC, WellDepthH, r, Kb, Temp
         
    end type 

    contains 

    subroutine AtomReading(filename, Molecule, AtomNumbers)
        character(len=*), INTENT(IN)           :: filename
        integer, INTENT(INOUT)                 :: AtomNumbers
        type(Atom), allocatable, INTENT(INOUT) :: Molecule(:)
        integer                                :: i, iu
        
        open(newunit = iu, file = filename, action = 'read') 
        read(iu,*) AtomNumbers                     
        
        allocate(Molecule(AtomNumbers))
    
            do i = 1,AtomNumbers
                read(iu,*) Molecule(i)%Element, Molecule(i)%x, Molecule(i)%y, Molecule(i)%z
            enddo

        close(iu)
      
    end subroutine AtomReading

    subroutine AssignParameters(Variables)
        type(Parameters), INTENT(INOUT) :: Variables

        Variables%CCBondLength  = 1.507
        Variables%CHBondLength  = 1.094
        Variables%CCBondStretch = 1.326*(10.0**6)
        Variables%CHBondStretch = 1.255*(10**6)
        Variables%ForceConstantAngle = 475720.8
        Variables%EquiAngle = 2.0
        Variables%V1 = 5857.6
        Variables%n = 2.0
        Variables%Y = 0
        Variables%ChargeH = 0.078 * (1.602176/(10.0**19))
        Variables%ChargeC = -0.344 * (1.602176/(10.0**19))
        Variables%Coulombs = 8987551787.0
        Variables%Avogadro = 6.022141*(10.0**23)
        Variables%SigmaH = 0.600
        Variables%SigmaC = 1.9080
        Variables%WellDepthH = 65.688
        Variables%WellDepthC  = 457.7296
        Variables%r = 1.0/(10.0**5)
        Variables%Kb = 1.38064852/(10.0**23)
        Variables%Temp = 300
        
    end subroutine 

end module  
