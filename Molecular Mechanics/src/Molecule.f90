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
        Variables%CCBondStretch = 317
        Variables%CHBondStretch = 300
        Variables%ForceConstantAngle = 150
        Variables%EquiAngle = 2
        Variables%V1 = 50
        Variables%n = 1.0
        Variables%Y = 0
        Variables%ChargeH = 1.602176/(10.0**19)
        Variables%ChargeC = 1.602176/(10.0**19)
        Variables%Coulombs = 8987551787.0
        Variables%Avogadro = 6.022141*(10.0**23)
        Variables%SigmaH = 1
        Variables%SigmaC = 1
        Variables%WellDepthH = 0.066
        Variables%WellDepthC  = 0.033
        Variables%r = 0.00000001
        Variables%Kb = 1.38064852/(10.0**23)
        Variables%Temp = 300

    end subroutine 

end module  