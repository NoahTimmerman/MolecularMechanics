module AtomData

    implicit none

    PRIVATE
    PUBLIC ::  Atom, Parameters, AssignParameters, AtomReading 

    type Atom
        character(1) :: Element
        real*8 :: x,y,z
    end type Atom

    type Parameters
        real*8 :: CCBondLength, CHBondLength, CCBondStretch, CHBondStretch, ForceConstantAngle, EquiAngle, V1, n, Y, ChargeH, &
        ChargeC, Coulombs, Avogadro, SigmaH, SigmaC, WellDepthC, WellDepthH
         
    end type 

    contains 

    subroutine AtomReading(filename, MoleculeData, AtomNumbers)
        character(len=*), INTENT(IN):: filename
        integer, INTENT(INOUT) :: AtomNumbers
        type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
        integer :: i, iu
        
        open(newunit = iu, file = filename, action = 'read') 
        read(iu,*) AtomNumbers                     
        
        allocate(MoleculeData(AtomNumbers))
    
            do i = 1,AtomNumbers
                read(iu,*) MoleculeData(i)%Element, MoleculeData(i)%x, MoleculeData(i)%y, MoleculeData(i)%z
            enddo

        close(iu)
      
    end subroutine AtomReading

    subroutine AssignParameters(Variables)
        type(Parameters), INTENT(INOUT) :: Variables

        Variables%CCBondLength  = 1.507
        Variables%CHBondLength  = 1.094
        Variables%CCBondStretch = 317
        Variables%CHBondStretch = 300
        Variables%ForceConstantAngle = 339
        Variables%EquiAngle = 2.3
        Variables%V1 = 1.50
        Variables%n = 1.0
        Variables%Y = 0
        Variables%ChargeH = 1.602176/(10.0**19)
        Variables%ChargeC = 1.602176/(10.0**19)
        Variables%Coulombs = 8987551787.0
        Variables%Avogadro = 6.022141*(10.0**23)
        Variables%SigmaH = 2
        Variables%SigmaC = 2
        Variables%WellDepthH = 2
        Variables%WellDepthC  = 2

    end subroutine 

end module 