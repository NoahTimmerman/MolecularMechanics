module AtomDataMod

    implicit none

    PRIVATE
    PUBLIC ::  AtomReading, Atom, AssignParameters, Parameters

    type Atom
        character(1) :: Element
        real*8 :: x,y,z
    end type Atom

    type Parameters
        real*8 :: CCBondLength, CHBondLength, CCBondStretch, CHBondStretch
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

    end subroutine 

end module  