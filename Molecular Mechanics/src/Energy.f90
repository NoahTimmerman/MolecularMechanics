module EnergyCalculation 

    USE AtomDataMod
    USE CalculationMod
    implicit none

    PRIVATE
    PUBLIC :: StretchEnergy
    
    contains

    !subroutine TotalEnergy(Variables, MoleculeData, BondingArray)
     !   type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
      !  type(Parameters), INTENT(IN) :: Variables
      !  type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
       ! real*8 :: totalen
        
       

    !end subroutine
    

    real*8 function StretchEnergy(Variables, BondingArray)
        type(Parameters), INTENT(INOUT) :: Variables
        type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
        integer :: i

        StretchEnergy = 0.0

        do i=1,size(BondingArray)
            if(BondingArray(i)%Elements == 'CH') then
                StretchEnergy = StretchEnergy + Variables%CHBondStretch * ((BondingArray(i)%BondLength - Variables%CHBondLength)**2)
            else
                StretchEnergy = StretchEnergy + Variables%CCBondStretch * ((BondingArray(i)%BondLength - Variables%CCBondLength)**2)
            endif
        enddo

    end function StretchEnergy



end module 