module EnergyCalculation 

    USE AtomData
    USE CalculationMod
    implicit none

    PRIVATE
    PUBLIC :: StretchEnergy, BendingEnergy, TorsionalEnergy, NonBondedEnergy, TotalEnergy
    
    contains

real*8 function TotalEnergy(Variables, BondingArray, AngleArray, TorsionalAngles, MoleculeData)
    type(Parameters), INTENT(INOUT) :: Variables
    type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
    type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
    real*8, INTENT(IN), ALLOCATABLE :: AngleArray(:)
    real*8, INTENT(IN), ALLOCATABLE :: TorsionalAngles(:)
  
    TotalEnergy = StretchEnergy(Variables, BondingArray) + BendingEnergy(Variables, AngleArray) + &
    TorsionalEnergy(Variables, TorsionalAngles) + NonBondedEnergy(Variables, MoleculeData)

end function 

real*8 function StretchEnergy(Variables, BondingArray)
    type(Parameters), INTENT(INOUT) :: Variables
    type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
    integer :: i

    StretchEnergy = 0.0

    do i=1,size(BondingArray)
        if(BondingArray(i)%Elements == 'CH')then
            StretchEnergy = StretchEnergy + Variables%CHBondStretch * ((BondingArray(i)%BondLength - Variables%CHBondLength)**2)
        else
            StretchEnergy = StretchEnergy + Variables%CCBondStretch * ((BondingArray(i)%BondLength - Variables%CCBondLength)**2)
        endif
    enddo

end function StretchEnergy

real*8 function BendingEnergy(Variables, AngleArray)
    type(Parameters), INTENT(INOUT) :: Variables
    real*8, INTENT(IN), ALLOCATABLE :: AngleArray(:)
    integer :: i

    BendingEnergy = 0.0
    do i = 1,size(AngleArray)
        if (abs(acos(AngleArray(i)) - Variables%EquiAngle) < 1)then
        BendingEnergy = BendingEnergy + Variables%ForceConstantAngle * (((acos(AngleArray(i)) - Variables%EquiAngle))**2)
        endif
    enddo

end function BendingEnergy

real*8 function TorsionalEnergy(Variables, TorsionalAngles)
    type(Parameters), INTENT(INOUT) :: Variables
    real*8, INTENT(IN), ALLOCATABLE :: TorsionalAngles(:)
    integer :: i

    do i = 1,size(TorsionalAngles)
        TorsionalEnergy = TorsionalEnergy + (0.5 * Variables%V1 * (1 + Cos((Variables%n * TorsionalAngles(i)) - Variables%Y)))
    enddo
end function TorsionalEnergy

real*8 function NonBondedEnergy(Variables, MoleculeData)
    type(Parameters), INTENT(INOUT) :: Variables
    type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
    integer :: i,j

    do i = 1,size(MoleculeData)
        do j = 1,size(MoleculeData)
            if (i /= j)then
                if (MoleculeData(i)%Element == 'H' .and. MoleculeData(j)%Element =='H')then
                    if (j > i)then 
                        NonBondedEnergy = NonBondedEnergy + (Variables%Coulombs * Variables%Avogadro * (Variables%ChargeH**2)  / &
                        (CalculateBondLength(i,j,MoleculeData)/(10.0**10)))
                    endif
                    NonBondedEnergy = NonBondedEnergy + ( (4 * Variables%WellDepthH) * &
                    ((Variables%SigmaH / CalculateBondLength(i,j,MoleculeData)**12) - &
                    ((Variables%SigmaH / CalculateBondLength(i,j,MoleculeData)**6))) )
                
                elseif (MoleculeData(i)%Element == 'C' .and. MoleculeData(j)%Element == 'C')then
                    if (j > i)then 
                        NonBondedEnergy = NonBondedEnergy + (Variables%Coulombs * Variables%Avogadro * (Variables%ChargeH**2)  / &
                        (CalculateBondLength(i,j,MoleculeData)/(10.0**10)))
                    endif
                    NonBondedEnergy = NonBondedEnergy + ( (4 * Variables%WellDepthC) * &
                    ((Variables%SigmaC / CalculateBondLength(i,j,MoleculeData)**12) - &
                    ((Variables%SigmaC / CalculateBondLength(i,j,MoleculeData)**6))) )
                else
                    if (j > i)then 
                        NonBondedEnergy = NonBondedEnergy + (Variables%Coulombs * Variables%Avogadro * (Variables%ChargeH**2)  / &
                        (CalculateBondLength(i,j,MoleculeData)/(10.0**10)))
                    endif
                    NonBondedEnergy = NonBondedEnergy + ( (4 * ((Variables%WellDepthC/2) + (Variables%WellDepthH/2)) ) * &
                    (( ((Variables%SigmaC/2)+(Variables%SigmaH/2)) / CalculateBondLength(i,j,MoleculeData)**12) - &
                    (( ((Variables%SigmaC/2)+(Variables%SigmaH/2)) / CalculateBondLength(i,j,MoleculeData)**6))) )
                endif
            endif
        enddo
    enddo
end function NonBondedEnergy                      

end module 