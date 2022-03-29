module EnergyCalculation 

    USE MoleculeMod
    USE CalculationMod
    implicit none

    PRIVATE
    PUBLIC :: StretchEnergy, BendingEnergy, TorsionalEnergy, NonBondedEnergy, TotalEnergy
    
    contains

real*8 function TotalEnergy(Variables, BondingArray, AngleArray, TorsionalAngles, Molecule)
    type(Parameters), INTENT(INOUT)        :: Variables
    type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
    type(Atom), allocatable, INTENT(INOUT) :: Molecule(:)
    real*8, INTENT(IN), ALLOCATABLE        :: AngleArray(:)
    real*8, INTENT(IN), ALLOCATABLE        :: TorsionalAngles(:)
  
    TotalEnergy = StretchEnergy(Variables, BondingArray) + BendingEnergy(Variables, AngleArray) + &
    TorsionalEnergy(Variables, TorsionalAngles) + NonBondedEnergy(Variables, Molecule)

end function 

real*8 function StretchEnergy(Variables, BondingArray)
    type(Parameters), INTENT(INOUT)        :: Variables
    type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
    integer                                :: i

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
    integer                         :: i

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
    integer                         :: i

    do i = 1,size(TorsionalAngles)
        TorsionalEnergy = TorsionalEnergy + (0.5 * Variables%V1 * (1 + Cos((Variables%n * TorsionalAngles(i)) - Variables%Y)))
    enddo
end function TorsionalEnergy

real*8 function NonBondedEnergy(Variables, Molecule)
    type(Parameters), INTENT(INOUT)        :: Variables
    type(Atom), allocatable, INTENT(INOUT) :: Molecule(:)
    integer                                :: i,j

    do i = 1,size(Molecule)
        do j = 1,size(Molecule)
            if (i /= j)then
                if (Molecule(i)%Element == 'H' .and. Molecule(j)%Element =='H')then
                    if (j > i)then 
                        NonBondedEnergy = NonBondedEnergy + (Variables%Coulombs * Variables%Avogadro * (Variables%ChargeH**2)  / &
                        (CalculateBondLength(i,j,Molecule)/(10.0**10)))
                    endif
                    NonBondedEnergy = NonBondedEnergy + ( (4 * Variables%WellDepthH) * &
                    ((Variables%SigmaH / CalculateBondLength(i,j,Molecule)**12) - &
                    ((Variables%SigmaH / CalculateBondLength(i,j,Molecule)**6))) )
                
                elseif (Molecule(i)%Element == 'C' .and. Molecule(j)%Element == 'C')then
                    if (j > i)then 
                        NonBondedEnergy = NonBondedEnergy + (Variables%Coulombs * Variables%Avogadro * (Variables%ChargeH**2)  / &
                        (CalculateBondLength(i,j,Molecule)/(10.0**10)))
                    endif
                    NonBondedEnergy = NonBondedEnergy + ( (4 * Variables%WellDepthC) * &
                    ((Variables%SigmaC / CalculateBondLength(i,j,Molecule)**12) - &
                    ((Variables%SigmaC / CalculateBondLength(i,j,Molecule)**6))) )
                else
                    if (j > i)then 
                        NonBondedEnergy = NonBondedEnergy + (Variables%Coulombs * Variables%Avogadro * (Variables%ChargeH**2)  / &
                        (CalculateBondLength(i,j,Molecule)/(10.0**10)))
                    endif
                    NonBondedEnergy = NonBondedEnergy + ( (4 * ((Variables%WellDepthC/2) + (Variables%WellDepthH/2)) ) * &
                    (( ((Variables%SigmaC/2)+(Variables%SigmaH/2)) / CalculateBondLength(i,j,Molecule)**12) - &
                    (( ((Variables%SigmaC/2)+(Variables%SigmaH/2)) / CalculateBondLength(i,j,Molecule)**6))) )
                endif
            endif
        enddo
    enddo
end function NonBondedEnergy                      

end module 