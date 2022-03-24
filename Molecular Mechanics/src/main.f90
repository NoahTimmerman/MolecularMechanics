program MolecularMechanics

    USE AtomDataMod
    USE CalculationMod
    USE EnergyCalculation

    IMPLICIT NONE

    type (Atom), ALLOCATABLE :: MoleculeData(:)
    type (Parameters) :: variables
    type (Bond) :: Bonding
    type (Bond), ALLOCATABLE :: BondingArray(:)
    integer :: AtomNumbers
    
   
    call AtomReading('c4h10.xyz', MoleculeData, AtomNumbers)
    call AssignParameters(Variables)
    !print *, MoleculeData
    
    call AssignBonds(Bonding, MoleculeData, Variables, BondingArray)
    print *, BondingArray
    !call planes(BondingArray, MoleculeData)
    print *, StretchEnergy(Variables, BondingArray)
    call CalculateBending(MoleculeData, BondingArray, variables)

end program 