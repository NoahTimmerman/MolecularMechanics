program MolecularMechanics

    USE MoleculeMod
    USE Calculation
    USE EnergyCalculation
    USE Minimise

    IMPLICIT NONE

    type (Atom), ALLOCATABLE :: MoleculeData(:)
    type (Parameters) :: variables
    type (Bond) :: Bonding
    type (Bond), ALLOCATABLE :: BondingArray(:)
    integer :: AtomNumbers
    real *8, ALLOCATABLE :: AngleArray(:)
    real*8, ALLOCATABLE :: TorsionalAngles(:)
    real*8 :: Angle
    
   
    call AtomReading('c4h10.xyz', MoleculeData, AtomNumbers)
    call AssignParameters(Variables)
    call AssignBonds(Bonding, MoleculeData, Variables, BondingArray)
    call CalculateAngle(MoleculeData, BondingArray, Angle, AngleArray)
    call planes(BondingArray, MoleculeData)
    
    print *, NonBondedEnergy(Variables, MoleculeData)
    print *, TotalEnergy(Variables, BondingArray, AngleArray, TorsionalAngles, MoleculeData)

    call MovingMolecule(MoleculeData)

end program 