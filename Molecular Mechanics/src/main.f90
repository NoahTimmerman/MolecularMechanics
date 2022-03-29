program MolecularMechanics

    USE MoleculeMod
    USE CalculationMod
    USE EnergyCalculation
    USE Minimise

    IMPLICIT NONE

    type (Atom), ALLOCATABLE :: Molecule(:)
    type (Parameters)        :: Variables
    type (Bond)              :: Bonding
    type (Bond), ALLOCATABLE :: BondingArray(:)
    integer                  :: AtomNumbers
    real *8, ALLOCATABLE     :: AngleArray(:)
    real*8, ALLOCATABLE      :: TorsionalAngles(:)
    real*8                   :: Angle
    real*8 :: InitialEnergy
    
   
    call AtomReading('c4h10.xyz', Molecule, AtomNumbers)
    call AssignParameters(Variables)
    call AssignBonds(Bonding, Molecule, Variables, BondingArray)
    call AssignAngles(Molecule, BondingArray, Angle, AngleArray)
    call AssignTorsional(BondingArray, Molecule)
 
    InitialEnergy = TotalEnergy(Variables, BondingArray, AngleArray, TorsionalAngles, Molecule)

    call MinimisingEnergy(Variables, Molecule, BondingArray, AngleArray, TorsionalAngles, InitialEnergy)

end program 