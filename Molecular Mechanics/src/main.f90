!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Author          : Noah Timmerman
!!  StudentNumber   : 14052687
!!  University      : University of Amsterdam
!!  Date            : 29-03-2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  
!!  Project         : Molecular Mechanics
!!
!!  This code calculates the minimised energy of a coumpound which only 
!!  exists of Carbons and hydrogens, and must be saturated. This code 
!!  exists of 4 modules: MoleculeMod, AssignMod, EnergyCalculationmod,
!!  MinimiseMod, and this main program.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program MolecularMechanics

    USE MoleculeMod
    USE AssignMod
    USE EnergyCalculationMod
    USE MinimiseMod

    implicit none

    type (Atom), ALLOCATABLE :: Molecule(:)
    type (Parameters)        :: Variables
    type (Bond)              :: Bonding
    type (Bond), ALLOCATABLE :: BondingArray(:)
    integer                  :: AtomNumbers
    real *8, ALLOCATABLE     :: AngleArray(:)
    real*8, ALLOCATABLE      :: TorsionalAngles(:)
    real*8                   :: Angle
    real*8                   :: InitialEnergy
    
   
    call AtomReading('c4h10.xyz', Molecule, AtomNumbers)
    call AssignParameters(Variables)
    call AssignBonds(Bonding, Molecule, Variables, BondingArray)
    call AssignAngles(Molecule, BondingArray, Angle, AngleArray)
    call AssignTorsional(BondingArray, Molecule)
 
    InitialEnergy = TotalEnergy(Variables, BondingArray, AngleArray, TorsionalAngles, Molecule)

    call MinimisingEnergy(Variables, Molecule, BondingArray, AngleArray, TorsionalAngles, InitialEnergy)
    print *, 'Minimised energy is:', InitialEnergy

end program 
