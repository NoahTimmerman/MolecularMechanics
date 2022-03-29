!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  This is the minimisation module. In this module two subroutines are
!!  used to minimise the energy. The subroutine: 'MovingMolecule' moves
!!  the x,y,z coordinates of 1 random atom of the total molecule and
!!  stores this information in a new molecule array. The second subroutine
!!  'MinimisingEnergy', calls this subroutine and checks if the TotalEnergy
!!  of the new coordinated molecule is lower than the initial energy. If so,
!!  the initial energy is set to this new energy en the inital molecule
!!  coordinates are updated to the new coordinated molecule. This is repeated
!!  till the do loop reaches the exit statement of 500.000 rejections.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module MinimiseMod

    USE MoleculeMod
    USE AssignMod
    USE EnergyCalculationMod
    
    implicit none

    PRIVATE
    PUBLIC :: MinimisingEnergy

contains

subroutine MovingMolecule(Variables, MoleculeMoved)
    type(Atom), INTENT(INOUT), ALLOCATABLE  :: MoleculeMoved(:)
    type(Parameters), INTENT(IN)            :: Variables
    real*8, DIMENSION(3)                    :: NewCoord, Rand
    real*8                                  :: Random
    integer                                 :: i, j

    call random_number(NewCoord(1))
    call random_number(NewCoord(2))
    call random_number(NewCoord(3))

    do i = 1,3
        call random_number(Rand(i))
        if(Rand(i) >= 0.5)then
            NewCoord(i) = -NewCoord(i)
        endif
    enddo

    do 
        call random_number(Random)
        j = int(Random * size(MoleculeMoved))
        if (j > 0) exit
    enddo
  
    MoleculeMoved(j)%x = MoleculeMoved(j)%x + Variables%r*NewCoord(1)
    MoleculeMoved(j)%y = MoleculeMoved(j)%y + Variables%r*NewCoord(2)
    MoleculeMoved(j)%z = MoleculeMoved(j)%z + Variables%r*NewCoord(3)

end subroutine

subroutine MinimisingEnergy(Variables, Molecule, BondingArray, AngleArray, TorsionalAngles, InitialEnergy)
    type(Atom), INTENT(INOUT), ALLOCATABLE  :: Molecule(:)
    type(Atom), ALLOCATABLE                 :: MoleculeMoved(:)
    type(Parameters), INTENT(INOUT)         :: Variables
    type(Bond), INTENT(INOUT), ALLOCATABLE  :: BondingArray(:)
    real*8, INTENT(IN), ALLOCATABLE         :: AngleArray(:)
    real*8, INTENT(IN), ALLOCATABLE         :: TorsionalAngles(:)
    real*8, INTENT(INOUT)                   :: InitialEnergy
    real*8                                  :: NewEnergy, DeltaEnergy, Pa, q
    integer                                 :: AcceptedEnergies, RejectedEnergies

    ALLOCATE(MoleculeMoved(size(Molecule)))
    MoleculeMoved = Molecule

    AcceptedEnergies = 0
    RejectedEnergies = 0
print *, InitialEnergy
    do
        call MovingMolecule(Variables, MoleculeMoved)
        NewEnergy = TotalEnergy(Variables, BondingArray, AngleArray, TorsionalAngles, MoleculeMoved)
        DeltaEnergy = NewEnergy - InitialEnergy
        if(DeltaEnergy < 0)then 
            Molecule = MoleculeMoved
            InitialEnergy = NewEnergy
            AcceptedEnergies = AcceptedEnergies + 1
        else
            Pa = exp((-DeltaEnergy)/(Variables%Kb * Variables%Temp))
            call random_number(q)
            if(q < Pa)then 
                Molecule = MoleculeMoved
                InitialEnergy = NewEnergy
                AcceptedEnergies = AcceptedEnergies + 1
            else
                RejectedEnergies = RejectedEnergies + 1
            endif
        endif
        if(RejectedEnergies >= 1000000)then
            exit
        endif  
    enddo
print *, InitialEnergy, AcceptedEnergies

end subroutine MinimisingEnergy

end module MinimiseMod