module Minimise

    use MoleculeMod
    use Calculation
    use EnergyCalculation
implicit none
    
contains

subroutine MovingMolecule(MoleculeData)
    type(Atom), INTENT(INOUT), ALLOCATABLE  :: MoleculeData(:)
    type(Atom), ALLOCATABLE                 :: MoleculeMoving(:)
    real*8, DIMENSION(3)                    :: NewCoord, Rand
    real*8                                  :: Random
    integer                                 :: i, j

    allocate(MoleculeMoving(size(MoleculeData)))
    MoleculeMoving = MoleculeData
    
    call random_number(NewCoord(1))
    call random_number(NewCoord(2))
    call random_number(NewCoord(3))

    do j = 1,3
        call random_number(Rand(j))
        if(Rand(j) >= 0.5)then
            NewCoord(j) = -NewCoord(j)
        endif
    enddo

    do 
        call random_number(Random)
        i = int(Random * size(MoleculeData))
        if (i > 0) exit
    enddo
    
    MoleculeMoving(i)%x = MoleculeMoving(i)%x + 0.0001*NewCoord(1)
    MoleculeMoving(i)%y = MoleculeMoving(i)%y + 0.0001*NewCoord(2)
    MoleculeMoving(i)%z = MoleculeMoving(i)%z + 0.0001*NewCoord(3)

end subroutine

end module Minimise