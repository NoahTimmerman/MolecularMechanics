!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  This is he assign module. In this module multiple subroutines and 
!!  functions are used. The first subroutine: 'AssignBonds', determines 
!!  which atoms are bonded to eachother, using the the function 
!!  'CalculateBondLength'. The bonds are then stored in a bondingarray
!!  by using the subroutine 'MakeBondList'. The next subroutine:
!!  'AssignAngles', calculate what the angle is of three bonded atoms.
!!  This subroutine also uses the function 'CalculateBondLength' and
!!  then store the angles in an anglearray, by using the subroutine:
!!  'MakeAngleList'. The next subroutine: 'AssignTorsional', calculates
!!  the torsional angles, by calculating the crossproduct using the 
!!  subroutine: 'CrossProduct'. At last the torsional angles are stored
!!  in an arrayb by using the subroutine: 'MakeTorsionalAngle'.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module AssignMod

    USE MoleculeMod
    implicit none

    PRIVATE
    PUBLIC :: Bond, AssignBonds, AssignAngles, AssignTorsional, CalculateBondLength

    type Bond
        Character(2) :: Elements
        Real*8       :: BondLength
        integer      :: FirstAtom, SecondAtom
    end type 

    contains 

real*8 function CalculateBondLength(FirstAtom, SecondAtom, MoleculeData) 
    type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
    integer, INTENT(IN)                    :: FirstAtom, SecondAtom
    real*8                                 :: Deltax, Deltay, Deltaz
      
        Deltax = MoleculeData(FirstAtom)%x - MoleculeData(SecondAtom)%x
        Deltay = MoleculeData(FirstAtom)%y - MoleculeData(SecondAtom)%y
        Deltaz = MoleculeData(FirstAtom)%z - MoleculeData(SecondAtom)%z

        CalculateBondLength = sqrt(Deltax**2 + Deltay**2 + Deltaz**2)

end function CalculateBondLength


subroutine AssignBonds(Bonding, MoleculeData, Variables, BondingArray)
    type(Bond), INTENT(INOUT)              :: Bonding
    type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
    type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
    type(Parameters), INTENT(IN)           :: Variables
    integer, allocatable                   :: AlreadyBonded(:,:)
    integer                                :: FirstAtom, SecondAtom

    ALLOCATE(AlreadyBonded(size(MoleculeData), size(MoleculeData)))
    AlreadyBonded = 0

    do FirstAtom = 1, size(MoleculeData)
        do SecondAtom = 1, size(MoleculeData)
            if(FirstAtom /= SecondAtom) then
                if (MoleculeData(FirstAtom)%Element == 'C' .and. MoleculeData(SecondAtom)%Element == 'C' .and.   &
                    abs(CalculateBondLength(FirstAtom, SecondAtom, MoleculeData) - Variables%CCBondLength) < 0.2 .and. &
                    AlreadyBonded(FirstAtom, SecondAtom) == 0) then
                            
                    Bonding%Elements = 'CC'
                    Bonding%BondLength = CalculateBondLength(FirstAtom, SecondAtom, MoleculeData)
                    Bonding%FirstAtom = FirstAtom
                    Bonding%SecondAtom = SecondAtom
                    AlreadyBonded(FirstAtom, SecondAtom) = 1
                    AlreadyBonded(SecondAtom, FirstAtom) = 1

                    call MakeBondList(BondingArray, Bonding)
                            
                elseif (MoleculeData(FirstAtom)%Element == 'C' .and. MoleculeData(SecondAtom)%Element == 'H' .and.  &
                        abs(CalculateBondLength(FirstAtom, SecondAtom, MoleculeData) - Variables%CHBondLength) < 0.2) then
                            
                    Bonding%Elements = 'CH'
                    Bonding%BondLength = CalculateBondLength(FirstAtom, SecondAtom, MoleculeData)
                    Bonding%FirstAtom = FirstAtom
                    Bonding%SecondAtom = SecondAtom

                    call MakeBondList(BondingArray, Bonding)
                endif
            endif
        enddo
    enddo
end subroutine AssignBonds


Subroutine MakeBondList(BondingArray, Bonding)
    type(Bond), ALLOCATABLE, INTENT(INOUT) :: BondingArray(:)
    type(Bond), INTENT(IN)                 :: Bonding
    type(Bond), ALLOCATABLE                :: DummyBondingArray(:)
    integer                                :: i, j

    if (allocated(BondingArray)) then
        j = size(BondingArray)
        ALLOCATE(DummyBondingArray(j+1))

        do i = 1,j
            DummyBondingArray(i) = BondingArray(i)
        enddo
        
        DummyBondingArray(j+1) = Bonding
        DEALLOCATE(BondingArray)
        call move_alloc(DummyBondingArray, BondingArray)
        
    else
        ALLOCATE(BondingArray(1))
        BondingArray(1) = Bonding
    endif

end subroutine MakeBondList

Subroutine AssignAngles(MoleculeData, BondingArray, Angle, AngleArray)
    type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
    type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
    real*8, INTENT(INOUT), ALLOCATABLE     :: AngleArray(:)
    integer                                :: i,j,k
    real*8, INTENT(INOUT)                  :: Angle

    do i = 1,size(MoleculeData)
        if (MoleculeData(i)%element == 'C') then
            do j = 1,size(BondingArray)
                do k = j,size(BondingArray)
                    if (j /= k .and. BondingArray(j)%FirstAtom == i .and. BondingArray(k)%FirstAtom == i) then   
                        Angle = -((CalculateBondLength(BondingArray(j)%SecondAtom, BondingArray(k)%SecondAtom, MoleculeData))**2 - &
                        (BondingArray(j)%Bondlength)**2 - (BondingArray(k)%Bondlength)**2) &
                        / (2*(BondingArray(j)%Bondlength) * (BondingArray(k)%Bondlength))
                        Call MakeAngleList(Angle, AngleArray)
                                           
                    elseif (j /= k .and. BondingArray(j)%SecondAtom == i .and. BondingArray(k)%FirstAtom == i) then
                        Angle = -((CalculateBondLength(BondingArray(j)%SecondAtom, BondingArray(k)%SecondAtom, MoleculeData))**2 - &
                                ((BondingArray(j)%Bondlength)**2 - (BondingArray(k)%Bondlength)**2) / &
                                (2*(BondingArray(j)%Bondlength) * (BondingArray(k)%Bondlength)))
                        Call MakeAngleList(Angle, AngleArray)
                    endif
                enddo
            enddo
        endif
    enddo
end subroutine AssignAngles

Subroutine MakeAngleList(Angle, AngleArray)
    real*8, ALLOCATABLE, INTENT(INOUT) :: AngleArray(:)
    real*8, INTENT(IN)                 :: Angle
    real*8, ALLOCATABLE                :: DummyAngleArray(:)
    integer                            :: i, j

    if (allocated(AngleArray)) then
        j = size(AngleArray)
        ALLOCATE(DummyAngleArray(j+1))

        do i = 1,j
            DummyAngleArray(i) = AngleArray(i)
        enddo
        
        DummyAngleArray(j+1) = Angle
        DEALLOCATE(AngleArray)
        call move_alloc(DummyAngleArray, AngleArray)
        
    else
        ALLOCATE(AngleArray(1))
        AngleArray(1) = Angle
    endif

end subroutine MakeAngleList

subroutine AssignTorsional(BondingArray, MoleculeData)
    type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
    type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
    integer                                :: i,j,k
    real*8, DIMENSION(3)                   :: VectorA, VectorB, VectorC, x, y, n1, n2, m
    real*8, ALLOCATABLE                    :: TorsionAngles(:), phi(:)

    do i = 1,size(BondingArray)
        if (BondingArray(i)%elements == 'CC') then
            do j = 1,size(BondingArray)
                do k = j,size(BondingArray)
                    if (j /= k .and. BondingArray(j)%FirstAtom == bondingArray(i)%FirstAtom &
                    .and. BondingArray(j)%SecondAtom /= bondingArray(i)%SecondAtom &
                    .and. BondingArray(k)%FirstAtom == BondingArray(i)%SecondAtom & 
                    .and. BondingArray(k)%SecondAtom /= BondingArray(i)%FirstAtom)then 

                        VectorA(1) = MoleculeData(BondingArray(j)%SecondAtom)%x - MoleculeData(BondingArray(i)%FirstAtom)%x
                        VectorB(1) = MoleculeData(BondingArray(i)%FirstAtom)%x - MoleculeData(BondingArray(i)%SecondAtom)%x
                        VectorC(1)  = MoleculeData(BondingArray(k)%SecondAtom)%x - MoleculeData(BondingArray(i)%SecondAtom)%x

                        VectorA(2) = MoleculeData(BondingArray(j)%SecondAtom)%y - MoleculeData(BondingArray(i)%FirstAtom)%y
                        VectorB(2) = MoleculeData(BondingArray(i)%FirstAtom)%y - MoleculeData(BondingArray(i)%SecondAtom)%y
                        VectorC(2)  = MoleculeData(BondingArray(k)%SecondAtom)%y - MoleculeData(BondingArray(i)%SecondAtom)%y

                        VectorA(3) = MoleculeData(BondingArray(j)%SecondAtom)%z - MoleculeData(BondingArray(i)%FirstAtom)%z
                        VectorB(3) = MoleculeData(BondingArray(i)%FirstAtom)%z - MoleculeData(BondingArray(i)%SecondAtom)%z
                        VectorC(3)  = MoleculeData(BondingArray(k)%SecondAtom)%z - MoleculeData(BondingArray(i)%SecondAtom)%z

                        n1 = CrossProduct(VectorA, VectorB)
                        n2 = CrossProduct(VectorB, VectorC)
                        m = CrossProduct(n1, VectorB)
                        x = DOT_PRODUCT(n1, n2)
                        y = DOT_PRODUCT(m, n2)
                        Phi = atan2(y,x)
                        call MakeTorsionAngle(Phi, TorsionAngles)

                    elseif(j /= k .and. BondingArray(j)%SecondAtom == bondingArray(i)%FirstAtom &
                    .and. BondingArray(j)%SecondAtom /= bondingArray(i)%SecondAtom &
                    .and. BondingArray(k)%FirstAtom == BondingArray(i)%SecondAtom & 
                    .and. BondingArray(k)%SecondAtom /= BondingArray(i)%FirstAtom)then   

                        VectorA(1) = MoleculeData(BondingArray(j)%SecondAtom)%x - MoleculeData(BondingArray(i)%FirstAtom)%x
                        VectorA(2) = MoleculeData(BondingArray(i)%FirstAtom)%x - MoleculeData(BondingArray(i)%SecondAtom)%x
                        VectorA(3)  = MoleculeData(BondingArray(k)%SecondAtom)%x - MoleculeData(BondingArray(i)%SecondAtom)%x

                        VectorB(1) = MoleculeData(BondingArray(j)%SecondAtom)%y - MoleculeData(BondingArray(i)%FirstAtom)%y
                        VectorB(2) = MoleculeData(BondingArray(i)%FirstAtom)%y - MoleculeData(BondingArray(i)%SecondAtom)%y
                        VectorB(2)  = MoleculeData(BondingArray(k)%SecondAtom)%y - MoleculeData(BondingArray(i)%SecondAtom)%y

                        VectorC(3) = MoleculeData(BondingArray(j)%SecondAtom)%z - MoleculeData(BondingArray(i)%FirstAtom)%z
                        VectorC(3) = MoleculeData(BondingArray(i)%FirstAtom)%z - MoleculeData(BondingArray(i)%SecondAtom)%z
                        VectorC(3)  = MoleculeData(BondingArray(k)%SecondAtom)%z - MoleculeData(BondingArray(i)%SecondAtom)%z

                        n1 = CrossProduct(VectorA, VectorB)
                        n2 = CrossProduct(VectorB, VectorC)
                        m = CrossProduct(n1, VectorB)
                        x = DOT_PRODUCT(n1, n2)
                        y = DOT_PRODUCT(m, n2)
                        Phi = atan2(y, x)
                        call MakeTorsionAngle(Phi, TorsionAngles)

                    endif   
                enddo
            enddo
        endif
    enddo
end subroutine AssignTorsional

function CrossProduct(x,y)
    real*8, DIMENSION(3) :: CrossProduct, x, y

    CrossProduct(1) = (x(2)*y(3)) - (x(3)*y(2))
    CrossProduct(2) = (x(1)*y(3)) - (x(3)*y(1))
    CrossProduct(3) = (x(1)*y(2)) - (x(2)*y(1))
end function CrossProduct

Subroutine MakeTorsionAngle(phi, TorsionAngles)
    real*8, ALLOCATABLE, INTENT(INOUT) :: phi(:), TorsionAngles(:)
    real*8, ALLOCATABLE                :: DummyAngleArray(:)
    integer                            :: i, j

    if (allocated(TorsionAngles)) then
        j = size(TorsionAngles)
        ALLOCATE(DummyAngleArray(j+1))

        do i = 1,j
            DummyAngleArray(i) = TorsionAngles(i)
        enddo
        
        DummyAngleArray(j+1) = phi(j+1)
        DEALLOCATE(TorsionAngles)
        call move_alloc(DummyAngleArray, TorsionAngles)
        
    else
        ALLOCATE(TorsionAngles(1))
        TorsionAngles(1) = phi(1)
    endif

end subroutine MakeTorsionAngle

end module 