Module CalculationMod

    USE AtomDataMod
    implicit none

    PRIVATE
    PUBLIC :: Bond, AssignBonds, CalculateBondLength, MakeBondList, CalculateBending

    type Bond
        Character(2) :: Elements
        Real*8 :: BondLength
        integer :: FirstAtom, SecondAtom
    end type 

    contains 

    real*8 function CalculateBondLength(FirstAtom, SecondAtom, MoleculeData) 
        type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
        integer, INTENT(IN) :: FirstAtom, SecondAtom
        real*8 :: Deltax, Deltay, Deltaz
      
            Deltax = MoleculeData(FirstAtom)%x - MoleculeData(SecondAtom)%x
            Deltay = MoleculeData(FirstAtom)%y - MoleculeData(SecondAtom)%y
            Deltaz = MoleculeData(FirstAtom)%z - MoleculeData(SecondAtom)%z

            CalculateBondLength = sqrt(Deltax**2 + Deltay**2 + Deltaz**2)

    end function CalculateBondLength


    subroutine AssignBonds(Bonding, MoleculeData, Variables, BondingArray)
        type(Bond), INTENT(INOUT) :: Bonding
        type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
        type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
        type(Parameters), INTENT(IN) :: Variables
        integer, allocatable :: AlreadyBonded(:,:)
        integer :: FirstAtom, SecondAtom

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
        type(Bond), INTENT(IN) :: Bonding
        type(Bond), ALLOCATABLE :: DummyBondingArray(:)
        integer :: i, j

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

    Subroutine CalculateBending(MoleculeData, BondingArray, Variables)
        type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
        type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
        type(Parameters), INTENT(IN) :: Variables
        integer :: FirstAtom,SecondAtom,k

        do FirstAtom = 1,size(MoleculeData)  
            do SecondAtom = FirstAtom,size(MoleculeData)
                if(SecondAtom /= FirstAtom .and. MoleculeData(FirstAtom)%Element == 'C' & 
                .and. MoleculeData(SecondAtom)%Element == 'C'  &
                .and. abs(CalculateBondLength(FirstAtom, SecondAtom, MoleculeData) - Variables%CCBondLength) < 0.2)then
                    do k = SecondAtom,size(BondingArray)
                        if(k /= SecondAtom .and. k /= FirstAtom .and. BondingArray(k)%FirstAtom == FirstAtom) then
                            print *, Moleculedata(FirstAtom)%element, MoleculeData(SecondAtom)%Element, MoleculeData(k)%Element
                        endif
                    enddo
                elseif(SecondAtom /= FirstAtom .and. MoleculeData(FirstAtom)%Element == 'C' &
                .and. MoleculeData(SecondAtom)%Element == 'H' .and.  &
                abs(CalculateBondLength(FirstAtom, SecondAtom, MoleculeData) - Variables%CHBondLength) < 0.2) then
                    do k = SecondAtom,size(BondingArray)
                        if(k /= SecondAtom .and. k /= FirstAtom .and. BondingArray(k)%FirstAtom == FirstAtom) then
                            print *, Moleculedata(FirstAtom)%element,MoleculeData(SecondAtom)%Element,MoleculeData(k)%Element
                        endif
                    enddo
                        
                endif
            enddo  
        enddo
    end subroutine



   ! subroutine Planes(BondingArray,MoleculeData)
       ! type(Bond), INTENT(INOUT), ALLOCATABLE :: BondingArray(:)
       ! type(Atom), allocatable, INTENT(INOUT) :: MoleculeData(:)
        !integer :: i,j,k
        !real*8 :: VectorX1

        !do i = 1,size(MoleculeData)
         !   if (MoleculeData(i)%element == 'C') then
          !     do j = 1,size(BondingArray)
           !       do k = j,size(BondingArray)
             !           if (j /= k .and. BondingArray(j)%FirstAtom == i .and. BondingArray(k)%FirstAtom == i ) then
            !               VectorX1 = MoleculeData(BondingArray(j)%SecondAtom)%x - MoleculeData(i)%x
             !               print *, VectorX1
              !          endif
               !     enddo
                !enddo
            !endif
        !enddo

    !end subroutine 
        
end module 