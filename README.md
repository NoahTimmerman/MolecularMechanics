# MolecularMechanics
Final project of Scientific Computing & Programming

This code calculates the minimised energy of a coumpound which only exists of Carbons and hydrogens, and must be saturated. This code 
exists of 4 modules: MoleculeMod, AssignMod, EnergyCalculationmod, MinimiseMod, and the main program. 

In the main program the first subroutine called is the subroutine AtomReading, from the molecule module. In this subroutine the text file is 
red and the coordinations and atoms are stored in type 'atom'. In this module the parameters are saved to the type 'parameters'. 

The next subroutine called by the main is the assignbonds from the Assign module. In the Assignbond subroutine the atoms which are bonded in 
the molecule are stored in a 'bondingarray'. This is done by using a type 'bond'. Using the function 'CalculateBondLength' the bond length 
between two atoms is calculated. If this value is equal or lower than the treshold (there are two tresholds used, CH and CC bond) then the two 
atoms are indeed bonded. If the atoms are bonded the subroutine 'makebondlist' is called and the the two atoms are stored in an array. 

The next subroutine in the Assign module is the Assignangles. The first thing we need to know is which pairs of three atoms we have. this done
by using 3 do loops and 2 if statements. If three atoms are bonded in a row, than the middle atom must be a carbon. So, the i must b a 'C'. 
Than j and k must be a different atom and must be both bonded to 'C'. When a pair of three is founded this is stored in an array using the 
'MakeAngleList' subroutine. 

The last thing needs to be determined is the pairs of 4 atoms. This is done in the subroutine 'AssignTorsional'. If 4 atoms are bonded, than the 
two atoms in the middle must be 'C-C'. So again 3 do loops, i,j,and k, are used. If in the bonded array place i is 'CC' then the two other loops 
will be executed. Than j and k must be different atoms again and using the two different if statements the pairs of 4 are made. Then by making 
vectors the crossproduct is calculated using the function 'crossproduct'. And the 'MakeTorsionAngle' makes an array of the pairs of 4 atoms.

Next using these arrays the energy can be calculated. This is done by using the module EnergyCalculationMod. The first function: 
'TotalEnergyfunction' takes as input the 4 other functions. The first function taken as input is 'StretchEnergy' function, this function 
calculates the stretchenergy using the BondingArray. The second function: 'BendingEnergy'calculates the bending energy which uses the angle array.
The third energy function: 'TorsionalEnergy' calculates the torsionalenergy using the torsional angles array, and the last function: 
'NonBondedEnergy', calculates the non-bonded energy, using the electrostatic and van der Waal's interactions. With these 4 functions the total
energy is then calculated using the function: 'TotalEnergy'.

The energy minimisation is done by making small difference in the coordination of the molecule. This module works as follows: first the accepted- 
and rejectedenergies are set to '0'. Next the subroutine 'MinimisingEnergy' calls the subroutine 'MovingMolecule'. In this subroutine a single 
atom on the molecule is changed in coordinations. this conformation is then saved as 'NewMolecule'. With this 'NewMolecule' the totalenergy is 
calculated. this newenergy is compared to the initial energy and if the new energy is lower the new conformation is accepted. Else, the chance
is calculated to accept it anyway, by using the formula: -DeltaEnergy)/(Variables%Kb * Variables%Temp. If the Pa is greaterthan a random called
number then the conformation is accepted. If not, the old molecule conformation will be changed again. This process is done till there are 500.000 
rejections. 
