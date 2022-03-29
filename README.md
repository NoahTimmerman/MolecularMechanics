# MolecularMechanics
Final project of Scientific Computing & Programming

This code calculates the minimised energy of a coumpound which only 
exists of Carbons and hydrogens, and must be saturated. This code 
exists of 4 modules: MoleculeMod, AssignMod, EnergyCalculationmod,
MinimiseMod, and the main program. The minimisation is done by making
small difference in the coordination of the molecule. The delta Energy
is calculated and if the new coordinated molecule is lower in energy
this conformation is accepted. After 500.000 rejections the final 
conformation is found to be the lowest in energy. 
