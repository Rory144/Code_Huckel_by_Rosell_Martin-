Hückel Polyenes Homework 

# Project Overview

In this project, two different types of polymerization are carried out: alternation of bonds and alternation of atoms, taking into account two different types of 1D polyene chain: linear and cyclic. The results are generalized for an alternation of atoms to more than two types of atoms in the chain (assuming a regular periodic distribution, given as input. 


# Installation
1. Clone this repository: git clone https://github.com/Rory144/Code_for_Huckel_by_Rosell_Martin.git
2. Navigate to the project directory: cd Code_for_Huckel_by_Rosell_Martin
3. Install dependencies: pip install -r requirements.txt

# Usage

This project presents different codes: 

1. data.txt: This file contains all the input parameters required to specify before executing the program, called by the main program.  

2. main_programe.f90 -> This is the main program that use the module and calls the two subroutines. Also, it defines parameters to construct the hamiltonian matrix. 

3. mod.f90 -> This module contains two subroutine to calculate the diagonalized matrix, the eigenvalues, the GAP and Fermi Level.

4. Results directory -> Contains the files with the results for different studied cases. 

5. makefile: this file contains all the necessary computational instructions. To execute the code, simply type make in the terminal.

6. main_program.exe: the executable. To execute the code just you should: ./main_program.exe after write make.

# License
This project is licensed under the MIT License - see the LICENSE.md file for details.

Remark : the file generated when the program is run must be deleted if you want to run the same calculation with the same name.  
# Code_Huckel_by_Rosell_Martin-
