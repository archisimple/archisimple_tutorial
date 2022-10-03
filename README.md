# archisimple_tutorial
Repo with exemple and explanations how to run the model

List of parameters : https://docs.google.com/spreadsheets/d/1wg_lvLkCnzfpeOvSQWQK7hsWt8fl-IZrAEiSfNwA1yA/edit?usp=sharing


## Installing ArchiSimple
Before running ArchiSimple, the source code should be compiled from C++, so it could be run as an executable on your local machine. 
The code should compile with any c++11 compiler, e.g. for g++: MinGWa dn cygwin have been tested.

### 0. Install compiler (windows only)

For cygwin installation, follow the instruction here https://www.optimizationcore.com/coding/compile-cpp-files-with-g-cygwin-windows/ and then here https://www.optimizationcore.com/free-solutions/cygwin-linux-tools-windows-install-dependency-check/

Check gpp installation : 

    gpp --version

### 1. Go to directory
Open the terminal or command invite (windows):
Go to the needed directory where the source code is
  
    cd src/archisimple93/
 
### 2.A Compile and run (Mac and Linux)
Compile the source code
  
    g++ *.cpp -std=c++11 -o archisimple   
  
Run the model
  
    archisimple   


### 2.B Compile and run (Windows)
Compile the source code

    g++ *.cpp external/tinyxml2/tinyxml2.cpp -std=c++11 -o archisimple.exe   

Run the model:

    archisimple.exe
    
    
If everything works, you should be able to run the R script in the tutorial.
