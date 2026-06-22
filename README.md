# README for catchyCFDEM
catchyFOAM: Laurien Vandewalle, August 2020
catchyCFDEM-v1.0 : Florian Wéry (), January 2024 

This is an adapted version of the original catchyCFDEM (catchyCFDEM (CATCHY = CATalytic CHemistrY) is an extension of the catchyFOAM framework) directory, authored by Hamed Hoorijani.  Up to date installation instructions to be written.
This update includes the dynamic particle electric conduction and joule heating within electrified particulate systems by integerating search algorithms within the previously released version of catchyCFDEM by Florian Wery. 

## Installation prerequisites
catchyFOAM (CATCHY = CATalytic CHemistrY) provides a set of OpenFOAM libraries, solvers and utilities for problems involving heterogeneous catalytic chemistry.<br/><br/>
The most important developments are made in the thermophysicalModels library, which is at the basis of many other libraries (such as combustionModels, ThermophysicalTransportModels, ...). The necessary libraries depending on thermophysicalModels are therefore copied from the local OpenFOAM installation, linked to the catchyFOAM thermo and recompiled. These *new* libraries are then used in the catchyFOAM solvers and utilities.

## Installation prerequisites
### OpenFOAM
catchyFOAM is compatible with OpenFOAM-8 (OpenFOAM Foundation). For download instructions, see [https://openfoam.org/version/8]. 

### (optional) Cantera
catchyFOAM comes with a utility to convert Cantera mechanism files into OpenFOAM format, <span style="font-family:Courier;">canteraToFoam</span> (as an alternative to the <span style="font-family:Courier;">chemkinToFoam</span> utility available by default). To install <span style="font-family:Courier;">canteraToFoam</span>, Cantera needs to be installed.

Get the Cantera source code from Github and switch to the correct branch via
```
    git clone https://github.com/lavdwall/cantera.git --branch gas_transport_access
```

Detailed installation instructions can be found [here](https://cantera.org/install/compiling-install.html). 
For Ubuntu, the following commands should do the trick
```
    sudo apt-get install g++ gfortran scons libboost-dev python3 python3-dev python3-setuptools python3-numpy python3-ruamel.yaml python3-pip
    pip3 install Cython
    
    cd cantera
    scons build prefix='$HOME/.local' python_cmd=/usr/bin/python3
    scons install
    cd ..
```

## Install catchyFOAM
Get the catchyFOAM source code from Github 
```
    git clone https://github.com/lavdwall/catchyFOAM.git
```
By default, it is assumed that the catchyFOAM repository is located in *$HOME/OpenFOAM*. If this is not the case, adjust *catchyFOAM/etc/bashrc* accordingly. Also, if Cantera is installed in a non-default location, it can be set in *catchyFOAM/etc/bashrc*.

Navigate to the catchyFOAM repository and run
```
    ./Allwmake
```
This will install all libraries, solvers and utilities.

To use catchyFOAM, the catchyFOAM environment should be set by sourcing the *catchyFOAM/etc/bashrc* file. This can be done by adding the following line at the end of the user's *.bashrc* file (to be adjusted accordingly if catchyFOAM is not located in *$HOME/OpenFOAM*)
```
    source $HOME/OpenFOAM/catchyFOAM/etc/bashrc
```

## Uninstall catchyFOAM
To clean the installation, run
```
    ./Allwclean
```
## Install catchyCFDEM
Get the catchyCFDEM source code from Github

    git clone https://github.com/FlorianWery/catchyCFDEM.git
By default, it is assumed that the catchyCFDEM repository is located in $HOME/OpenFOAM. If this is not the case, adjust catchyCFDEM/etc/bashrc accordingly. Also, if Cantera is installed in a non-default location, it can be set in catchyCFDEM/etc/bashrc.

To use catchyCFDEM, the catchyCFDEM environment should be set by sourcing the catchyCFDEM/etc/bashrc and catchyCFDEM/CFDEM/CFDEMcoupling/etc/bashrc files. This can be done by adding the following line at the end of the user's .bashrc file (to be adjusted accordingly if catchyCFDEM is not located in $HOME/OpenFOAM)

    source $HOME/OpenFOAM/catchyCFDEM/etc/bashrc
    source $HOME/OpenFOAM/catchyCFDEM/CFDEM/CFDEMcoupling/etc/bashrc
Navigate to the catchyCFDEM repository and run

    ./Allwmake
This will install all libraries, solvers and utilities of catchyFOAM.
