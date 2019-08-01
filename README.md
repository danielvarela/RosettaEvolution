
# AI methods applied to protein folding with Rosetta Software Suite

My PhD research project focuses on the computational modeling of the protein folding process. The goal is to obtain
the three-dimensional structure of the protein only with the amino acid sequence information (ab initio) by using evolutionary computing, artificial life techniques and complex systems theory methods. Since the structure defines the function of a protein, this would allow a computational drug design. More information about the protein structure prediction problem can be found at Wikipedia. More information about my research can be found at my publications or, in case of more detailed information, please, feel free to contact me.

This project is a C++ implementation using the Rosetta Software Suite tool. A simplificated version is implemented at python implementation the same methods but with a simple 2D lattice and known python libraries. ![2D_protein_AI](https://github.com/danielvarela/2D_protein_AI).

### Dependencies

* icc 2018.5.274 ( Intel Compiler)
* openmpi
* boost 1.68.0 (with mpi module compiled)
* Rosetta software suite 2018
* python 3.7.0
* cmake 3.10.3

### Install

A CMakeLists.txt file is provided in order to compile the project.

First, set the $ROSETTA_PATH environment variable to your Rosetta Software Suite compiled library.

Then

```
cd $RosettaEvolution/build
cmake ..
build mpi_app
```

to run the code, modify the options.ini file with your input files and exec the app with:

```
$RosettaEvolution/bin/mpi_app options.ini
```

or


```
mpirun -np $RosettaEvolution/bin/mpi_app options.ini
```


# Protein Structure Prediction

ToDo

# List of Publications


ToDo

