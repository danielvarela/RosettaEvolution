
# Evolutionary Computation methods applied to protein structure prediction and protein folding modeling with the Rosetta Software Suite

My research focuses on the protein structure prediction problem and on computational modeling of the protein folding process. In the formerm, the goal is to obtain
the three-dimensional structure of the protein only with the amino acid sequence information (ab initio prediction) by using evolutionary computation methods. Since the structure defines the function of a protein, this would allow a computational drug design. More information about the protein structure prediction problem can be found at Wikipedia. More information about my research can be found at my publications or, in case of more detailed information, please, feel free to contact me.

This project is a C++ implementation using the Rosetta Software Suite tool. A simplified version is implemented in python implementation the same methods but with a simple 2D lattice and known python libraries. [2D_protein_AI](https://github.com/danielvarela/2D_protein_AI).

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


## Basic Usage

First, you should include the protein input files in the code. Modify the file src/Controller/ProteinInfoArchive.cpp following one of the existing examples. For example:

```cpp
Protinfo c8c;
c8c.name = "1c8cA";
c8c.pdb_file = "./input_files/info_1c8cA/vf_1c8c.pdb";
c8c.frag_3 = "./input_files/info_1c8cA/boinc_vf_aa1c8cA03_05.200_v1_3";
c8c.frag_9 = "./input_files/info_1c8cA/boinc_vf_aa1c8cA09_05.200_v1_3";
c8c.ss_file = "./input_files/info_1c8cA/vf_1c8cA.psipred_ss2";
prot_selection["1c8cA"] = c8c;
```

Compile the code and modify the options.ini in order to indicate the input protein that you selected. 

Use the command line:

```
$RosettaEvolution/bin/mpi_app options.ini
```

or


```
mpirun -np X $RosettaEvolution/bin/mpi_app options.ini
```

where *X* is the number of processes to run the application.


# Protein Structure Prediction

n this work Differential Evolution (DE) was combined with fragment replacement for improving the search of protein structure conformations with minimum energy. The Rosetta environment was used, employing some of its phases for the ab initio prediction in the initialization of the genetic population, as well as its fragment-assembly technique. DE provides a global search in the multimodal energy landscape whereas fragment replacement based on the Monte-Carlo procedure provides a useful local search that locally refines protein conformations and that accelerates the DE search.

## Differential Evolution

Differential Evolution [Price05] is a population-based search method. DE creates new candidate solutions by combining existing ones according to a simple formula of vector crossover and mutation, and then keeping whichever candidate solution has the best score or fitness on the optimization problem at hand.

![Differential Evolution Scheme](https://github.com/danielvarela/RosettaEvolution/blob/master/images/DE_scheme_improved.PNG)
 
## List of Publications

In order to know more about the methods used at this project, please, find the detailed information in the following publications:

* Crowding Differential Evolution for Protein Structure Prediction
From Bioinspired Systems and Biomedical Applications to Machine Learning
2019 | book-chapter
DOI: 10.1007/978-3-030-19651-6_19


* Automatically obtaining a cellular automaton scheme for modeling protein folding using the FCC model
Natural Computing
2018-08 | journal-article
DOI: 10.1007/s11047-018-9705-y


* A Hybrid Evolutionary Algorithm for Protein Structure Prediction Using the Face-Centered Cubic Lattice Model
Neural Information Processing
2017 
DOI: 10.1007/978-3-319-70087-8_65

* A protein folding model using the face-centered cubic lattice model
Proceedings of the Genetic and Evolutionary Computation Conference Companion on - GECCO '17
2017 
DOI: 10.1145/3067695.3082543

* Protein Folding Modeling with Neural Cellular Automata Using the Face-Centered Cubic Model
Natural and Artificial Computation for Biomedicine and Neuroscience
2017 
DOI: 10.1007/978-3-319-59740-9_13

* Protein Folding Modeling with Neural Cellular Automata Using Rosetta
Proceedings of the 2016 on Genetic and Evolutionary Computation Conference Companion - GECCO '16 Companion
2016 
DOI: 10.1145/2908961.2931720

* Combination of Differential Evolution and Fragment-based Replacements for Protein Structure Prediction
Proceedings of the Companion Publication of the 2015 on Genetic and Evolutionary Computation Conference - GECCO Companion '15
2015
DOI: 10.1145/2739482.2768437
