
# Evolutionary Computation methods applied to protein structure prediction and protein folding modeling with the Rosetta Software Suite

My research focuses on the protein structure prediction problem and on the computational modeling of the protein folding process. In the former, the goal is to obtain the three-dimensional structure of the protein only with the amino acid sequence information (ab initio prediction) by using evolutionary computation methods. Since the structure defines the function of a protein, this would allow a computational drug design. More information about the protein structure prediction problem can be found at Wikipedia. More information about my research can be found at my publications or, in case of more detailed information, please, feel free to contact me. 

 This project is a C++ implementation using the Rosetta Software suite tool.  
 
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

In this work Differential Evolution (DE) was combined with fragment replacement for improving the search of protein structure conformations with minimum energy. The Rosetta environment was used, employing its coarse-grained representation and its fragment-assembly technique. DE provides a global search in the multimodal energy landscape whereas the fragment replacement technique, with a Metropolis Monte-Carlo procedure, provides a useful local search that locally refines protein conformations and that accelerates the DE search.

## Differential Evolution

Differential Evolution [1] is a population-based search method. DE creates new candidate solutions by combining existing ones according to a simple formula of vector crossover and mutation, and then keeping whichever candidate solution has the best score or fitness on the optimization problem at hand.

![Differential Evolution Scheme](https://github.com/danielvarela/RosettaEvolution/blob/master/images/DE_scheme_improved.PNG)
 
## Code description 
 
The code includes the integration of classic niching methods (crowding, fitness sharing and speciation) into a hybrid version of Differential Evolution (DE) for protein structure prediction. For protein representation, the Rosetta coarse-grained representation model was used. Rosetta is one of the most successful software environments for protein design [5]. The hybrid DE version incorporates the Rosetta fragment replacement technique as a local search operator.  

Given the inaccuracies of the Rosetta energy model, the inclusion of niching allows the simultaneous search in different areas of the energy landscape that correspond to different minima, with the aim to obtain a diversified set of optimized (native-like) folds.  

## Source code tree: 

In the “src” directory, the code used can be found, with the following structure:  

* ***main.cpp*** 

* ***Algorithm***: Differential Evolution and niching algorithms were integrated.  It includes the classes for *CrowdingDE, FitnessSharingDE* and *Species-BasedDE (SDE)*. Also, DE_types.hpp defines the types used in the differential evolution algorithm versions, for example, *Population* or *Individual*. 

* ***Controller***: Contains the classes that start the program, read the input files and runs the program with the options indicated in the option.ini file. 

* ***Movers***: Classes with the different parts of the main differential evolution algorithm, for example, initialize population, print the information, insert fragments at a vector of individuals, etc...  

* ***MpiFiles***: Source files for the MPI implementations (master-slave and gather-scatter schemes). 

* ***Extra***: additional functions that are commonly used in different files. 

* ***test_main.cpp***: main file for running the tests. 

* ***Test***: Test files for the main operations. 

 
## Short notes about the different algorithms: 

Differential evolution and niching methods integration 

The implementations of Differential Evolution and the niching algorithms (crowding, fitness sharing and speciation) can be located at the file: 
```
    src/Algorithm/DifferentialEvolutionMover.hpp
```

### Differential Evolution 

The classic implementation of Differential Evolution [3] is used, following the scheme of the previous figure. The source code can be found in class **MoverDE**, which is used as a base class for the rest of the implementations. The main function of the file is called apply. Those names follow the code guidelines of the Rosetta code.  

 
```cpp
MoverDE::apply() 

std::vector<Individual> popul; // initial population 
while (gen_count < Gmax) { 
     reset_stat();
     trial_popul.resize(0);
     
     for (int i = 0; i < NP; ++i) { 
            select_parents(i, parents); 
            Individual ind = sample_new_individual(i, parents, popul); 
            trial_popul.push_back(ind); 
     } 
    
     trial_popul = mpi_calculator->run(trial_popul); // evaluate the fitness of the population with mpi parallelized version 
     new_best_found = select_population(trial_popul); 
     print_generation_information(gen_count, new_best_found); 
     ++gen_count; 
} 
```
 

### Differential Evolution with fragments replacement 

As it is defined in our preliminary work [2], we combined the global search of Differential Evolution with Rosetta’s fragment replacement.  Fragment replacement techniquecan be considered as a local search since it locally refines the dihedral angles of a conformation. More information can be found in [3].  

The corresponding class is named as **HybridMoverDE**. The fragment replacement operation is controlled with a boolean variable "frags_at_popul", which, in this class, is always set to True.  The code of the HybridMoverDE is the same as the **MoverDE**, but using the appropriate function of the MPI parallelized version for evaluating the fitness of the individuals of the population.  

```cpp
HybridMoverDE::apply()
      trial_popul = mpi_calculator->run(trial_popul, frags_at_popul); 
```

### Crowding 

The integration of the crowding niching method and DE, defined by Thomsen [6], was followed, extending DE with the classic crowding scheme.  The implementation is in class **CrowdingMoverDE**, which replaces the base function "select_population". Now, instead of comparing each trial individual with its corresponding individual in the population, the trial individual is compared against its nearest neighbour (most similar) among a subset of the population. For this reason, a distance function is needed in order to differentiate protein folds and find the most similar individual.  

 
```cpp
CrowdingMoverDE::select_population 

for (int i = 0; i < trial_popul.size(); i++) { 
      int nearest_ind = calculate_distances_popul->find_nearest(trial_popul[i], popul); 
      if ( trial_popul[i].score < popul[nearest_ind].score ) { 
            popul[nearest_ind] = trial_popul[i]; 
            trial_sucess_n++; 
      } 
} 
```
 

The distance functions used are located in file "src/Movers/CalculateDistancePopulation.cpp". Three functions are used: 1. Root Mean Square Deviation (RMSD), 2. The structural diversity measure (SDM) defined by Garza-Fabre et al [4], which describes the relative position of each pair of secondary structure elements with respect of each other and 3. The Template Modeling score (TM-score) defined by Zhang et al [8].  


### Fitness Sharing 

The fitness sharing (FS) niching method was also integrated into DE to define a SharingDE version, following the implementation defined by Thomsen [6]. The basic idea of FS is to punish individuals that occupy the same area of the search space, rescaling the fitness of each encoded solution considering the number of individuals in its neighbourhood. The implementation in class **SharedMoverDE** modifies the "select_population" function.   

The implementation uses a **SharedFitnessIndividual** class defined in "src/Algorithm/DifferentialEvolutionMover.hpp". 


### Speciation 

The species-based DE (SDE) defined by Li [7] was used as base in our approach. In SDE, each of the "species" is built around a dominating species’ seed. All individuals that fall within the radius from the species seed are identified as the same species. Since DE mutation is carried out within each species, the technique has the ability to maintain high diversity and stable niches over generations.  

The class **SeedsMoverDE** defines the code used for this SDE version. As all the other classes, the apply function carries out the main operations of the algorithm. In this case, it makes use of a new function "create_seeds", which is responsible for creating the different populations of species. 


## References 

[1] C. Rohl, C. Strauss, K. Misura, D. Baker, Protein structure prediction using Rosetta, Methods in Enzymology 383 (2004), 66-93 

[2] D. Varela, J. Santos, Combination of Differential Evolution and Fragment-based Replacements for Protein Structure Prediction
Proceedings of the Companion Publication of the 2015 on Genetic and Evolutionary Computation Conference - GECCO Companion '15
2015

[3] K. Price, R. Storn, J. Lampinen, Differential evolution. A practical approach to global optimization, Springer – Natural Computing Series, 2005. 

[4] M. Garza-Fabre, S. Kandathil, J. Handl, J. Knowles, S. Lovell, Generating, maintaining, and exploiting diversity in a memetic algorithm for protein structure prediction, Evolutionary Computation 24(4) (2016). 577-607  

[5] Rosetta system, www.rosettacommons.org 

[6] R. Thomsen, Multimodal optimization using crowding-based differential evolution, in: Proceedings IEEE Congress on Evolutionary Computation, 2004, pp. 1382-1389 

[7] X. Li, Efficient differential evolution using speciation for multimodal function optimization, in: Proceedings GECCO 2005 - Conference on Genetic and Evolutionary Computation, 2005, pp. 873 – 880 

[8] Y. Zhang and J. Skolnick. Scoring function for automated assessment of protein structure template quality. Proteins: Structure, Function, and Bioinformatics, 57(4):702–710, 2004. 


## List of Publications

In order to know more about the methods used at this project, please, find the detailed information in the following publications: 

* Varela, D. and Santos, J. (2019), "Crowding differential evolution for protein structure prediction", *Proceedings International Work-Conference on the Interplay between Natural and Artificial Computation - IWINAC 2019*, *Lecture Notes in Computer Science* 11487:193-203. DOI: 10.1007/978-3-030-19651-6_19

* Varela, D. and Santos, J. (2018), "Automatically obtaining a cellular automaton scheme for modeling protein folding using the FCC model", *Natural Computing*, DOI: 10.1007/s11047-018-9705-y.

* Varela, D. and Santos, J. (2017), "A hybrid evolutionary algorithm for protein structure prediction using the Face-Centered Cubic lattice model", *Proceedings International Conference on Neural Information Processing ICONIP*, *Lecture Notes in Computer Science*, 10634:628-638. DOI: 10.1007/978-3-319-70087-8_65

* Varela, D. and Santos, J. (2017), "A protein folding model using the Face-Centered Cubic lattice model", *GECCO 2017 ACM Proceedings Companion, Workshop Evolutionary Computation in Computational Biology, Genetic and Evolutionary Computation Conference*, 1674-1678. 
DOI: 10.1145/3067695.3082543

* Varela, D. and Santos, J. (2017), "Protein folding modeling with neural cellular automata using the Face-Centered Cubic model", *Proceedings International Work-Conference on the Interplay between Natural and Artificial Computation*, *Lecture Notes in Computer Science*, 10337:125-134. DOI: 10.1007/978-3-319-59740-9_13

* Varela, D. and Santos, J. (2016), "Protein folding modeling with neural cellular automata using Rosetta", *GECCO 2016 ACM Proceedings Companion, Workshop Evolutionary Computation in Computational Structural Biology, Genetic and Evolutionary Computation Conference*, 1307-1312. DOI: 10.1145/2908961.2931720

* Varela, D. and Santos, J. (2015), "Combination of differential evolution and fragment-based replacements for protein structure prediction", *GECCO 2015 ACM Proceedings Companion, Workshop Evolutionary Computation in Computational Structural Biology, Genetic and Evolutionary Computation Conference*, 911-914. DOI: 10.1145/2739482.2768437
