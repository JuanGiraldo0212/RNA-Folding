# Project description

RNAs are nucleotide polymers fundamental to a diversity of elementary biological functions, including the (de)coding and regulation of genes, protein construction, cellular signaling, and catalysis [1]. Integral to these roles is the capacity and propensity of RNAs to self-interact through hydrogen-bond base-pairing between nucleotides and fold into specific, stable structures [2]. This folded structure of an RNA along with its primary sequence chemistry combine to dictate its interactions with other biomolecules [3]. Understanding and predicting RNA folding is thus a pressing interest of the biological sciences, basic and applied [4, 5].

RNA folding prediction from primary sequence information alone remains challenging classically, viewed from both the standpoints of chemical dynamics and combinatorial optimization of free energy. Quantum approaches, with their demonstrable advantage in both of these realms, therefore lend themselves well to the RNA folding problem. To our knowledge, only one attempt has been made to map RNA folding to quantum computing, via quantum annealing [6]. With this project, we seek to: (1) modify the Hamiltonian presented therein to better reflect the underlying chemistry of base-pairing and pseudoknots, and (2) optimize the free parameters of the Hamiltonian against a suite of RNAs with known structures. Expanding on (2), we aim to implement and test our Hamiltonian with the quantum annealing hardware of D-Wave, and demonstrate a parallel approach with gate-based hardware via QAOA, using the very same Hamiltonian.

# Project organization

The annealing implemention of this project is located in the file: qrna_annealing.ipynb, where we use D-wave's SDK, Ocean, to build the annealing solution for the RNA folding problem.

The QAOA implementation of this project is located in the file: qrna_circuit.ipynb where we use the QAOA implementation of Qiskit to find the optimal solution for the proposed Hamiltonian. 

The training implemention of this project is located in the file: training.ipynb, where we use D-wave's simulated annealing sampler and Qiskit's SPSA optimizer to optimize the parameters of our Hamiltonian.

The `models` folder contains standalone scripts for each of model 1, model 2, and our model (model 3).


The `data` folder contains a dataset with multiple RNA from [bpRNA](http://bprna.cgrb.oregonstate.edu/index.html) sequences of different size (from 0 to 500 bases) and type (known to contain pseudoknots, pseudoknot-free). These are used as inputs for each implementation of our methodology. 

The `results` folder contains the raw results of out annealing implementation in .csv form, including: structure bpRNA ID, known stems, known energy, predicted stems, predicted energy, and comparison metrics. The `cts` folder contains all actual and predicted connectivity tables. The `figures` folder contains [Circle Compare](https://rna.urmc.rochester.edu/RNAstructureWeb/Servers/CircleCompare/CircleCompare.html) diagrams and other plots.

The `.old` folder contains older scripts (H-model, N-model, and plotting scripts) and older results/data.

The `misc` folder contains the circuit model implementation, as well as a script which plots the number of potential stems against sequence length to compare with the theoretical maximum. 

The linked [presentation](https://docs.google.com/presentation/d/1r9f51t5CkaSBdAiReTvXFvJ-J_PL5B8CKlnhDdwZ0fM/edit?usp=sharing) contains a walk-through of our methodology and results.

# References

[1] J. Li and C. Liu, “Coding or noncoding, the converging concepts of RNAs,” Frontiers in
Genetics, vol. 10, May 2019.    
[2] G.L.Conn and D.E.Draper,“RNAstructure,”Current Opinion in Structural Biology,vol.8,
no. 3, pp. 278–285, Jun. 1998.    
[3] S. R. Holbrook, “RNA structure: the long and the short of it,” Current Opinion in Structural
Biology, vol. 15, no. 3, pp. 302–308, Jun. 2005.     
[4] M. D. Disney, “Targeting RNA with small molecules to capture opportunities at the intersec-
tion of chemistry, biology, and medicine,” Journal of the American Chemical Society, vol.
141, no. 17, pp. 6776–6790, Mar. 2019.   
[5] N. G. Walter and L. E. Maquat, “Introduction—RNA: From single molecules to medicine,”
Chemical Reviews, vol. 118, no. 8, pp. 4117–4119, Apr. 2018.    
[6] D. M. Fox, C. M. MacDermaid, A. M. Schreij, M. Zwierzyna, and R. C.
Walker, “RNA folding using quantum computers,” May 2021.

# Bibliography

 - Fox et al., 2021 (RNA folding with quantum computers)
 - Lewis et al., 2021 (A new modelling paradigm for the RNA folding problem)
 - Bhatia and Zheng, 2020 (Evolving two-way quantum finite automata with classical states for RNA secondary structures)
 - Fox et al., 2021 (mRNA codon optimization with quantum computers)
 - Emani et al., 2019 (Quantum computing at the frontiers of biological sciences)
 - Fingerhuth et al., 2018 (A quantum alternating operator ansatz with hard and soft constraints for lattice protein folding)
 - Lu and Li, 2019 (Quantum approach to fast protein-folding time)
 - Robert et al., 2020 (Resource-efficient quantum algorithm for protein folding)
 - Casares et al., 2022 (Quantum walks and deep learning to solve protein folding)
 - Thirumalai and Hyeon, 2005 (RNA and protein folding: Common themes and variations)
 - Perdomo-Ortiz et al., 2012 (Finding low-energy conformations of lattice protein models by quantum annealing)
 - Li and Liu, 2019 (Codong or noncoding: The converging concepts of RNAs)
 - Conn and Draper, 1998 (RNA structure)
 - Hollbrook, 2005 (RNA structure: The long and short of it)
 - Disney, 2019 (Targeting RNA with small molecules to capture opportunities at the intersection of chemistry, biology, and medicine)
 - Walter and Maquat, 2018 (Introduction - RNA: From single molecules to medicine)
 - Waldsich, 2014 (RNA folding: Methods and protocols)
 - Turner and Mathews, 2009 (The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure)
 - Mathews, 2019 (How to benchmark RNA secondary structure prediction accuracy)
 - Bhatia and Zheng, 2020 (Evolving two-way quantum finite automata with classical states for RNA secondary structure)
 - Amin et al., 2011 (Quantum Boltzmann machine)
 - Matsumoto et al., 2021 (Calculation of Gibbs partition function with imaginary time evolution on near-term quantum computers)
 - McCaskill,1990 (The equilibrium partition function and base pair binding probabilities for RNA secondary structure)
 - Temme et al., 2011 (Quantum metropolis sampling)
 - Schuld et al., 2019 (Evaluating analytic gradients on quantum hardware)
 - Jordan, 2005 (Fast quantum algorithm for numerical gradient estimation)
 - Bittel, 2021 (Training variational quantum algorithms is NP-hard, even for logarithmically many qubits and free fermionic systems)
 - Wang et al., 2021 (Noise-induced barren plateaus in variational quantum algorithms)
 - Campbell et al., 2019 (Applying quantum algorithms to constraint satisfaction problems)
 - King et al., 2017 (Quantum annealing amid local ruggedness and global frustration)
 - Morita and Nishimori, 2008 (Mathematical foundation of quantum annealing)
 - Biswas et al., 2017 (A NASA perspective on quantum computing: Opportunities and challenges)
 - Susa et al., 2018 (Exponential speedup of quantum annealing by inhomogeneous driving of the transverse field)
 - Hauke et al., 2020 (Perspectives of quantum annealing)
 - Barkoutsos et al., 2020 (Improving Variational Quantum Eigensolver with CVaR)
 - Jorg et al., 2010 (Energy gaps in quantum first-order mean-field–like transitions: The problems that quantum annealing cannot solve)
 - Kelley and Didulo, 2018 (Computational biology - a hypertextbook)
 - Hajdin et al., 2013 (Accurate shape-directed RNA secondary structure modelling including pseudoknots)
 - Andronescu et al., 2010 (Improved free energy parameters for RNA pseudoknotted secondary structure prediction)
 - Bapst and Semerjian, 2012 (On quantum mean-field models and their quantum annealing)