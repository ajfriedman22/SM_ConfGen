.. _doc_basic_idea:

1. Basic idea
=============
There are many ways for which you can generate ligand conformers for the purpose of performing small molecule docking or utilizing 
other protein design tools. The primary method by which these ligand conformers are traditionally generated is through the use of 
fixed torsion drives which perform incremental rotations on all bonds within the ligands. Various scoring metrics are then utilized 
to minimize the number of conformers to a number set by the user. This method produces a large number of ligand conformers which have 
little garuntee to correspond to the sampled conformers of the ligand either in solution or in complex with the protein of interest. 
This method of generating ligand conformers is significantly more computationally expensive, but has been found to be useful for 
generating ligand conformers for protein design efforts. This method is though to be particularly useful for more flexible ligands 
which contain a large number of rotatable torsions. In this situation traditional ligand conformer generation methods would either 
produce a large number of conformers for which the user would be unable or unwilling to utilize all of them or the user would need 
to rely heavily on the ranking metric utilized by the conformer generation method which is often inaccurate. This conformer generation 
method provides a relatively small number of distinct conformers which are shown to be sampled by the ligand in solution and are thus 
more likely to be similar to those conformations sampled in complex with the target protein (or other biomolecule).

This method utilizes TREMD to generate conformers. For details on the implementation of replica exchange in GROMACS please visit the 
`GROMACS documentation`. The basic principle of TREMD is that a series of MD simualtions are performed in parallel to one another at 
.. _`GROMACS documentation`: https://manual.gromacs.org/2024.0/reference-manual/algorithms/replica-exchange.html
varying temperatrures and at a set interval conformations are able to exchange between two adjacent temperature replicas (Figure 1). 
Typically the temerature of interest is the lowest temperature sampled and the higher energy temperature replicas are utilized in 
order to surrmount large free energy barriers in configurational space. The increase in kinetic energy in these higher temperature 
replicas increase the probability of crossing a large free energy barrier. The exchange between replicas allows for conformations 
visited in higher temperature replicas to be sampled in lower temperture replicas without needing to run simulations longer.

.. figure:: _static/images/TREMD.png
   :name: Fig. 1
   :width: 800
   :align: center
   :figclass: align-center

   **Figure 1.** Schematic representation of TREMD simulations in which various replica simulation are completed at increasing temperature 
   and complete conformational swaps in order to share configurational sampling with other replicate simulations to increase sampling.


.. _doc_selecting_temperature_windows:

2. Selecting Temperature Windows
================================
The selection of temperature windows is determined by the ligand of interest. There are defualt values provided in the example MDP, but 
custom temperature windows can easily be provided in the `input.yaml` file. These temperature windows should be selected so as to ensure 
the swap acceptance ratio is not too low. A general rule of thumb is to keep the swap acceptance ratio between 0.2-0.3. A short (1-2 ns) 
`sm_confgen` simulation can be performed in order to tune the temperature replica spacing. The maximum and minimum temperature are also 
system dependent. In general you should make your minimum temperture to be a physiological temperature which generally range between 300 
to 310 K. The maximum temperture should generally be around 400-450 K in order to have replicas at significantly higher temperatures, but 
not too high such that the solvent begins to move into a gas phase for which solvent parameters are not typically well suited.

.. _doc_analyzing_conformers:

2. Analyzing Conformers
=======================
The analysis of the TREMD simulations is automated by the ::code::`sm_confgen` package. The analysis has three basic steps: trajectory processing, 
conformer sorting, and conformer clustering. 
1. The trajectory processing centers the ligand in the solvent box, removes rotational and translational motion, and eliminates breaks across 
the periodic boundary condition which may occur during the simualtion.

2. The confomer sorting analyzes all flexible dihedrals within the molecule excluding those involving hydrogen atoms. From these distributions we 
can determine which flexible dihedrals are multi-modal meaning that there are multiple peaks in the probability distribution. Then we sort each 
sampled conformer by which combination of dihedrals are sampled in the lowest energy simulation. 
    * We then compute the relative free energy of each confomer using the Boltzmann equation $\Delta G = kT*ln(\frac{p_i}{p_j})$ where 
    $p_j$ is the probability of the lowest energy conformer. 
    * We also select a representative conformer from within each diedral group which is output to a PDB file. 

3. We then perform hierarchical clustering on all representative conformers from each unique conforer group determined in step 2. This step further 
reduces the conformers to those which are structurally distinct from one another. This is useful if the intended use for the output ligand conformers 
is computationally expensive and therefore requires the minimum number of ligand conformers. Otherwise we suggest using all dihedral conformers 
determined in step 2.

Below we show an example of the ligand conformers determined in step 2 for the ligand mandipropamid.


.. figure:: _static/images/example_output.png
   :name: Fig. 1
   :width: 800
   :align: center
   :figclass: align-center

   **Figure 2.** We show the 5 lowest energy conformers which were generated using ::code::`sm_confgen` and plot the reltive free energy 
   difference between the lowest energy confomer for all sampled conformers. The star represents the conformer which is sampled in the 
   Mandipropamid -- PYR1 complex conformation. 