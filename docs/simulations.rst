.. _doc_cli:

1. Command-line interface (CLI)
===============================
:code:`sm_confgen` provides two command-line interfaces (CLI), including :code:`prep_ConfGen`, :code:`run_ConfGen` and :code:`analyze_ConfGen`.
:code:`prep_ConfGen` parameterizes your ligand and runs energy minimzation and equilibration at each specified temperature. :code:`run_ConfGen` 
runs the TREMD simulation in GROMACS and :code:`analyze_ConfGen` performs the analysis necessary to segment out the unique conformers which were 
sampled in the TREMD simulation.

.. _doc_prep_ConfGen:

1.1. CLI :code:`prep_ConfGen`
------------------------------
::
    usage: prep_ConfGen [-h] [-y YAML] [-o OUTPUT] [-m MAXWARN]

    This code prepares a ConfGen TREMD simulation given necessary inputs.

    options:
      -h, --help            show this help message and exit
      -y YAML, --yaml YAML  The input YAML file that contains REXEE parameters. (Default: params.yaml)
      -o OUTPUT, --output OUTPUT
                        The output file for logging how replicas interact with each other. (Default: ConfGen_log.txt)
      -m MAXWARN, --maxwarn MAXWARN
                        The maximum number of warnings in parameter specification to be ignored.

1.2. CLI :code:`run_ConfGen`
------------------------------
::
    usage: run_ConfGen [-h] [-y YAML] [-o OUTPUT] [-m MAXWARN]

    This code runs a runs a set of TREMD simulations given necessary inputs.

    options:
      -h, --help            show this help message and exit
      -y YAML, --yaml YAML  The input YAML file that contains REXEE parameters. (Default: params.yaml)
      -o OUTPUT, --output OUTPUT
                        The output file for logging how replicas interact with each other. (Default: ConfGen_log.txt)
      -m MAXWARN, --maxwarn MAXWARN
                        The maximum number of warnings in parameter specification to be ignored.

1.3. CLI :code:`analyze_ConfGen`
------------------------------
::
    usage: analyze_ConfGen [-h] [-y YAML] [-o OUTPUT] [-m MAXWARN]

    This code analyzes the output from TREMD to generate conformers given necessary inputs.

    options:
      -h, --help            show this help message and exit
      -y YAML, --yaml YAML  The input YAML file that contains REXEE parameters. (Default: params.yaml)
      -o OUTPUT, --output OUTPUT
                        The output file for logging how replicas interact with each other. (Default: ConfGen_log.txt)
      -m MAXWARN, --maxwarn MAXWARN
                        The maximum number of warnings in parameter specification to be ignored.

.. _doc_input_yaml_parameters:
2. Input YAML parameters
========================
  - :code:`gmx_executable`: (Required)
        This must be the executable for a non-MPI enabeled GROMACS for :code:`prep_ConfGen` and nn MPI enabeled GROMACS for :code:`run_ConfGen`. Either will 
        work for :code:`analyze_ConfGen`, but the non-MPI enabeled installation is recommended.
  - :code:`input_structure` : (Required)
        The filename and path for the input structure for which you want to run analysis. We recommend `sdf` file format to ensure accurate parameterization.
  - :code:`output_name` : (Required)
        The name prefix which will be applied to all output files.
  - :code:`runtime_args` : (Optional : Default=None)
        The arguments which will be passed to GROMACS during the :code:`gmx mdrun` step.
  - :code:`grompp_args` : (Optional : Default=None)
        The arguments which will be passed to GROMACS during the :code:`gmx grompp` step.        
  - :code:`temp_range` : (Optional : Default=[300, 306, 311.08, 316.23, 321.45, 326.74, 332.10, 337.54, 343.05, 348.63, 354.29, 360.02, 365.84, 371.73, 377.71, 383.76, 389.90, 396.12, 402.43, 408.83, 415.31, 421.88, 428.55, 435.30, 442.14, 449.07, 450])  
        The temperature for each temperature replica. Both the number of the value of the temperature for each replica can be changed. These values are generic and may not work well for all systems.
  - :code:`nst_sim` : (Optional : Default=None)
        The number of simulation steps to run for the production simulation. This number will be mul;tiplied by the time step to determine the simulation length. The default production simulation length is 100 ns.
  - :code:`df` : (Optional : Default=None)
        The time step to run at if a time step other than the 0.002 ps or 2 fs is desired.
  - :code:`water_model` : (Optional : Default=TIP3P)
        The water model used for the simulation. Only water models accepted by OpenFF Interchange can be input.
  - :code:`replex_rate` : (Optional : Default=1000)
        The rate for which exchanges will be attempted between adjacent temperature replicas. This rate is in units of time steps.
  - :code:`conformer_threshold` : (Optional : Default=1)            
        The number of frames for which a confomer must be sampled in order to create a new conformer group.
  - :code:`peak_threshold` : (Optional : Default=0.0005) 
        The proability minima which must be reached in order to differentiate a peak in the torsion proability distributions.
  - :code:`min_probability` : (Optional : Default=0.001)   
        The minimum probabilty that a conformer must have in order to distingush a new conformer group.