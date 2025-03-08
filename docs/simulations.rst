.. _doc_cli:

1. Command-line interface (CLI)
===============================
:code:`sm_confgen` provides two command-line interfaces (CLI), including :code:`prep_ConfGen`, :code:`run_ConfGen` and :code:`analyze_ConfGen`.
:code:`prep_ConfGen` parameterizes your ligand and runs energy minimzation and equilibration at each specified temperature. :code:`run_ConfGen` 
runs the TREMD simulation in GROMACS and :code:`analyze_ConfGen` performs the analysis necessary to segment out the unique conformers which were 
sampled in the TREMD simulation.

.. _doc_explore_REXEE:

1.1. CLI :code:`explore_REXEE`
------------------------------

.. _doc_implemted_workflow:
2. Implemented workflow
=======================

.. _doc_input_yaml_parameters:
3. Input YAML parameters
========================