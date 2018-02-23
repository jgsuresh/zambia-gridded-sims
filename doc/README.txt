This code base runs DTK simulations of malaria transmission using nodes which correspond to 1km x 1km grid cells.

Code structure:
- experiment_setup.py is the workhorse, containing:
    - Class GriddedConfigBuilder which generates a config-builder for a spatial sim
    - Class GriddedInputFilesCreator which takes this config-builder, as well as a gridded population CSV
      file and generates the input files for the simulation

Instructions:
- For a new setting, write a new class which inherits from GriddedConfigBuilder and contains site-specific functions
  (see zambia_experiments.py or mozambique_experiments.py for examples)
- Write a run script which calls your new class, uses GriddedInputFilesCreator to create the corresponding input files,
  and then submits the job for execution.
