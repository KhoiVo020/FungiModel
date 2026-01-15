# Fungal Structure Code

## Before attempting to run the code:
- Go to the line "cwd_path ='/Users/libra/FBI_project/fungalGrowthModel_singleNutrient_py/0521'" in the driver_fungalGrowth_singleNutrient.py and change the path to where you  place the folder containing all files in this repository.
- Go to "def get_filepath(params):" in help_functions.py to change the name of folder string and file string.

## To Run the code:
- **Option 1:** Just run the driver_fungalGrowth_singleNutrient.py file 
    - The number of iterations of runs is determined by the variable ‘num_runs’ (set around line 269 - later will be an input if this file changes to a function)
    - If num_runs==1, then it will run just one iteration of a fungal mycelia growth
    - If num_runs>1, then it will run multiple iterations using the same set of parameters in parallel using the python Parallel function from the joblib package. Currently we do not need Parallel function yet.
- **Option 2:** run the function (from within the py file listed in Option 1) using the command ‘driver_singleNutrient(1)’ to run one iteration of fungal mycelia growth.

## About the various files:
- **driver_fungalGrowth_singleNutrient.py:** contains the main driver file that executes the steps of diffusion in the external domain, elongation, branching, fusion, translocation, and uptake.
- **parameters.ini:** where all the parameters values are stored. If you want to change a parameter, it is most likely listed in this file.
- **helper_functions.py:** contains functions that do the following
    - convert parameters from the ini file to a usable form 
    - generates filenames used for saving data specific to simulation with given parameters
    - plotting funtions (the fungal structure, the external domain, various output stats)
- **setup_functions.py:** contains functions that do the following
    - set up the dictionary for storing info related to the mycelia structure
    - sets up the initial fungal structure and it’s properties
    - sets up the external grid shape and amount of nutrient in each cell
- **growth_functions.py:** contains functions that pertain to
    - elongation of hyphae at the tip
    - branching 
    - anastomosis (fusion) of hyphae
- **nutrient_functions2.py:** contains functions that pertain to
    - translocation update
    - uptake of nutrients
    - also contain diffusion method: currently using finite_volume_method
      
  ## To Run the code on HPCC:
    - change the directory in **nutrient_functions2.py** and **driver_fungalGrowth_singleNutrient.py**, then run the command dos2unix your_script_name.sh to ensure environmental match. Then run sbatch submission_script_2.sh to start the simulation.
