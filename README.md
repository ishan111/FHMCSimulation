# Flat Histogram Monte Carlo Simulation

Nathan A. Mahynski

---

Status

[![Build Status](https://travis-ci.org/mahynski/FHMCSimulation.svg?branch=master)](https://travis-ci.org/mahynski/FHMCSimulation) [![Open Source Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badge/) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/f5b0edf4e77e4902b871d7f1faeabc6f)](https://www.codacy.com/app/nathan-mahynski/FHMCSimulation?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=mahynski/FHMCSimulation&amp;utm_campaign=Badge_Grade)

---

## Installation:

```
$ git clone https://github.com/mahynski/FHMCSimulation.git
$ cd FHMCSimulation/build
$ cmake CMakeLists.txt; make
```

---

## Unittests

```
$ cd FHMCSimulation/tests/src
$ cmake CMakeLists.txt; make
$ cd ../bin; ./runUnitTests
```

---

## Dependencies

1. GCC
2. C++11 support required

---

## Running

* See helper/window_helper.py for python functions to construct windows
* $ ./binary_name input.json 2>> err >> log

### Input File

Example input file(s) may be generated with helper/window_helper.py, and can be modified by the user as needed.  A sample input.json file is given below:

```python
{
"num_species" : 2, # 2 species
"beta" : 0.89, # 1/T* = 0.89
"box" : [9.0, 
  9.0, 
  9.0], # Box = [9, 9, 9] with a corner at (0,0,0)
"mu" : [1.234, 
  2.345], # (Excess) Chemical potential of each species
"seed" : -10, # RNG seed
"max_N" : [200, 
  200], # Maximum number of each species type to allow
"min_N" : [0, 
  0], # Minimum number of each species type to allow
"window" : [10, 
  30], # Window to sample, [N_tot,min , N_tot,max]
"restart_file" : "", # File to read configuration from
"num_crossover_visits" : 50, # Number of times for the TMMC to sample each point in C matrix before taking over from WL
"tmmc_sweep_size" : 1e3, # Number of times each point in C matrix must be visited per sweep
"total_tmmc_sweeps" : 300, # Total number of sweeps to perform during TMMC stage
"wala_sweep_size" : 10e4, # Number of Monte Carlo moves per sweep in WL stage (after each sweep) flatness is checked for
"wala_g": 0.5, # WL reduction factor
"wala_s" : 0.8, # WL flatness criterion
"lnF_start" : 1.0, # Where to start WL
"lnF_end" : 1.0e-5, # Where to end WL
"prob_pr_ins_del_1" : 0.1, # Relative weight to apply to insert/delete moves during production (TMMC) for species 1
"prob_pr_swap_1_2" : 0.2, # Relative weight to apply to swap moves during production (TMMC) for species 1 and 2
"prob_pr_displace_1" : 0.7, # Relative weight to apply to displacement moves during production (TMMC) for species 1
"prob_pr_ins_del_2" : 0.7, # Relative weight to apply to insert/delete moves during production (TMMC) for species 2
"prob_pr_displace_2" : 0.1, # Relative weight to apply to displacement moves during production (TMMC) for species 2
"prob_eq_ins_del_1" : 0.25, # Relative weight to apply to insert/delete moves during equilibration (WL/Crossover) for species 1
"prob_eq_swap_1_2" : 0.15, # Relative weight to apply to swap moves during equilibration (WL/Crossover) for species 1 and 2
"prob_eq_displace_1" : 0.6, # Relative weight to apply to displacement moves during equilibration (WL/Crossover) for species 1
"prob_eq_ins_del_2" : 0.4, # Relative weight to apply to insert/delete moves during equilibration (WL/Crossover) for species 2
"prob_eq_displace_2" : 0.6, # Relative weight to apply to displacement moves during equilibration (WL/Crossover) for species 2
"max_pr_displacement_1" : 0.5, # Maximum displacement for species 1 during production (TMMC)
"max_eq_displacement_1" : 0.6, # Maximum displacement for species 1 during equilibration (WL/Crossover)
"max_pr_displacement_2" : 0.4, # Maximum displacement for species 2 during production (TMMC)
"max_eq_displacement_2" : 0.7, # Maximum displacement for species 2 during equilibration (WL/Crossover)
"ppot_1_1" : "square_well", # 1-1 interaction, in this case, square well
"ppot_1_1_params" : [1.0, 
  0.3, 
  1.0], # 1-1 interaction parameters
"ppot_1_1_use_cell_list" : true, # use cell list for this potential to accelerate
"ppot_1_2" : "square_well",
"ppot_1_2_params" : [1.1, 
  0.4, 
  1.1],
"ppot_1_2_use_cell_list" : true,
"ppot_2_2" : "square_well",
"ppot_2_2_params" : [1.2, 
  0.5, 
  1.2],
"ppot_2_2_use_cell_list" : true
}
```
