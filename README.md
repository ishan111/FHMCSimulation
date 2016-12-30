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

1. C++11 support required

---

## Running

* See helper/window_helper.py for python functions to construct windows
* $ ./binary_name input.json 2>> err >> log

### Output

The simulation produces checkpoints which allow the simulation to be restarted automatically.  During the "production" TMMC phase, extensive properties, lnPI, as well as energy and particle number histograms collected at each value of the scalar order parameter used.  These are later automatically identified and analyzed by [FHMCAnalysis](https://mahynski.github.io/FHMCAnalysis/).

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
    "moves" : {
        "ins_del_1" : 0.1, # Relative weight to apply to insert/delete moves of species 1
        "swap_1_2" : 0.2, # Relative weight to apply to swap moves between species 1 and 2
	    "displace_1" : 0.7, # Relative weight to apply to displacement moves of species 1
	    "ins_del_2" : 0.7, # Relative weight to apply to insert/delete moves of species 2
        "displace_2" : 0.1, # Relative weight to apply to displacement moves of species 2
        "max_displacement_1" : 0.5, # Maximum displacement for species 1
        "max_displacement_2" : 0.4 # Maximum displacement for species 2
        },
    "ppot_1_1" : "square_well", # 1-1 interaction, in this case, square well
    "ppot_1_1_params" : {
        "sigma" : 1.0,
        "width" : 0.3,
        "depth" : 1.0,
        "cell_list"  : true
    },
    "ppot_1_2" : "square_well", # 1-2 interaction, in this case, square well
    "ppot_1_2_params" : {
        "sigma" : 1.1,
        "width" : 0.4,
        "depth" : 1.1,
        "cell_list"  : true
    },
    "ppot_2_2" : "square_well", # 2-2 interaction, in this case, square well
    "ppot_2_2_params" : {
        "sigma" : 1.2,
        "width" : 0.5,
        "depth" : 1.2,
        "cell_list"  : true
    },
    "barriers" : {
        "mybarrier1" : {
            "type" : "hard_wall_z", # A barrier that interacts with species 1 as a hard wall in the z-direction
            "species" : 1,
            "ub" : 6.0,
            "lb" : 1.0,
            "sigma" : 1.0
        },
        "mybarrier2" : {
            "type" : "hard_wall_z" # A barrier that interacts with species 2 as a hard wall in the z-direction
            "species" : 2
            "ub" : 6.0,
            "lb" : 1.0,
            "sigma" : 1.0
	    }
    }
}
```
