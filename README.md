# Grand-canonical MD simulation with H4D method

H4D, Hybrid 4D NEMD/MC method in conjunction with LAMMPS


- Author        : Jeongmin Kim, a former postdoc at PHENIX with Benjamin Rotenberg
- Contact       : jeongmin0658 AT gmail DOT com


The sources and example folders are included:

Each folder contains its own ```readme``` file too.

```H4D```         : LAMMPS c++ and head files

```ver0.0```      : H4D Python scripts

```examples```    : Example systems, including 

- "LJ_salt100"  : Example of LJ electrolytes   
  - "gcmd"      : Grand-canonical MD simulation
  - "sampling"  : Simulation to compute chemical potential (no accept of trial moves)
    - "insertion" : Insertion trial moves
    - "deletion"  : Deletion trial moves
- "JC_SPCE_1m"  : Example of aqueous NaCl electrolytes 

---
# References

[J. Kim, L. Belloni, B. Rotenberg, Arxiv](https://arxiv.org)

[L. Belloni, J. Chem. Phys. 151, 021101 (2019)](https://pubs.aip.org/aip/jcp/article/151/2/021101/197798/Non-equilibrium-hybrid-insertion-extraction) - original H4D method

---
# LAMMPS installation with H4D
Our H4D method works well with LAMMPS of a version of 27 Oct 2021.
You can download [lammps-27Oct2021.tar.gz](https://download.lammps.org/tars/index.html).

```H4D``` folder includes all the cpp and head files for lammps installation

- *_neq.cpp
- *_neq.h
 
How to use:
1. Copy this folder to the ```scr``` folder for LAMMPS installation
2. Type ```make yes-H4D```, before you install LAMMPS
3. Then, install LAMMPS as a shared library and Python module (```make install-python```) See [Manual 2.2.1](https://docs.lammps.org/Python_install.html#installing-the-lammps-python-module-and-shared-library)
  - Installation option should be ```mode=shared```. See [Manual 3.4.3](https://docs.lammps.org/Build_basics.html#build-the-lammps-executable-and-library)
  - For example, ```make mode=shared machine```

---
# H4D run
Before running LAMMPS with H4D, you must define where the "ver0.0" folder is located. You may use ```setup.sh``` in the folder where you execute the simulation.

```run_gcmc_sample.py``` contains the parameters that control the H4D.

Simply run:

```python3 run_gcmc_sample.py```
or in case you use multiple CPUs, 

```mpirun -np number_of_cpus python3 run_gcmc_sample.py```
You can see the ```example_job_submission.sh```.

---
# LAMMPS outputs
All the outputs can be moved to the "output" folder. Use ```copy_output.sh```.
- ```logfile_hneqMDMC.*```: LAMMPS output of each "n_MC" (0 - n_MC-1).
- ```restart_after_MC_*.data```: LAMMPS restart file (frequency determined by "save_freq").
- ```restart_end_run.data```: LAMMPS restart file at the end of simulation.
- ```mc_statistics.out```: contains H4D results such as energy changes, acceptance, etc.

---
# Notes
* The code is not tested for Python versions other than Python3

---
# Dependencies
H4D requires the following packages to run:
* numpy
* scipy
* math
* os
* sys
* random
* ctypes
