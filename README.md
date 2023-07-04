# Grand-canonical MD simulation with H4D method
==========

H4D, Hybrid 4D NEMD/MC method in conjunction with LAMMPS



Author        : Jeongmin Kim, a former postdoc at PHENIX with Benjamin Rotenberg

Contact       : jeongmin0658 AT gmail DOT com


The folder includes all the sources and example folders:

Each folder has its own readme file too.

"H4D"         : LAMMPS c++ and head files

"ver0.0"      : H4D Python scripts

"examples"    : Example systems, including 

              - "LJ_salt100"  : Example of LJ electrolytes   
              
              - "JC_SPCE_1m"  : Example of aqueous NaCl electrolytes 

---
# References

[J. Kim, L. Belloni, B. Rotenberg, Arxiv](https://arxiv.org)

[L. Belloni, J. Chem. Phys. 151, 021101 (2019)](https://pubs.aip.org/aip/jcp/article/151/2/021101/197798/Non-equilibrium-hybrid-insertion-extraction) - original H4D method

---
# LAMMPS installation
Our H4D method works well with LAMMPS of a version of 27 Oct 2021.
You can download [lammps-27Oct2021.tar.gz](https://download.lammps.org/tars/index.html).

---
# H4D parameters

* Put the location of source codes.
* XXX
* XXX



---
# H4D run
Before run LAMMPS with H4D, you must define where the "ver0.0" folder is located. Use ```setup.sh```.

Simply run:

```python3 run_gcmc_sample.py```
or in case you use multiple cpus, 

```mpirun -np number_of_cpus python3 run_gcmc_sample.py```

---
# LAMMPS outputs
All the outputs can be moved to the "output" folder. Use ```copy_output.sh```.


---
# Notes
* All the output
* Parameters such as the temperature, duration of time steps and the list of desired frequencies are hard-coded, for the moment.
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
