# h4d_lammps

H4D method
==========

H4D, , Hybrid 4D NEMD/MC method in conjunction with LAMMPS

Grand-canonical MD simulation with H4D

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

[L. Belloni, J. Chem. Phys. 151, 021101 (2019)](https://pubs.aip.org/aip/jcp/article/151/2/021101/197798/Non-equilibrium-hybrid-insertion-extraction)

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
The electrode charge time series should be in the same folder and named ```total_charges.out```. An example is provided in ```Example/total_charges.out``` for testing purposes only.

Simply run:

```python3 XXX```

---
# LAMMPS outputs



---
# Notes
* All the output
* Parameters such as the temperature, duration of time steps and the list of desired frequencies are hard-coded, for the moment.
* The code is not tested for Python versions other than Python3

---
# Dependencies
Q2Z requires the following packages to run:
* numpy
* scipy
