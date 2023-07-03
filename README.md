# h4d_lammps

H4D method
==========

H4D, , Hybrid 4D NEMD/MC method in conjunction with LAMMPS

XXX for YYY.

This method is for grand-canonical molecular dynamics simulations.

---
# Reference

Will be added

---
# LAMMPS installation

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
Our H4D method works well with LAMMPS of a version of 27Oct2021.
You can download the particular version here: https://download.lammps.org/tars/index.html.


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
