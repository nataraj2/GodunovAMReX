# Godunov PPM in amr-wind

This is a python framework for determining the order of accuracy of 
the Godunov piecewise parabolic method implemented in amr-wind. The 
linear advection equation is solved for a Gaussian profile on a domain of 
size (0,1) with the number of grid points ranging from 16 to 512, and the order 
of accuracy is computed.
Just run

```python Godunov.py_constantCFL``` - This runs the advection equation for all grid sizes

```python Accuracy.py_constCFL``` - This computes and plots the order of accuracy

For constant dt, use the scripts appended with `constdt` for the Godunov and Accuracy.

A numerical spectral analysis code is also available which plots the transfer function of the scheme 
which shows the resolving power of the scheme.


