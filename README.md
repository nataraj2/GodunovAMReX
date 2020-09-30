# Godunov PPM in amr-wind

This is a python framework for determining the order of accuracy of 
the Godunov piecewise parabolic method implemented in amr-wind. The 
linear advection equation is solved for a Gaussian profile on a domain of 
size (0,1) with the number of grid points ranging from 16 to 512, and the order 
of accuracy is computed.
Just run

```python Godunov.py``` - This runs the advection equation for all grid sizes

```python Accuracy.py``` - This computes and plots the order of accuracy


