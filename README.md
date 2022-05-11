Backstepping Tracking Control Using Gaussian Processes with  Event-Triggered Online Learning
=====================================

This code accompanies the paper [1] and implements the code for simulating a one-link planar manipulator with motor dynamics under the control laws proposed in [1] and [2].


[1] J. Jiao, A. Capone, S. Hirche, Backstepping Tracking Control Using Gaussian Processes with  Event-Triggered Online Learning," Under Review for IEEE Control Systems Letters.

[2] A. Capone, S. Hirche, "Backstepping for Partially Unknown Nonlinear Systems Using Gaussian Processes," in IEEE Control Systems Letters, vol. 3, no. 2, pp. 416-421, April 2019, doi: 10.1109/LCSYS.2018.2890467.

Getting started
---------------

This library is tested based on Matlab 2022a together with the GPML toolbox, provided here.

To run the simulation just run the script run_olrt_experiment on Matlab. System parameters and Gaussian process hyperparameters can be changed in the function get_parameters.
