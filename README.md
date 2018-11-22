PyORB (Python Optimized Reduced Basis) is a set of Python scripts for the use of reduced basis method. The code is designed to rely on a third party library for the Full Order Model (FOM) approoximation of the PDE, which can be either implemented in MATLAB or C (C++). In order to employ any third party library (TPL), an external engine is implemented in pyorb_core/tpl_managers. Such engine requires that the actual implementations of functions required for the offline phase of the RB method are found in the TPL which is provided in the main through a path.

Current tests included in the library in the folder example/ take into account either the MATLAB TPL feamat or the C++ TPL LifeV.

