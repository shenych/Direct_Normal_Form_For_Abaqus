This program is to generate a nonlinear reduced order model (ROM) from a full-order finite element
model by applying the DNF approach in ABAQUS. The ROM can be used to compute nonlinear dynamics, e.g. time
response, backbone curves, and frequency response functions. Most importantly, thanks to invariant manifold theory, the
ROM is able to accurately capture the nonlinear properties and characteristics with only a small approximation error as
compared to the full-order model.

One should install MATLAB and ABAQUS in order to launch the program. The program is able to generate the
ROM from any mesh files in *.inp format which ABAQUS can read. The program is divided into two parts, one for the
Main code.m where to launch the full program and change input parameters, and the folder SRC DNF including all the
functions.
