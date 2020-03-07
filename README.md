# PGIexample
examples from PGI installation folder



     Copyright (c) 2018, NVIDIA CORPORATION.  All rights reserved.

NVIDIA CORPORATION and its licensors retain all intellectual property
and proprietary rights in and to this software, related documentation
and any modifications thereto.  Any use, reproduction, disclosure or
distribution of this software and related documentation without an 
express license agreement from NVIDIA CORPORATION is strictly 
prohibited.


This directory contains all the example programs which are included in 
the PGI package.  The examples are organized by feature, with examples 
for each feature contained in a seperate subdirectory. 

Following is a list of the files in this directory:

    README             This file

    AutoPar            Directory containing auto parallelization examples

    CUDA-Fortran       Directory containing CUDA Fortran examples

    F2003              Directory containing Fortran 2003 examples

    MPI                Directory containing MPI examples

    OpenACC            Directory containing OpenACC examples

    OpenMP             Directory containing OpenMP examples


 
To build and run the examples within one of these subdirectories using 
the 64-bit compilers for example, do the following:

    % cd /my/work/dir
    % cp -r $PGI/win64/2018/examples/AutoPar .
    % cd AutoPar

Now build/run all the tests within the AutoPar subdirectory using the 
following command.

    % make NTHREADS=4 all

To get further instructions about building/running the tests within one
of the subdirectories type "make" without any arguments.

