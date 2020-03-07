################################################################################
#
# Copyright (c) 2019, NVIDIA Corporation.  All rights reserved.
#
# Please refer to the NVIDIA end user license agreement (EULA) associated
# with this source code for terms and conditions that govern your use of
# this software. Any use, reproduction, disclosure, or distribution of
# this software and related documentation outside the terms of the EULA
# is strictly prohibited.
#
################################################################################

TESTNO=1
FC=pgfortran
OBJ=o
EXE=out
RUN=
UNAME:=$(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
	OBJ=obj
	EXE=exe
endif

FCFLAGS=-Mcuda
MAIN=t$(TESTNO).CUF
TEST=wmma$(TESTNO)
FILE=wmma$(TESTNO).CUF

ifeq ($(findstring ppc64le, $(UNAME)), ppc64le)
	WMMAUTILS=real2_host_functions.$(OBJ) check_mod.$(OBJ) vector_types.$(OBJ)
	WMMAPP=-DOPENPOWER
else
	WMMAUTILS=real2_host_functions.$(OBJ) check_mod.$(OBJ) vector_types.$(OBJ) ../Utils/avx.s
	WMMAPP=-DX86
endif

all: build run verify

build: $(MAIN) $(FILE)
	$(FC) $(FCFLAGS) -c $(WMMAPP) ../Utils/real2_host_functions.F90
	$(FC) $(FCFLAGS) -c $(WMMAPP) ../Utils/vector_types.F90
	$(FC) $(FCFLAGS) -c ../Utils/check_mod.f90
	$(FC) $(FCFLAGS) -I../Utils -o $(TEST).$(EXE) $(MAIN) $(FILE) $(WMMAUTILS)

run: $(TEST).$(EXE)
	$(RUN) ./$(TEST).$(EXE)

verify:


clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.dwf *.pdb *.mod *.$(OBJ)
