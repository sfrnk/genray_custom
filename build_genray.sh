#!/bin/bash
cd mpi
make
cd ..
make -f makefile_mpi_gnu.perl
