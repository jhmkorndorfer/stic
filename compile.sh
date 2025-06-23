export OS=Linux
export CPU=x86_64

export OMPI_CC=gcc
export OMPI_CXX=g++
export OMPI_FC=gfortran

module load netCDF-C++4/4.3.1-gompi-2023a FFTW.MPI/3.3.10-gompi-2023a libtirpc/1.3.3-GCCcore-12.3.0
module load GCCcore/11.3.0
ml Eigen/3.4.0-GCCcore-11.3.0
ml CMake/3.24.3-GCCcore-11.3.0

cd ./src/rh
#make clean
make -j 8

cd rh_1d/
#make clean
make -j 8

cd ../../
#make clean
make -j 8