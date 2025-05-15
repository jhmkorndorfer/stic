# Daily Progress Log

## Date: 15/05/2025

### Accomplishments
- [x] Task 1: create new branch and rsynch files with Jonas Z stic version. 


### Challenges
- NA

### Learnings
- Interestinly JZ version of STiC is VERY different from the newest version available on the STiC reposity. We will have to investigate this further another day.

---


## Date: 14/05/2025

### Accomplishments
- [x] Task 1: Run stic with problematic input provided by Jonas Z.
    - It worked while we expected it to fail!!!!!!!


### Challenges
- STiC seems to be producing different outputs for the same input files depending on how it was compiled...

### Learnings
- NA

---

## Date: 13/05/2025

### REMEMBER, CHECK PREVIOUS COMMITS FOR THE SMALL CHANGES IN MAKEFILES TAKEN FROM Jonas

### Accomplishments
- [x] Task 1: Fix compilation issues.
- [x] Task 2: Finally compile stic on UBELIX.
- [x] Task 3: Trying to prepare input for stic.
    - **Using ml Python/3.11.3-GCCcore-12.3.0 and virtual env namely sticenv. python3 -m venv sticenv. This is inside example folder of stic.**
        - source sticenv/bin/activate
        - Requirements: numpy, Matplotlib, scipy, netCDF4, astropy
- [x] Task 4: Test run stic on UBELIX.
- [x] Task 5: Update stic README to include some of the changes required to compile it which are also described here.


### Challenges
- Error in the first make for rh:
    - g++ -I../ -I../ -I../eigen3 -O3 -march=native -I/usr/include/tirpc -Ieigen3/ -std=c++11 -g3 -c  readAtomFile.cc -o readAtomFile.o g++ -I../ -I../ -I../eigen3 -O3 -march=native -I/usr/include/tirpc -Ieigen3/ -std=c++11 -g3 -c  solveLinearCXX.cc -o solveLinearCXX.o solveLinearCXX.cc:3:10: fatal error: eigen3/Eigen Dense: No such file or directory 3 | #include <eigen3/Eigen/Dense>
- Needed to figure it out why the includes were wrong... Turns out the old code was including some libraries locally, assuming that the system modules were not available and requiring additional local installations. Fix described below (Learnings).
- Errors in the last make for STiC itself:
    - clm.cc:51:10: fatal error: eigen3/Eigen/Dense: No such file or directory
   51 | #include <eigen3/Eigen/Dense> compilation terminated.
    - clm.h:42:10: fatal error: eigen3/Eigen/Dense: No such file or directory
   42 | #include <eigen3/Eigen/Dense>


### Learnings
- To fix the issues above: 
    - **ml Eigen/3.4.0-GCCcore-11.3.0**
    - Modify the include statement on **solveLinearCXX.cc:3:1, clm.cc:51:10, clm.h:42:10** to avoid including from a local compilation of Eigen and use the actual module. Otherwise one needs to install Eigen Localy.

---

## Date: 12/05/2025

### Accomplishments
- [x] Task 1: Forking repository.
- [x] Task 2: Copy everything to UBELIX sever.
- [x] Task 3: Copy everything to miniHPC sever.
- [x] Task 4: Install on UBELIX. Modules:
    - ml netCDF-C++4/4.3.1-gompi-2023a
    - ml FFTW.MPI/3.3.10-gompi-2023a
    - ml libtirpc/1.3.3-GCCcore-12.3.0
    - ml GCCcore/11.3.0
    - ml Eigen/3.4.0-GCCcore-11.3.0
- Example bash script to run STIC
```bash
#!/bin/bash
#SBATCH --job-name="STIC FRAME0 first cycle recalibrated first test map"
#SBATCH --time=6:00:00
#SBATCH --partition=epyc2
# SBATCH --qos=job_icpu-aiub
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1000
#SBATCH --cpus-per-task=1
##SBATCH --mail-user=
#SBATCH --mail-type=end,fail
#SBATCH --verbose
# SBATCH --nodes=20
# SBATCH --nodelist=bnode[001-011]
# SBATCH --exclude=bnode[001-028] # bnode002,bnode003,bnode004,bnode005,bnode006,bnode007,bnode008,bnode009,bnode010,bnode011
# Your code below this line
# HPC_WORKSPACE=aiub_sml_ws module load Workspace
module load netCDF-C++4/4.3.1-gompi-2023a FFTW.MPI/3.3.10-gompi-2023a libtirpc/1.3.3-GCCcore-12.3.0
module load GCCcore/11.3.0
ml Eigen/3.4.0-GCCcore-11.3.0
# module load Python

# Set the batchhost variable
# export BATCHHOST=$(cnode03)

# module load Anaconda3
# eval "$(conda shell.bash hook)"
# conda activate stic

#cp input_RF.cfg input.cfg

mpiexec ../src/STiC.x
# srun ../src/STiC.x
```

### Challenges
- NA

### Learnings
- NA

---
