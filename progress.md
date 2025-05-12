# Daily Progress Log

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
    - Or a in single line... ml netCDF-C++4/4.3.1-gompi-2023a FFTW.MPI/3.3.10-gompi-2023a libtirpc/1.3.3-GCCcore-12.3.0 GCCcore/11.3.0
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
    # module load Python
    
    # Set the batchhost variable
    # export BATCHHOST=$(cnode03)
    
    # module load Anaconda3
    # eval "$(conda shell.bash hook)"
    # conda activate stic
    
    #cp input_RF.cfg input.cfg
    
    mpiexec ../src/STiC.x
    # srun ../src/STiC.x

### Challenges
- NA

### Learnings
- NA

### Next Steps
- Install it on miniHPC
- Trace stic

---
