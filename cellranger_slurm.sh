#!/bin/bash
#SBATCH --job-name=test_slurm       # Job name
#SBATCH --output=test.txt           # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --partition=exacloud    # Partition or queue name
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --mem=64G                     # amount of memory
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=24             # Number of CPU cores per task
#SBATCH --time=7:00:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=hamiltsy@ohsu.edu    # Email address for notifications




/home/exacloud/gscratch/NikolovaLab/software/cellranger-8.0.0/cellranger multi --id GRCh38_multi_test --csv Final_Config.csv --localmem 64 --localcores 24
