#! /bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH --job-name="cleanup"
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=500Mb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu


cd /projects/b1052/mckenna/wssc_meta/results

rm -r ./megahit/*/intermediate_contigs
rm -r ./metabat_checkm/*/bins
rm -r ./metabat_checkm/*/storage