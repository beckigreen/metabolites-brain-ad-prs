#============================#
#  6.Derive PGS with PRSice  #
#============================#

filepath=""
cd $filepath

module load apps/R

# APOE region included
Rscript program/PRSice.R \
--dir . \
--prsice program/PRSice_linux \
--base dat/genetics/base/Kunkle_Stage1_post-qc_semi_clumped.txt \
--target dat/genetics/target/NSHD.QC \
--chr CHR \
--bp BP \
--snp  SNP \
--A1 A1 \
--A2 A2 \
--no-regress \
--clump-kb 250 \
--clump-p 1 \
--clump-r2 0.001 \
--bar-levels 5e-8,0.1 \
--print-snp \
--fastscore \
--no-full \
--out analysis/pgs/NSHD_PRS

# APOE region excluded
Rscript program/PRSice.R \
--dir . \
--prsice program/PRSice_linux \
--base dat/genetics/base/Kunkle_Stage1_post-qc_semi_clumped.txt \
--target dat/genetics/target/NSHD.QC_noAPOE  \
--chr CHR \
--bp BP \
--snp  SNP \
--A1 A1 \
--A2 A2 \
--no-regress \
--clump-kb 250 \
--clump-p 1 \
--clump-r2 0.001 \
--bar-levels 5e-8,0.1 \
--print-snp \
--fastscore \
--no-full \
--out analysis/pgs/NSHD_PRS_noAPOE