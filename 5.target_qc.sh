#====================#
#  5.Target data QC  #
#====================#

filepath=""
cd $filepath

# Load plink
module load apps/plink

# Step 1: Update sex
plink \
--bfile NSHD_QCed \
--update-sex ID_sex.txt \
--make-bed \
--out NSHD.QCsex

# Step 2: Keep only unrelated individuals
plink \
--bfile NSHD.QCsex \
--keep rel.txt \
--make-bed \
--out NSHD.QCsex_rel

# Step 3: Remove individuals with mismatching genetic & self-reported sex
plink \
--bfile NSHD.QCsex_rel \
--remove sexcheck_neuro.txt \
--make-bed \
--out NSHD.QCsex_rel_sex 

plink \
--bfile NSHD.QCsex_rel_sex \
--remove sexcheck_drugdev.txt \
--make-bed \
--out NSHD.QCsex_rel_sex_sex 

# Step 4: ppt & SNP QC (APOE included)
plink \
--bfile NSHD.QCsex_rel_sex_sex \
--maf 0.05 \
--hwe 1e-5 \
--geno 0.02 \
--mind 0.02 \
--make-bed \
--out NSHD.QC

# Step 5: ppt & SNP QC (no APOE)
plink \
--bfile NSHD.QC \
--exclude APOE_snps_NSHD.txt \
--make-bed \
--out NSHD.QC_noAPOE
