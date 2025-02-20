module load R
Rscript --vanilla code/get_common_rs.r

mkdir -p work
plink2 --bfile PLINK_set2/ADNI3_PLINK_FINAL_2nd --extract adni3_common_var.txt --make-bed --out work/ADNI3_2nd_set
plink2 --bfile PLINK_set1/ADNI3_PLINK_Final --extract adni3_common_var.txt --make-bed --out work/ADNI3_1st_set
plink --bfile work/ADNI3_1st_set --bmerge work/ADNI3_2nd_set --make-bed --geno 0.01 --out ADNI3_merged

cp ADNI3_merged* work
cd work/
FILENAME=ADNI3_merged
# Remove palindrome SNPs
plink --bfile $FILENAME --write-snplist --out all_snps
awk '{ if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")) print $1 }' all_snps.snplist > palindromic_snps.txt
plink --bfile $FILENAME --exclude palindromic_snps.txt --make-bed --out ${FILENAME}_no_palindromic
FILENAME=${FILENAME}_no_palindromic

# Heterozygosity
plink --bfile $FILENAME --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile $FILENAME --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet
awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > HETEROZYGOSITY_OUTLIERS.txt
cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt
plink --bfile $FILENAME --remove all_outliers.txt --make-bed --out ${FILENAME}_after_heterozyg # no heterozygosity rate > 0.03
FILENAME=${FILENAME}_after_heterozyg
plink --bfile $FILENAME --mind 0.05 --make-bed --out ${FILENAME}_after_call_rate # no missing rate > 0.05
FILENAME=${FILENAME}_after_call_rate

# update variant IDs
plink --bfile $FILENAME --update-name ../adni3_common_ID_rsID_mapping.txt --make-bed --out ${FILENAME}_rename

# sex check skip because of missingness
# Population stratification
FILENAME=${FILENAME}_rename
FILENAME3=premerge
FILENAME4=premerge2
plink --bfile $FILENAME --bmerge ../../../HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile $FILENAME --flip hapmap3_bin_snplis-merge.missnp --make-bed --out $FILENAME3
plink --bfile $FILENAME3 --bmerge ../../../HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile $FILENAME3 --exclude hapmap3_bin_snplis-merge.missnp --out $FILENAME4 --make-bed
plink --bfile $FILENAME4 --bmerge ../../../HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile hapmap3_bin_snplis --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out hapmap3_bin_snplis_pruning
plink --bfile hapmap3_bin_snplis --extract pruning.prune.in --make-bed --out hapmap3_bin_snplis_pruned
plink --bfile hapmap3_bin_snplis_pruned --geno 0.01 --out pca --make-bed --pca 10

# Ancestry comparison
python3 ../code/check_adni_pop.py

# Europeans only
plink --bfile $FILENAME --keep genetic_ancestry_EUROPE.txt --make-bed --out  ${FILENAME}_eur
cat genetic_ancestry_ASIA.txt genetic_ancestry_AFRICA.txt genetic_ancestry_ADMIX.txt > hapmap_outliers.txt

# relatedness
FILENAME=${FILENAME}_eur
module load GCTA
gcta --bfile $FILENAME --make-grm --out GRM_matrix --autosome --maf 0.05 
gcta --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0125 --make-grm
plink --bfile $FILENAME --keep GRM_matrix_0125.grm.id --make-bed --out ${FILENAME}_relatedness


# Missing by haplotype
FILENAME=${FILENAME}_relatedness
plink --bfile $FILENAME --test-mishap --out missing_hap
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\n/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt
plink --bfile $FILENAME --exclude missing_haps_1E4_final.txt --make-bed --out ${FILENAME}_mishap

# Hardy-Weinberg equilibrium
FILENAME=${FILENAME}_mishap
plink --bfile $FILENAME --hwe 1e-4 --out ${FILENAME}_hwe --make-bed




# eur only pcs
FILENAME=${FILENAME}_hwe
plink --bfile $FILENAME --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning_eur
plink --bfile $FILENAME --extract pruning_eur.prune.in --make-bed --out pruned_data_eur
plink --bfile pruned_data_eur --pca 10 --out pca_eur
python ../code/check_adni_pop_eur.py




# Final Prep for MIchigan Imputation Server
plink --bfile ${FILENAME} --freq --out ${FILENAME} 
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
unzip HRC-1000G-check-bim-v4.2.7.zip

perl HRC-1000G-check-bim.pl -b $FILENAME.bim -f $FILENAME.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
sh Run-plink.sh

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  plink --bfile ${FILENAME}-updated-chr$chnum --recode vcf --chr $chnum --out ${FILENAME}$chnum 
done

module load vcftools

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
	vcf-sort ${FILENAME}$chnum.vcf | bgzip -c >  pre_impute_${FILENAME}_$chnum.vcf.gz
done