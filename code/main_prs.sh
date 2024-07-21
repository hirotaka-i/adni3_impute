source export_bw.sh
mkdir -p prs
cd prs

# Download the valiant list
VARIANT_FILE=$RESOURCE_DIR/PRS/rsid_of_interest_chrpos.txt
python ../code/make_bed_from_variant_list.py --variant_list $VARIANT_FILE

# extract variants
plink2 \
--pfile ../work2/ADNI3_eur_impute_hg38_split \
--extract bed1 variants.bed1 \
--var-filter \
--make-pgen \
--rm-dup force-first \
--out variants_eur

python ../code/detect_dup_TYPED_IMPUTED.py --pvar variants_eur.pvar # 0

plink2\
 --pfile variants_eur\
 --exclude variants_eur_exclude.txt\
 --normalize --ref-from-fa force --fa $RESOURCE_DIR/hg38.fa.gz\
 --recode vcf id-delim=^\
 --out variants_eur_norm

bcftools norm -m +any variants_eur_norm.vcf -Ov -o variants_eur_norm_multi.vcf

plink2\
 --vcf variants_eur_norm_multi.vcf\
 --set-all-var-ids chr@:#\
 --sort-vars\
 --vcf-half-call m\
 --make-pgen\
 --out variants_eur_norm_multi

STUDIES=(AD_Bellenguez AD_Kunkle PD_Nalls PD_Nalls_noGBA PD_Nalls_noLRRK2 PD_Nalls_noLRRK2noGBA)
for PRS in ${STUDIES[@]}; do
    SCORE_FILE=$RESOURCE_DIR/PRS/PRS_rsid_${PRS}.txt
    python ../code/adjust_score_file.py --pvar variants_eur_norm_multi.pvar --rsid_ref $VARIANT_FILE --score_file $SCORE_FILE --out variants_eur_norm_multi_${PRS}.score
done

# missing snps are actually not in the imputed data
# ['rs1140239' 'rs1160871' 'rs139643391' 'rs149080927' 'rs35048651'
#  'rs60755019' 'rs616338' 'rs7157106']
# grep "chr16:30010081" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs1140239
# grep "chr7:28129131" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs1160871
# grep "chr2:202878717" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs139643391
# grep "chr19:185425" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs149080927
# grep "17:1728056"   ../work2/ADNI3_eur_impute_hg38_split.pvar # rs35048651
# grep "chr6:41181270" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs60755019
# grep "chr17:49219935" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs616338
# grep "chr14:105761758" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs7157106
# ['rs1160871']
# ['rs9271058']
# grep 'chr6:32607629' ../work2/ADNI3_eur_impute_hg38_split.pvar # rs9271058
# ['rs76763715']
# grep "chr1:155235843" ../work2/ADNI3_eur_impute_hg38_split.pvar # rs76763715


for PRS in ${STUDIES[@]}; do
    SCORE_FILE=variants_eur_norm_multi_${PRS}.score
    plink2 --pfile variants_eur_norm_multi --score $SCORE_FILE list-variants --out variants_eur_norm_multi_${PRS}
done

plink2 --pfile variants_eur_norm_multi --export A include-alt --out variants_eur_norm_multi