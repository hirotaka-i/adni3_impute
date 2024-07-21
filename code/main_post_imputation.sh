source export_bw.sh # download links and password,

# After imputation
mkdir -p impute_eur
cd impute_eur
curl -sL $link_imputation_server | bash
curl -sL $link_imputation_server2 | bash
curl -sL $link_imputation_server3 | bash # download the chromosome files

# check the md5 of the *zip files (Optional)
md5sum impute_eur/*.zip > calculated_checksums.txt
cut -d' ' -f1 impute_eur/results.md5 | sort > sorted_provided_checksums.txt
cut -d' ' -f1 calculated_checksums.txt | sort > sorted_calculated_checksums.txt
diff sorted_provided_checksums.txt sorted_calculated_checksums.txt # should be empty
rm sorted_provided_checksums.txt sorted_calculated_checksums.txt calculated_checksums.txt

# unzip the files (password = Bwvxh9jRGMo4P[)
mkdir -p impute_eur_unzip
for file in impute_eur/*.zip; do unzip -P $password_imputation_server $file -d impute_eur_unzip; done

# Process1 (Standardize the dataset and lift-over hg38)
module load bcftools
module load plink/3.6
wget https://raw.githubusercontent.com/michael-ta/longitudinal-GWAS-pipeline/main/bin/process1.sh 
# Need some modificaton to process1.sh

# Working folder (work2)
mkdir -p work2
cd work2
# get liftOver
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
# hg19ToHg38.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
# bash ../process1.sh 2 ../impute_eur_unzip/chr22.dose.vcf.gz 0.3 hg19 22 chr22
# rewrite in swarm 
for i in {1..22}; do
    echo "bash ../process1.sh 10 ../impute_eur_unzip/chr${i}.dose.vcf.gz 0.3 hg19 $i chr${i}" >> process1.swarm
done
swarm -f process1.swarm -g 200 -t 10 --time=0:30:00 --module bcftools,plink/3.6
# combine all the pfile outputs
for i in {1..22}; do
    echo "chr${i}_p1out" >> merge_list.txt
done
module load plink/6
plink2 --pmerge-list merge_list.txt --make-pgen --out ADNI3_eur_impute_hg38_split
# clean up
rm chr* ADNI3_eur_impute_hg38_split-merge*