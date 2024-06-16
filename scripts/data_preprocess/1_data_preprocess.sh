# !/bin/bash

echo "----------------GWAS PREPROCESSING START----------------"

BASE_DATA_PATH = $1
TARGET_DATA_PATH = $2
PH_PATH = $3
LD_DATA_PATH = $4
OUTPUT_PATH = $5



echo "----------------Standard GWAS QC----------------"
awk 'BEGIN {
       OFS="\t";
       print "CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "OR", "INFO", "MAF", "BETA";
     }
     NR > 1 && ($11 > 0.001  && $10 > 0.3) {
       print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12;
     }' $BASE_DATA_PATH > $OUTPUT_PATH/gwas.a1.txt

echo "----------------SNP ID----------------"
awk 'BEGIN {FS=OFS="\t"} 
     NR==FNR {a[$1,$4]=$2; next} 
     FNR==1 {print; next} 
     ($1,$2) in a { $3=a[$1,$2]; print }' $TARGET_DATA_PATH.bim $OUTPUT_PATH/gwas.a1.txt > $OUTPUT_PATH/gwas.a2.txt

echo "----------------Removing Duplicate SNPs----------------"
awk '{seen[$3]++; if(seen[$3]==1){ print}}' $OUTPUT_PATH/gwas.a2.txt > $OUTPUT_PATH/gwas.a3.txt

echo "----------------Removing Ambiguous SNPs----------------"
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' $OUTPUT_PATH/gwas.a3.txt \
        > $OUTPUT_PATH/gwas.a4.txt

echo "----------------Removing Mismatch----------------"
Rscript mismatch.R $TARGET_DATA_PATH $OUTPUT_PATH/gwas.a4.txt 'target_data' $OUTPUT_PATH/
Rscript mismatch.R $LD_DATA_PATH $OUTPUT_PATH/gwas.a4.txt 'ld_data' $OUTPUT_PATH/
awk 'NR == FNR { mismatched[$1]; next } FNR == 1 || !($3 in mismatched)' $OUTPUT_PATH/target_data.tofilter.snplist $OUTPUT_PATH/gwas.a4.txt > $OUTPUT_PATH/gwas.a5.txt
awk 'NR == FNR { mismatched[$1]; next } FNR == 1 || !($3 in mismatched)' $OUTPUT_PATH/ld_data.tofilter.snplist $OUTPUT_PATH/gwas.a5.txt > $OUTPUT_PATH/gwas.a6.txt

echo "----------------Test missing----------------"
plink --bfile $TARGET_DATA_PATH --pheno $PH_PATH --test-missing midp --pfilter 1e-5 --out $OUTPUT_PATH/SNPLIST.MISSDIFF --memory 65536
tail -n +2 $OUTPUT_PATH/SNPLIST.MISSDIFF.missing | sed s/^\ *//g | tr -s ' ' ' ' | cut -f 2 -d ' ' > $OUTPUT_PATH/MISSDIFF_to_filter.txt
echo "----------------Update GWAS----------------"
awk 'BEGIN {FS=OFS="\t"} 
     NR==FNR {exclude[$1] = 1; next} 
     FNR==1 || !($3 in exclude) {print}' $OUTPUT_PATH/MISSDIFF_to_filter.txt $OUTPUT_PATH/gwas.a6.txt >  $OUTPUT_PATH/gwas.a7.txt
awk '{print $3}' $OUTPUT_PATH/gwas.a7.txt > $OUTPUT_PATH/gwas.snplist


echo "----------------Make LD data bed----------------"
plink \
    --bfile $LD_DATA_PATH \
    --out $OUTPUT_PATH/ld_data.SNP \
    --extract $OUTPUT_PATH/gwas.snplist \
    --a1-allele $OUTPUT_PATH/ld_data.a1 \
    --memory 65536 \
    --make-just-bim
echo "----------------Make Target data bed----------------"
plink \
    --bfile $TARGET_DATA_PATH \
    --out $OUTPUT_PATH/target_data.SNP \
    --extract $OUTPUT_PATH/gwas.snplist \
    --a1-allele $OUTPUT_PATH/target_data.a1 \
    --memory 65536 \
    --make-just-bim 

echo "----------------Final Construction----------------"
awk 'FNR==NR{snps[$2]; next} $2 in snps {print $2}' $OUTPUT_PATH/ld_data.SNP.bim $OUTPUT_PATH/target_data.SNP.bim > $OUTPUT_PATH/common_snps.txt
comm -12 <(sort $OUTPUT_PATH/common_snps.txt) <(sort $OUTPUT_PATH/gwas.snplist) > $OUTPUT_PATH/snplist

plink \
    --bfile $LD_DATA_PATH \
    --make-bed \
    --out $OUTPUT_PATH/ld_data.PH \
    --extract $OUTPUT_PATH/snplist \
    --a1-allele $OUTPUT_PATH/ld_data.a1 \
    --memory 65536 

plink \
    --bfile $TARGET_DATA_PATH \
    --make-bed \
    --out $OUTPUT_PATH/target_data.PH \
    --extract $OUTPUT_PATH/snplist \
    --a1-allele $OUTPUT_PATH/target_data.a1 \
    --memory 65536 
awk 'FNR==NR{snps[$1]; next} FNR==1 || $3 in snps' $OUTPUT_PATH/snplist $OUTPUT_PATH/gwas.a7.txt > $OUTPUT_PATH/gwas.QC.txt
awk '{print $3,$8}' $OUTPUT_PATH/gwas.QC.txt > $OUTPUT_PATH/SNP.pvalue
gzip -c $OUTPUT_PATH/gwas.QC.txt > $OUTPUT_PATH/out/gwas.QC.gz