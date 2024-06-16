


DATA_PATH = $1
PH_PATH = $2
OUTPUT_PATH = $3
EXTENTION = 10kb
R2 = 0.5

PPI_PATH = ../../data/gene_bed_files_10kb

echo "----------------Make up GWAS bed file----------------"
python make_up_gwas_bed.py $DATA_PATH
echo "----------------Generating gene snp bed file----------------"
for chr in {1..22}
do
    if [ ! -d "$DATA_PATH/gene_snp_bed_files_$EXTENTION/chr${chr}/" ]; then
        mkdir -p $$DATA_PATH/gene_snp_bed_files_$EXTENTION/chr${chr}/
    fi
    ls $PPI_PATH/chr${chr}/ | parallel -j 32 "bedtools intersect -a $DATA_PATH/gwas.bed -b $PPI_PATH/chr${chr}/{} > $DATA_PATH/gene_snp_bed_files_$EXTENTION/chr${chr}/{}"
done

echo "----------------Generating gene gwas----------------"
for chr in {1..22}
do
    if [ ! -d "$DATA_PATH/gene_gwas_$EXTENTION/chr${chr}/" ]; then
        mkdir -p $DATA_PATH/gene_gwas_$EXTENTION/chr${chr}/
    fi
done
python generate_gene_gwas.py $DATA_PATH $EXTENTION

echo "----------------Generating gene prs----------------"
for chr in {1..22}
do
    echo $chr-$R2
    if [ ! -d "$DATA_PATH/gene_clumps_$EXTENTION/$R2/chr${chr}" ]; then
        mkdir -p $DATA_PATH/gene_clumps_$EXTENTION/$R2/chr${chr}
    fi
    ls $DATA_PATH/gene_gwas_$EXTENTION/chr${chr} | parallel -j 32 "plink \
        --bfile $DATA_PATH/target_data.PH \
        --extract $DATA_PATH//gene_gwas_$EXTENTION/chr${chr}/{} \
        --chr $chr \
        --clump-p1 1 \
        --clump-r2 $R2 \
        --clump-kb 250 \
        --clump $DATA_PATH//gene_gwas_$EXTENTION/chr${chr}/{} $DATA_PATH/gwas.QC.txt \
        --clump-index-first \
        --clump-snp-field SNP \
        --clump-field P \
        --out $DATA_PATH/gene_clumps_$EXTENTION/$R2/chr${chr}/{}" \


    if [ ! -d "$DATA_PATH/gene_prs_$EXTENTION/$R2/chr${chr}" ]; then
            mkdir -p $DATA_PATH/gene_prs_$EXTENTION/$R2/chr${chr}
    fi
            find $DATA_PATH/gene_clumps_$EXTENTION/$R2/chr${chr}/ -type f -name "*clumped" | xargs -I {} basename {} | parallel -j 24 "plink2 \
                --bfile $DATA_PATH/target_data.PH \
                --chr $chr \
                --score $DATA_PATH/gwas.QC.txt 3 4 12 header \
                --q-score-range range_list_full $DATA_PATH/SNP.pvalue \
                --extract $DATA_PATH/gene_clumps_$EXTENTION/$R2/chr${chr}/{} \
                --out $DATA_PATH/gene_prs_$EXTENTION/$R2/chr${chr}/{}"
done
if [ ! -d "$OUTPUT_PATH/$EXTENTION\_$R2" ]; then
    mkdir -p $OUTPUT_PATH/$EXTENTION\_$R2
fi  

python generate_prsnet_data.py $DATA_PATH $PH_PATH $EXTENTION $R2 $OUTPUT_PATH