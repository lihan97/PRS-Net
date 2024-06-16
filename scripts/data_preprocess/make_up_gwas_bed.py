import sys
if __name__ == '__main__':
    ROOT_PATH = sys.argv[1]
    first_line=True
    with open(f"{ROOT_PATH}/gwas.bed", "w") as f:
        for line in open(f'{ROOT_PATH}/gwas.QC.txt'):
            if first_line:
                first_line = False
                continue
            line = line.strip().split()
            snpid = line[2]

            chrom = 'chr'+line[0]
            if line[1] == 'NA':
                continue
            pos = int(float(line[1]))
            f.write('\t'.join([chrom, str(pos-1), str(pos), snpid])+'\n')


