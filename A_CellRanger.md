## CellRanger

### cellranger mkfastq
module load cellranger/3.0.2
module load bcl2fastq/2.19.1

echo "START"
cellranger mkfastq --id=SLEEPSEQ_FASTQ \
                   --run=SLEEP_SNRNASEQ/00_FROM_CORE/201118_A00672_0044_AHMHJ2DSXY.RGRP \
                   --csv=samplesheet.csv \
                   --ignore-dual-index

echo "DONE"

*SampleSheet.csv*
Lane,Sample,Index
3,A1,SI-GA-A1
3,A2,SI-GA-A2
3,A3,SI-GA-A3
3,A4,SI-GA-A4
3,A5,SI-GA-A5
3,A6,SI-GA-A6
3,A7,SI-GA-A7
3,A8,SI-GA-A8



### fastqc
module load fastqc/0.11.5

FQDIR="SLEEPSEQ_FASTQ"

echo "START"
cd ${FQDIR}

for mylib in `ls`
 do
   cd ${mylib}
   echo ${mylib}
   ls *.fastq.gz | sed "s/_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; fastqc %_001.fastq.gz'
   cd ..
 done

echo "DONE"


### cellranger count
module load cellranger/3.0.2

echo "START"

echo "======================> A1"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A1/*.fastq.gz
cellranger count --id=A1 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A1 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "======================> A2"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A2/*.fastq.gz
cellranger count --id=A2 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A2 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "======================> A3"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A3/*.fastq.gz
cellranger count --id=A3 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A3 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "======================> A4"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A4/*.fastq.gz
cellranger count --id=A4 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A4 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "======================> A5"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A5/*.fastq.gz
cellranger count --id=A5 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A5 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "======================> A6"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A6/*.fastq.gz
cellranger count --id=A6 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A6 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "======================> A7"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A7/*.fastq.gz
cellranger count --id=A7 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A7 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "======================> A8"
fastqc /endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/A8/*.fastq.gz
cellranger count --id=A8 --fastqs=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/SLEEP_SNRNASEQ/01_DEMUX/SLEEPSEQ_DEC2020/outs/fastq_path/HMHJ2DSXY/ --sample=A8 --transcriptome=/endosome/work/Neuroinformatics_Core/s175848/DATA_BACKUP_AK/RESOURCES/MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER --expect-cells=10000

echo "DONE"

