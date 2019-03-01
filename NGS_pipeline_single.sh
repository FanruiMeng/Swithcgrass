##bwa mapping
##$1:sample name;
##$2:genome file;
##$3:fastq file;
##$4:ploidy
#module load BWA/0.7.17
#module load SAMtools/1.9
#module load picard/2.18.1-Java-1.8.0_152
#module load GATK/3.8-1-0-gf15c1c3ef-Java-1.8.0_112
mkdir $1
cd $1
/mnt/home/mengfanr/software/bwa-0.7.17/bwa mem -R "@RG\tID:$1\tSM:$1\tLB:$1\tPL:Illumina\tPU:$1" -t 8 $2  $3 > 01_$1.sam
##picard reorder
java -jar $EBROOTPICARD/picard.jar ReorderSam I=01_$1.sam O=02_$1_reorder.sam REFERENCE=$2
#samtools sam to bam
samtools view -bS 02_$1_reorder.sam -o 03_$1_reorder.bam
##picard sorted bam
java -jar $EBROOTPICARD/picard.jar  SortSam INPUT=03_$1_reorder.bam OUTPUT=04_$1_sorted.bam SORT_ORDER=coordinate
##mark duplication
#java -jar $PICARD/MarkDuplicates.jar REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=04_$1_sorted.bam OUTPUT=05_$1_dedup.bam METRICS_FILE=05_$1_dedup.metrics

##samtools index last step bam
samtools index 04_$1_sorted.bam

##First variants calling
#java -jar $GATK/GenomeAnalysisTK.jar -glm BOTH -l INFO -R /mnt/home/mengfanr/Switchgrass/00_New_Genome_Data/Panicum_virgatum_var_AP13.mainGenome.fasta -T UnifiedGenotyper -I 05_$1_dedup.bam  -o Variant_calling_first.vcf -metrics Variant_calling_first.metrics -stand_call_conf 40 -stand_emit_conf 10
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T HaplotypeCaller  -I 04_$1_sorted.bam -o Variant_calling_first.vcf -stand_call_conf 40 --max_alternate_alleles 8 -ploidy $4 -nct 8
##seperate snp an indel
java -Xmx2g -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T SelectVariants --variant Variant_calling_first.vcf -o Variant_calling_first_snp.vcf --selectTypeToExclude INDEL
java -Xmx2g -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T SelectVariants --variant Variant_calling_first.vcf -o Variant_calling_first_indel.vcf --selectTypeToExclude SNP

##RealignerTargetCreator
##RealignerTargetCreator first step
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T RealignerTargetCreator -I 04_$1_sorted.bam -o 06_$1.realign.intervals -known Variant_calling_first_indel.vcf

##RealignerTargetCreator second step
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T IndelRealigner -targetIntervals 06_$1.realign.intervals -I 04_$1_sorted.bam -o 06_$1_realn.bam -known  Variant_calling_first_indel.vcf

##Base quality score recalibration
##BQSR first step: BaseRecalibrator
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R $2 -I 06_$1_realn.bam -knownSites Variant_calling_first_indel.vcf -knownSites Variant_calling_first_snp.vcf -o 07_$1.recal.1.grp

##BQSR second step: BaseRecalibrator
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R $2 -I 06_$1_realn.bam -BQSR 07_$1.recal.1.grp -o 07_$1.recal.2.grp -knownSites Variant_calling_first_indel.vcf -knownSites Variant_calling_first_snp.vcf

##BQSR third step: PrintReads
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T PrintReads -R $2 -I 06_$1_realn.bam -BQSR 07_$1.recal.1.grp -o 07_$1_recal.bam

##BQSR fourth step: AnalyzeCovariates
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $2 -before 07_$1.recal.1.grp -after 07_$1.recal.2.grp -csv 07_$1.recal.1.grp.csv -plots 07_$1.recal.1.pdf

##Variants calling Second
#java -jar $GATK/GenomeAnalysisTK.jar -glm BOTH -l INFO -R /mnt/home/mengfanr/Switchgrass/00_New_Genome_Data/Panicum_virgatum_var_AP13.mainGenome.fasta -T UnifiedGenotyper -I 07_$1_recal.bam -D Variant_calling_first_snp.vcf -o Variant_calling_second.vcf -metrics Variant_calling_second.metrics -stand_call_conf 40 -stand_emit_conf 10
#java -jar $GATK/GenomeAnalysisTK.jar -R /mnt/home/mengfanr/Switchgrass/00_New_Genome_Data/Panicum_virgatum_var_AP13.mainGenome.fasta -T HaplotypeCaller  -I 07_$1_recal.bam -o Variant_calling_second.vcf -D Variant_calling_first_snp.vcf -stand_call_conf 40 -stand_emit_conf 10 --max_alternate_alleles 8 -ploidy $2
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T HaplotypeCaller  -I 07_$1_recal.bam -o Variant_calling_second.g.vcf -D Variant_calling_first_snp.vcf -stand_call_conf 40  --max_alternate_alleles 8 -ploidy $4 --emitRefConfidence GVCF

##combine samples
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CombineGVCFs -R /mnt/home/mengfanr/Switchgrass/00_New_Genome_Data/Genome_1.0_version/Pvirgatum_273_v1.0_stitch_Contigs.fa --variant Blackwell.g.vcf --variant Carthage.g.vcf --variant Cave-in-Rock.g.vcf --variant Dacotah.g.vcf --variant ECS-10.g.vcf --variant ECS-11.g.vcf --variant ECS-12.g.vcf --variant ECS-1.g.vcf --variant ECS-2.g.vcf --variant ECS-6.g.vcf --variant High-Tide.g.vcf --variant Kanlow.g.vcf --variant KY1625.g.vcf --variant Pathfinder.g.vcf --variant Shelter.g.vcf --variant Sunburst.g.vcf --variant SW102.g.vcf --variant SW109.g.vcf --variant SW110.g.vcf --variant SW112.g.vcf --variant SW114.g.vcf --variant SW115.g.vcf --variant SW116.g.vcf --variant SW122.g.vcf --variant SW123.g.vcf --variant SW124.g.vcf --variant SW127.g.vcf --variant SW128.g.vcf --variant SW129.g.vcf --variant SW31.g.vcf --variant SW33.g.vcf --variant SW38.g.vcf --variant SW40.g.vcf --variant SW43.g.vcf --variant SW46.g.vcf --variant SW49.g.vcf --variant SW50.g.vcf --variant SW51.g.vcf --variant SW58.g.vcf --variant SW63.g.vcf --variant SW64.g.vcf --variant SW65.g.vcf --variant SW781.g.vcf --variant SW782.g.vcf --variant SW786.g.vcf --variant SW787.g.vcf --variant SW788.g.vcf --variant SW789.g.vcf --variant SW790.g.vcf --variant SW793.g.vcf --variant SW795.g.vcf --variant SW796.g.vcf --variant SW797.g.vcf --variant SW798.g.vcf --variant SW799.g.vcf --variant SW802.g.vcf --variant SW803.g.vcf --variant SW805.g.vcf --variant SW806.g.vcf --variant SW808.g.vcf --variant SW809.g.vcf --variant SWG32.g.vcf --variant SWG39.g.vcf --variant Timber.g.vcf --variant WS4U.g.vcf --variant WS98-SB.g.vcf -o Switchgrass_GBS_Genome_version_1.0_all_groups.g.vcf
##Genotype
java -jar $EBROOTGATK/GenomeAnalysisTK.jar  -T GenotypeGVCFs -R /mnt/home/mengfanr/Switchgrass/00_New_Genome_Data/Genome_1.0_version/Pvirgatum_273_v1.0_stitch_Contigs.fa --variant Switchgrass_GBS_Genome_version_1.0_all_groups.g.vcf -o Switchgrass_GBS_Genome_version_1.0_all_groups_GT.g.vcf --max_num_PL_values 510

##seperate snp an indel

#java -Xmx2g -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T SelectVariants --variant Variant_calling_second.vcf -o Variant_calling_second_snp.vcf --selectTypeToExclude INDEL
#java -Xmx2g -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $2 -T SelectVariants --variant Variant_calling_second.vcf -o Variant_calling_second_indel.vcf --selectTypeToExclude SNP

##Hard filter
#java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $2 -V Variant_calling_second_snp.vcf --filterExpression "QD < 2.0" --filterName "QDfilters" --filterExpression "FS > 60.0" --filterName "FSfilters"   --filterExpression "ReadPosRankSum < -8.0" --filterName "RPRSfilters" --filterExpression "MQ < 40.0" --filterName "MQfilters" --filterExpression "MQRankSum < -12.5" --filterName "MQRSfilters" --clusterSize 3 --clusterWindowSize 10 -o Variant_calling_second_snp_filters.vcf
#java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $2 -V Variant_calling_second_indel.vcf --filterExpression "QD < 2.0" --filterName "QDfilters" --filterExpression "FS > 200.0" --filterName "FSfilters"   --filterExpression "ReadPosRankSum < -20.0" --filterName "RPRSfilters" -o Variant_calling_second_indel_filters.vcf

##Select PASS
#awk -F '\t' '($0~/^#/ || $7=="PASS"){print $0}' Variant_calling_second_snp_filters.vcf>Variant_calling_second_snp_PASS_only.vcf
#awk -F '\t' '($0~/^#/ || $7=="PASS"){print $0}' Variant_calling_second_indel_filters.vcf>Variant_calling_second_indel_PASS_only.vcf



















