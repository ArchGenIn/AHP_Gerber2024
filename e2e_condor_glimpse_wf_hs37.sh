#!/usr/bin/env bash

#threadid==$1; bamlist.txt==$2 output_directory==$3

samplerow=`expr $1 + 1`

raw_bamfile=`sed -n $samplerow"p" $2` #

samplename=`basename $(sed -n $samplerow"p" $2 | cut -d'.' -f1 )`

endir=$3

#

REFGENOME='/mnt/agi_gluster_volume/GLIMPSE/REFGENOMES/hs37d5.fa.gz'

REFPANEL_DIR='/mnt/agi_gluster_volume/GLIMPSE/ALI_REFPANEL/1KG_30X_hg19/'

RECMAP_DIR='/mnt/agi_gluster_volume/GLIMPSE/MAP/'

#

CHUNK_EXEC='/mnt/agi_gluster_volume/GLIMPSE/static_bins/GLIMPSE_chunk_static'

PHASE_EXEC='/mnt/agi_gluster_volume/GLIMPSE/static_bins/GLIMPSE_phase_static'

LIGATE_EXEC='/mnt/agi_gluster_volume/GLIMPSE/static_bins/GLIMPSE_ligate_static'

SAMPLE_EXEC='/mnt/agi_gluster_volume/GLIMPSE/static_bins/GLIMPSE_sample_static'


#

WORKDIR=$endir'/'$samplename'/'

mkdir -p $WORKDIR

echo $samplename > $WORKDIR'/'$samplename'.temp'

trim_bamfile=$WORKDIR'/'$samplename'.trim_3_bp.bam'

bam trimBam $raw_bamfile $trim_bamfile 3

samtools index $trim_bamfile

#

for ch in {1..22}; do
	bcftools mpileup -q 20 -Q 20 --ignore-RG -f $REFGENOME \
	 -E -a 'FORMAT/AD' -T $REFPANEL_DIR'/chr'$ch'.sites.tsv.gz' \
	 -r $ch $trim_bamfile -Ou \
	 | bcftools reheader -s $WORKDIR'/'$samplename'.temp' \
	 | bcftools norm -m+any -f $REFGENOME \
	 -Ob -o $WORKDIR'/'$samplename'.chr'$ch'.mpileup.bcf';
	#
	bcftools index -f $WORKDIR'/'$samplename'.chr'$ch'.mpileup.bcf';
	#
	bcftools view $WORKDIR'/'$samplename'.chr'$ch'.mpileup.bcf' \
	 -T $REFPANEL_DIR'/chr'$ch'.sites.tsv.gz' -Ou \
	 | bcftools call -Aim -C alleles \
	 -T $REFPANEL_DIR'/chr'$ch'.sites.tsv.gz' \
	 -Ob -o $WORKDIR'/'$samplename'.chr'$ch'.call.snps.bcf';
	#
	bcftools index -f $WORKDIR'/'$samplename'.chr'$ch'.call.snps.bcf';
	#
	$CHUNK_EXEC \
	 --input $WORKDIR'/'$samplename'.chr'$ch'.call.snps.bcf' \
	 --region $ch --window-size 2000000 --buffer-size 200000 \
	 --output $WORKDIR'/'$samplename'.chr'$ch'.chunk.txt';
	#
	#
	mkdir -p $WORKDIR'/CHR'$ch'/';
	while IFS="" read -r LINE || [ -n "$LINE" ];
		do
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3);
	ORG=$(echo $LINE | cut -d" " -f4);
	#IMPUTED_OUT=$WORKDIR'/'$samplename'.imputed.'$ch'.'${ID}'.vcf.gz';
	IMPUTED_OUT=$WORKDIR'/CHR'$ch'/'$samplename'.imputed.'$ch'.'${ID}'.bcf';
	$PHASE_EXEC \
	 --input $WORKDIR'/'$samplename'.chr'$ch'.call.snps.bcf' \
	 --reference $REFPANEL_DIR'/chr'$ch'.bcf' \
	 --map $RECMAP_DIR'/chr'$ch'.b37.gmap.gz' \
	 --input-region ${IRG} --output-region ${ORG} \
	 --output $IMPUTED_OUT
		done < $WORKDIR'/'$samplename'.chr'$ch'.chunk.txt';
	#
	#
	IMPUTED_VCF_LIST=$WORKDIR'/list.chr'$ch'.txt';
	ls -1 $WORKDIR'/CHR'$ch'/'*'.bcf' > ${IMPUTED_VCF_LIST};
	LIGATED_OUT=$WORKDIR'/'$samplename'.chr'$ch'.ligated.bcf';
	$LIGATE_EXEC --input ${IMPUTED_VCF_LIST} --output ${LIGATED_OUT};
	#
	#
	PHASED_OUT=$WORKDIR'/'$samplename'.chr'$ch'.ligated.phased.bcf';
	$SAMPLE_EXEC --input ${LIGATED_OUT} --solve --output ${PHASED_OUT};
	#
	#
	#-c FORMAT/GP --rename-chrs /mnt/md0/GLIMPSE/chr_name_conv_2nochr.txt -Oz
	GP_OUT=$WORKDIR'/'$samplename'.chr'$ch'.ligated.phased.GP.bcf'
	bcftools annotate -a $LIGATED_OUT -c FORMAT/GP \
 	 -Oz -o $GP_OUT $PHASED_OUT;
 	bcftools index -f $GP_OUT;
done

