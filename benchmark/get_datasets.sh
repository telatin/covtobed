#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

echo "Specify a parameter to also download the exome from 1000 Genomes Project (large)"

set -eo pipefail
command -v samtools || exit 2

mkdir -p "$DIR/../TE"
cd "$DIR/../TE"
echo "Get TE panels"
wget -O "example1.bam" "https://zenodo.org/record/3643386/files/16.bam?download=1"
wget -O "example2.bam" "https://zenodo.org/record/3643386/files/17.bam?download=1"

if  [ "NO$1" != "NO" ]; then
 echo "Get exome"
 wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00258/exome_alignment/HG00258.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram"
 wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00258/exome_alignment/HG00258.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram.crai"
 echo "Convert exome"
 samtools view HG00258.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram -b > HG00258.bam
 samtools index HG00258.bam
 echo "Get genome"
 wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/alignment/HG00114.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai
 wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/alignment/HG00114.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram
 samtools view HG00114.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram -b > genome.bam
 samtools index genome.bam
else
	echo "Skipping exome download"
fi
