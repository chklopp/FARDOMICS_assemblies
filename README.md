# FARDOMICS assembly procedure

This page describes the differents steps used to assemble and check assemblies of *Tecia* *solanivora*, *Cydalima* *perspectalis* and *Leptoglossus* *occidentalis* processed in the frame of the FARDOMICS project. It includes the command lines used. Because all these command lines were ran on the Genotoul Bioinformatic facility platform cluster all software packages used were already install and available as module. Before the command lines you will find the module load commands enabling to use the software package on this infrasctructure. To reproduce these analyses on another cluster the software packages have to be locally present or installed. 
When needed the procedure will be split depending on the available data per species. 

## read checking and genome metrics collection with genomescope2

```
module load bioinfo/jellyfish-2.2.10
module load system/Miniconda3-4.4.10
module load bioinfo/kat-2.4.1

jellyfish count -m 21 -C -s 100M -t 16 -o hifi.jf <(gzip -c -d hifi_reads.fastq.gz) 
jellyfish histo -t 16 -h 10000000 hifi.jf > hifi.histo

module load system/R-3.6.2
module load bioinfo/genomescope2.0-5034ed4

Rscript /usr/local/bioinfo/src/GenomeScope2.0/genomescope2.0-5034ed4/genomescope.R -i hifi.histo  -k 21 -p 2 -o genomescope2_hifi

```

## long read assembly with hifiasm

For *Leptoglossus* *occidentalis* and 
```
module load bioinfo/hifiasm-v0.18.8
hifiasm -t 24 --write-ec  -o hifiasm_0.18.8_Lo_no_HiC  hifi_reads.fastq.gz
```



## assembly kmer content checking with kat

```
module load bioinfo/jellyfish-2.2.10
module load system/Miniconda3-4.4.10
module load bioinfo/kat-2.4.1

jellyfish count -m 21 -C -s 100M -t 16 -o  assembly.jf hifiasm_0.18.8_Lo_no_HiC.bp.p_ctg.gfa.fa
kat comp -t 16 -o kat hifi.jf assembly.jf
kat plot spectra-cn -x 100 -o kat.png -p png kat_comp.mx
```

## assembly protein content checking with busco

```
module load system/Miniconda3-4.7.10
module load bioinfo/busco-5.2.2

busco -c 10 -m geno -i hifiasm_0.18.8_Lo_no_HiC.bp.p_ctg.gfa.fa -o busco_Ts -l insecta_odb10

```


## assembly purging with purge_dups 

```
module load bioinfo/purge_dups-1964aaa
module load bioinfo/minimap2-2.11

run_purge_dups.py -p slurm Lo.json /usr/local/bioinfo/src/Purge_Dups/purge_dups-1964aaa/bin/ purge_dups_out

```

## assembly checking with merqury
```
module load compilers/gcc/12.2.0
module load statistics/R/4.3.1
module load bioinfo/Meryl/1b9294d
module load bioinfo/Merqury/1ef2c36

meryl count threads=16 memory=512 k=31 hifi_reads.fastq.gz output hifi.meryl
merqury.sh hifi.meryl hifiasm_0.18.8_Lo_no_HiC.bp.p_ctg.gfa.fa merqury_out

```

## finding mito genome with MitoHIFI
```
module load containers/singularity/3.9.9
module load bioinfo/souporcell/2.5

singularity exec -B path_to_your_working_dir -B /usr/local/bioinfo/src  MITOHIFI_HOME/mitohifi.sif mitohifi.py -c hifiasm_0.18.8_Lo_no_HiC.bp.p_ctg.gfa.fa  -f Ref_mitochondrion.fasta -g Ref_mitochondrion.gb -t 4

```

## assembly scaffolding with juicer, 3D-dna and juicebox
```
module load bioinfo/samtools/1.18
module load bioinfo/bwa/0.7.17
module load bioinfo/Juicer/1.6

ln -s $SCRIPTS/.
mkdir fastq genome
cd fastq 
ln -s ~/HiC/*.fastq.gz .
cd ../genome 
fold ~/hifiasm_0.18.8_Lo_no_HiC.bp.p_ctg.gfa.fa genome.fa 

samtools faidx genome.fa
cut -f2 genome.fa.fai > genome.tsv # retrieve contig lengths
bwa index genome.fa 
$SCRIPTS/../../misc/generate_site_positions.py Arima "Leptoglossus_occidentalis" genome.fa
cd ..
export workDir=`pwd`
$SCRIPTS/juicer.sh -g Leptoglossus_occidentalis -s Arima -d $workDir -z $workDir/genome/genome.fa -p $workDir/genome/genome.fa.tsv -y $workDir/genome/Leptoglossus_occidentalis_Arima.txt -S early -D $workDir -q workq -l unlimitq -Q 96:00:00 -t 4

mkdir 3D-DNA
cd 3D-DNA

module load  bioinfo/LASTZ/1.04.22 devel/python/Python-3.6.3
module load bioinfo/3D-DNA/529ccf4

run-asm-pipeline.sh -r 0 $PWD/../genome/genome.fa $PWD/../aligned/merged_nodups.txt

```

## assembly annotation with repeatmodeler, repeatmasker, helixer and braker
```

module load bioinfo/RepeatModeler/2.0.4

BuildDatabase -name Scaphoideus_titanus_repeatmodeler -engine ncbi genome.fa
RepeatModeler -threads 32 -database Scaphoideus_titanus_repeatmodeler -engine ncbi 2>&1 | tee repeatmodeler.log

module load devel/python/Python-3.6.3 bioinfo/RepeatMasker/4.1.5

RepeatMasker -pa 12 -html -xsmall -gff -e rmblast -lib Scaphoideus_titanus_repeatmodeler-families.fa genome.fa

module load compilers/gcc/12.2.0 bioinfo/BRAKER/2.0.4
module load bioinfo/Bamtools/2.5.2 bioinfo/samtools/1.14

mkdir $PWD/braker_tmp
braker.pl --species=Scaphoideus_titanus --workingdir=$PWD/braker_tmp --gff3 --cores=12 --genome=genome.fa.masked --bam=RNA-Seq/SRR16109749.bam,RNA-Seq/SRR16109750.bam,RNA-Seq/SRR16109753.bam,RNA-Seq/SRR16109754.bam

singularity run --nv helixer-docker_helixer_v0.3.1_cuda_11.2.0-cudnn8.sif Helixer.py --fasta-path genome.fa --lineage invertebrate --gff-output-path helixer.gff3

agat_sp_complement_annotations.pl --ref braker.gtf --add helixer.gff3 --out complement_annotations.gff3

```

