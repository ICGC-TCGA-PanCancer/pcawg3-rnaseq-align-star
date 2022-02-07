# ICGC PCAWG Scripts for RNA-Seq alignments using STAR

All data pre-processing and interface choices for the STAR aligner are wrapped in the `star_align.py` wrapper function. This wrapper takes care of checking ICGC naming conventions and to gather the relevant read group information from the metadata. 

For generating the RNA-Seq alignments of the ICGC PCAWG project, the wrapper was called on each of the samples (which were provided as one tar-ball per analysis ID)using the following setup:

```
genome=/path/to/hg19_hs37d5_G19.overhang100_STAR    ### the STAR index built on hg19_hs37d5_G19 with an sj over of 100
genomeFasta=/path/to/hg19_hs37d5_G19.fa             ### Fasta version of the indexed genome
input=/path/to/input.tar.gz                         ### input tarball containing reads to be aligned
outfile=/path/to/output                             ### output prefix
metadata=/path/to/metadata_table                    ### metadata table
aaid=1234-5678-9012                                 ### analysis ID of current run
export TMPDIR=/path/to/temp                         ### path to folder to store temp data


python star_align.py --useTMP TMPDIR \
                     --keepJunctions \
                     --twopass1readsN -1 \
                     --genomeDir $genome \
                     --tarFileIn ${input} \
                     --analysisID $aaid \
                     --metaData $metadata \
                     --out $outfile \
                     --genomeFastaFiles $genomeFasta \
                     --limitBAMsortRAM 70000000000
```
