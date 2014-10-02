#!/usr/bin/env python

import os
import re
import string
import tempfile
import subprocess
import argparse
import shutil
from glob import glob

def scan_workdir(base):	
	fastq_files = glob(os.path.join(base, "*.fastq"))
	if len(fastq_files):
		o = {}
		for i in fastq_files:
			basename = re.sub(r'_[12].fastq$', '', i)
			o[basename] = o.get(basename, 0) + 1
		if not all( (i == 2 for i in o.values())):
			raise Exception("Missing Pair")
		return ( 'cat', list( (os.path.basename(i), "%s_1.fastq" % i,"%s_2.fastq" % i) for i in o.keys())) 
	
	fastq_gz_files = glob(os.path.join(base, "*.fastq.gz"))
	if len(fastq_gz_files):
		o = {}
		for i in fastq_gz_files:
			basename = re.sub(r'_[12].fastq.gz$', '', i)
			o[basename] = o.get(basename, 0) + 1
		if not all( (i == 2 for i in o.values())):
			raise Exception("Missing Pair")
		return ( 'zcat', list( (os.path.basename(i), "%s_1.fastq.gz" % i,"%s_2.fastq.gz" % i) for i in o.keys())) 

	raise Exception("Unable to determine input type")

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("--out", default="out.bam")
	parser.add_argument("--genomeDir")
	parser.add_argument("--tarFileIn", default=None)
	parser.add_argument("--dirIn", default=None)
	parser.add_argument("--workDir", default="./")
	parser.add_argument("--runThreadN", type=int, default=8)
	parser.add_argument("--outFilterMultimapScoreRange", type=int, default=1 )
	parser.add_argument("--outFilterMultimapNmax", type=int, default=20 )
	parser.add_argument("--outFilterMismatchNmax", type=int, default=10 )
	parser.add_argument("--alignIntronMax", type=int, default=500000 )
	parser.add_argument("--alignMatesGapMax", type=int, default=1000000 )
	parser.add_argument("--sjdbScore", type=int, default=2 )
	parser.add_argument("--alignSJDBoverhangMin", type=int, default=1)
	parser.add_argument("--genomeLoad", default="NoSharedMemory")
	parser.add_argument("--outFilterMatchNminOverLread", type=float, default=0.33)
	parser.add_argument("--outFilterScoreMinOverLread", type=float, default= 0.33)
	parser.add_argument("--outSAMstrandField", default="intronMotif" )
	parser.add_argument("--outSAMattributes", default=["NH", "HI", "NM", "MD", "AS", "XS"] )
	parser.add_argument("--outSAMunmapped", default="Within")
	parser.add_argument("--outSAMtype", default=["BAM", "Unsorted"])
	
	args = parser.parse_args()
	
	created_input_dir = False
	if args.tarFileIn is not None:
		workdir = tempfile.mkdtemp(dir=args.workDir, prefix="star_inputdir_")
		if args.tarFileIn.endswith(".gz"):
			tarcmd = "tar xvzf %s -C %s" % (args.tarFileIn, workdir)
		elif args.tarFileIn.endswith(".tar"):
			tarcmd = "tar xvf %s -C %s" % (args.tarFileIn, workdir)
		subprocess.check_call(tarcmd, shell=True)
		args.dirIn = os.path.abspath(workdir)
		created_input_dir = True
	
	align_sets = scan_workdir(args.dirIn)
	
	align_template_str = """STAR \
--genomeDir ${genomeDir} --readFilesIn ${fastq_left} ${fastq_right} \
--runThreadN ${runThreadN} \
--outFilterMultimapScoreRange ${outFilterMultimapScoreRange} \
--outFilterMultimapNmax ${outFilterMultimapNmax} \
--outFilterMismatchNmax ${outFilterMismatchNmax} \
--alignIntronMax ${alignIntronMax} \
--alignMatesGapMax ${alignMatesGapMax} \
--sjdbScore ${sjdbScore} \
--alignSJDBoverhangMin ${alignSJDBoverhangMin} \
--genomeLoad ${genomeLoad} \
--readFilesCommand ${readFilesCommand} \
--outFilterMatchNminOverLread ${outFilterMatchNminOverLread} \
--outFilterScoreMinOverLread ${outFilterScoreMinOverLread} \
--outSAMstrandField ${outSAMstrandField} \
--outSAMattributes ${outSAMattributes} \
--outSAMunmapped ${outSAMunmapped} \
--outSAMtype ${outSAMtype}
"""

	readhead_template_str="samtools view -H ${align_dir}/Aligned.out.bam > ${align_dir}/header.sam"

	#sort_template_str="samtools rehead ${align_dir}/header.sam ${align_dir}/Aligned.out.bam | samtools samtools sort -@ 8 - ${align_dir}/Aligned.out.sorted"
	sort_template_str="samtools sort -@ 8 ${align_dir}/Aligned.out.bam ${align_dir}/Aligned.out.sorted"

	
	out_dirs = []
	for pair in align_sets[1]:
		cmd = string.Template(align_template_str).substitute({
			'genomeDir' : os.path.abspath(args.genomeDir),
			'fastq_left' : os.path.abspath(pair[1]),
			'fastq_right' : os.path.abspath(pair[2]),
			'runThreadN' : args.runThreadN,
			'outFilterMultimapScoreRange' : args.outFilterMultimapScoreRange,
			'outFilterMultimapNmax' : args.outFilterMultimapNmax,
			'outFilterMismatchNmax' : args.outFilterMismatchNmax,
			'alignIntronMax' : args.alignIntronMax,
			'alignMatesGapMax': args.alignMatesGapMax,
			'sjdbScore': args.sjdbScore,
			'alignSJDBoverhangMin' : args.alignSJDBoverhangMin,
			'genomeLoad' : args.genomeLoad,
			'readFilesCommand' : align_sets[0],
			'outFilterMatchNminOverLread' : args.outFilterMatchNminOverLread,
			'outFilterScoreMinOverLread' : args.outFilterScoreMinOverLread,
			'outSAMstrandField' : args.outSAMstrandField,
			'outSAMattributes' : " ".join(args.outSAMattributes),
			'outSAMunmapped' : args.outSAMunmapped, 
			'outSAMtype' : " ".join(args.outSAMtype)
		})
		align_dir = os.path.abspath( tempfile.mkdtemp(dir=args.workDir, prefix="star_aligndir_") )
		print "Running", cmd
		subprocess.check_call(cmd, shell=True, cwd=align_dir)

		"""
		cmd = string.Template(readhead_template_str).substitute({
			'align_dir' : align_dir
		})
		print "Running", cmd		
		subprocess.check_call(cmd, shell=True, cwd=align_dir)
		rg_line = "@RG\tID:%s\tSM:%s\n" % (pair[0], pair[0])		
		with open( "%s/header.sam" % (align_dir), "a") as handle:
			handle.write(rg_line)
		"""		
		cmd = string.Template(sort_template_str).substitute({
			'align_dir' : align_dir
		})
		print "Running", cmd
		subprocess.check_call(cmd, shell=True, cwd=align_dir)

		out_dirs.append(align_dir)
	
	cmd = "samtools merge -r -@ 8 %s %s" % (os.path.abspath(args.out), " ".join(out_dirs))
	print "Running", cmd
	subprocess.check_call(cmd, shell=True, cwd=align_dir)
		
	#if created_input_dir:
	#	shutil.rmtree(args.dirIn)
	
