#!/usr/bin/env python

import os
import sys
import re
import string
import tempfile
import subprocess
import argparse
import shutil
import lxml.etree as etree
from glob import glob

def scan_workdir(base): 

    ### scan for paired-end files
    #############################

    ### unzipped input
    fastq_files = glob(os.path.join(base, "*_[12].fastq"))
    if len(fastq_files):
        o = {}
        for i in fastq_files:
            basename = re.sub(r'_[12].fastq$', '', i)
            o[basename] = o.get(basename, 0) + 1
        if not all( (i == 2 for i in o.values())):
            raise Exception("Missing Pair")
        return ( 'cat', list( (os.path.basename(i), "%s_1.fastq" % i,"%s_2.fastq" % i) for i in o.keys()), 'PE') 
    
    ### gzipped input
    fastq_gz_files = glob(os.path.join(base, "*_[12].fastq.gz"))
    if len(fastq_gz_files):
        o = {}
        for i in fastq_gz_files:
            basename = re.sub(r'_[12].fastq.gz$', '', i)
            o[basename] = o.get(basename, 0) + 1
        if not all( (i == 2 for i in o.values())):
            raise Exception("Missing Pair")
        return ( 'zcat', list( (os.path.basename(i), "%s_1.fastq.gz" % i,"%s_2.fastq.gz" % i) for i in o.keys()), 'PE') 

    ### bzipped input
    fastq_bz_files = glob(os.path.join(base, "*_[12].fastq.bz"))
    if len(fastq_gz_files):
        o = {}
        for i in fastq_gz_files:
            basename = re.sub(r'_[12].fastq.bz$', '', i)
            o[basename] = o.get(basename, 0) + 1
        if not all( (i == 2 for i in o.values())):
            raise Exception("Missing Pair")
        return ( 'bzcat', list( (os.path.basename(i), "%s_1.fastq.bz" % i,"%s_2.fastq.bz" % i) for i in o.keys()), 'PE') 

    ### scan for single-end files
    #############################

    ### unzipped input
    fastq_files = glob(os.path.join(base, "*.fastq"))
    if len(fastq_files):
        return ( 'cat', list( (os.path.basename(re.sub(r'.fastq$', '', i)), i) for i in fastq_files), 'SE') 

    ### gzipped input
    fastq_files = glob(os.path.join(base, "*.fastq.gz"))
    if len(fastq_files):
        return ( 'zcat', list( (os.path.basename(re.sub(r'.fastq.gz$', '', i)), i) for i in fastq_files), 'SE') 

    ### bzipped input
    fastq_files = glob(os.path.join(base, "*.fastq.bz"))
    if len(fastq_files):
        return ( 'bzcat', list( (os.path.basename(re.sub(r'.fastq.bz$', '', i)), i) for i in fastq_files), 'SE') 

    raise Exception("Unable to determine input type")


def spreadsheet2dict(spreadFile):
    """
    Takes the filename of the spreadsheet, loads the data and organizes
    it into a dictionary"""

    spreadDict = {}
    key2field = {}
    for l, line in enumerate(open(spreadFile)):
        sl = line.strip().split('\t')
        if l == 0:
            for k, key in enumerate(sl):
                key2field[key] = k
        else:
            spreadDict[sl[key2field['fastq_analysis_id']]] = sl
    
    return (spreadDict, key2field)


def spreadsheet2RGdict(spreadFile, analysisID):
    """Compiles a read group dictionary from the information
    in the spreadFile for the given analysisID."""

    sD, k2f = spreadsheet2dict(spreadFile)

    try:
        rec = sD[analysisID]
    except KeyError:
        raise Exception('Information for analysis ID %s could not be found in %s' % (analysisID, spreadFile))

    ### build dictionary
    RG_dict = { 'ID' : '%s:%s' % (rec[k2f['center_name']], analysisID),
                'CN' : rec[k2f['center_name']],
                'LB' : 'RNA-Seq:%s:%s' % (rec[k2f['center_name']], rec[k2f['lib_id']]),
                'PL' : rec[k2f['platform']],
                'PM' : rec[k2f['platform_model']],
                'SM' : rec[k2f['sample_id']],
                'SI' : rec[k2f['submitter_sample_id']],
                'RG' : rec[k2f['read_group_label']].split(',')}

    return RG_dict


def xml2RGdict(xmlfile):
    
    ### read xml in
    root = etree.parse(xmlfile)
    rtree = root.find('Result')

    ### analysis_id
    analysis_id = rtree.find('analysis_id').text
    center = rtree.find('center_name').text
    try:
        date_string = rtree.find('analysis_xml/ANALYSIS_SET/ANALYSIS').attrib['analysis_date']
    except KeyError:
        date_string = rtree.find('run_xml/RUN_SET/RUN').attrib['run_date']
    sample_id = rtree.find('sample_id').text
    submitter_id = rtree.find('legacy_sample_id').text
    library_id = rtree.find('experiment_xml/EXPERIMENT_SET/EXPERIMENT').attrib['alias']
    platform = rtree.find('experiment_xml/EXPERIMENT_SET/EXPERIMENT/PLATFORM').getchildren()[0].tag 
    instrument = rtree.find('experiment_xml/EXPERIMENT_SET/EXPERIMENT/PLATFORM/*/INSTRUMENT_MODEL').text

    ### build dictionary
    RG_dict = { 'ID' : '%s:%s' % (center, analysis_id),
                'CN' : center,
                'DT' : date_string,
                'LB' : 'RNA-Seq:%s:%s' % (center, library_id),
                'PL' : platform,
                'PM' : instrument,
                'SM' : sample_id,
                'SI' : submitter_id}

    ### collect read group labels and add them to dict
    RG_dict['RG'] = []
    for x in rtree.find('analysis_xml/ANALYSIS_SET/ANALYSIS/ANALYSIS_TYPE/REFERENCE_ALIGNMENT/RUN_LABELS').getchildren():
        RG_dict['RG'].append(x.attrib['read_group_label'])
        
    return RG_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="ICGC RNA-Seq alignment wrapper for STAR alignments.", formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage='%(prog)s [options]', add_help=False)
    required = parser.add_argument_group("Required input parameters")
    required.add_argument("--genomeDir", default=None, help="Directory containing the reference genome index", required=True)
    required.add_argument("--tarFileIn", default=None, help="Input file containing the sequence information", required=True)
    optional = parser.add_argument_group("optional input parameters")
    optional.add_argument("--out", default="out.bam", help="Name of the output BAM file")
    optional.add_argument("--workDir", default="./", help="Work directory")
    optional.add_argument("--metaDataTab", default=None, help="File containing metadata for the alignment header")
    optional.add_argument("--analysisID", default=None, help="Analysis ID to be considered in the metadata file")
    optional.add_argument("--keepJunctions", default=False, action='store_true', help="keeps the junction file as {--out}.junctions")
    optional.add_argument("--useTMP", default=None, help="environment variable that is used as prefix for temprary data")
    optional.add_argument("--weakRGcheck", action='store_true', default=False, help="only perform weak RG record check and generate generic RG ID in case of a single alignment file with multiple RG records present. Use with caution!")
    optional.add_argument("-h", "--help", action='store_true', help="show this help message and exit")
    star = parser.add_argument_group("STAR input parameters")
    star.add_argument("--runThreadN", type=int, default=4, help="Number of threads")
    star.add_argument("--outFilterMultimapScoreRange", type=int, default=1, help="outFilterMultimapScoreRange")
    star.add_argument("--outFilterMultimapNmax", type=int, default=20, help="outFilterMultimapNmax")
    star.add_argument("--outFilterMismatchNmax", type=int, default=10, help="outFilterMismatchNmax")
    star.add_argument("--alignIntronMax", type=int, default=500000, help="alignIntronMax")
    star.add_argument("--alignMatesGapMax", type=int, default=1000000, help="alignMatesGapMax")
    star.add_argument("--sjdbScore", type=int, default=2, help="sjdbScore")
    star.add_argument("--alignSJDBoverhangMin", type=int, default=1, help="alignSJDBoverhangMin")
    star.add_argument("--genomeLoad", default="NoSharedMemory", help="genomeLoad")
    star.add_argument("--genomeFastaFiles", default=None, help="genome sequence in fasta format to rebuild index")
    star.add_argument("--outFilterMatchNminOverLread", type=float, default=0.33, help="outFilterMatchNminOverLread")
    star.add_argument("--outFilterScoreMinOverLread", type=float, default=0.33, help="outFilterScoreMinOverLread")
    star.add_argument("--twopass1readsN", type=int, default=-1, help="twopass1readsN (-1 means all reads are used for remapping)")
    star.add_argument("--sjdbOverhang", type=int, default=100, help="sjdbOverhang (only necessary for two-pass mode)")
    star.add_argument("--outSAMstrandField", default="intronMotif", help="outSAMstrandField")
    star.add_argument("--outSAMattributes", default=["NH", "HI", "NM", "MD", "AS", "XS"], help="outSAMattributes")
    star.add_argument("--outSAMunmapped", default="Within", help="outSAMunmapped")
    star.add_argument("--outSAMtype", default=["BAM", "SortedByCoordinate"], help="outSAMtype")
    star.add_argument("--outSAMheaderHD", default=["@HD", "VN:1.4"], help="outSAMheaderHD")
    star.add_argument("--outSAMattrRGline", default=None, help="RG attribute line submitted to outSAMattrRGline")
    star.add_argument("--outSAMattrRGfile", default=None, help="File containing the RG attribute line submitted to outSAMattrRGline")
    star.add_argument("--outSAMattrRGxml", default=None, help="XML-File in TCGA format to compile RG attribute line")
    
    ### check completeness of command line inputs
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    ### some sanity checks on command line parameters
    if args.metaDataTab is not None:
        if not os.path.exists(args.metaDataTab):
            raise Exception("File provided via --metaDataTab does not exist\nFile: %s" % args.metaDataTab)
        if args.analysisID is None:
            raise Exception("When providing information in a metadata file, a value for --analysisID is required")
    if args.outSAMattrRGxml is not None and not os.path.exists(args.outSAMattrRGxml):
        raise Exception("File provided via --outSAMattrRGxml does not exist\nFile: %s" % args.outSAMattrRGxml)
    if args.outSAMattrRGfile is not None and not os.path.exists(args.outSAMattrRGfile):
        raise Exception("File provided via --outSAMattrRGfile does not exist\nFile: %s" % args.outSAMattrRGfile)

    ### handling of input file (unpacking, etc. )
    if args.useTMP is not None:
        workdir = tempfile.mkdtemp(dir=os.environ[args.useTMP], prefix="star_inputdir_")
    else:
        workdir = tempfile.mkdtemp(dir=args.workDir, prefix="star_inputdir_")
    if args.tarFileIn.endswith(".gz"):
        tarcmd = "tar xvzf %s -C %s" % (args.tarFileIn, workdir)
    elif args.tarFileIn.endswith(".bz"):
        tarcmd = "tar xvjf %s -C %s" % (args.tarFileIn, workdir)
    elif args.tarFileIn.endswith(".tar"):
        tarcmd = "tar xvf %s -C %s" % (args.tarFileIn, workdir)
    else:
        raise Exception("Unknown input file extension for file %s" % (args.tarFileIn))
    subprocess.check_call(tarcmd, shell=True)
    
    ### collect fastq information from extraction dir
    align_sets = scan_workdir(os.path.abspath(workdir))
    
    if align_sets[2] == 'PE':
        read_str = '${fastq_left} ${fastq_right}'
    else:
        read_str = '${fastq_left}'
    
    ### simulate two pass alignment until STAR fully implements this
    if args.twopass1readsN != 0:
        
        ### run first round of alignments and only record junctions
        align_template_str_1st = """STAR \
--genomeDir ${genomeDir} --readFilesIn %s \
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
--sjdbOverhang ${sjdbOverhang} \
--outSAMstrandField ${outSAMstrandField} \
--outSAMtype None \
--outSAMmode None""" % read_str

        if args.twopass1readsN > 0:
            align_template_str_1st += " --readMapNumber %i" % args.twopass1readsN

        cmd = string.Template(align_template_str_1st).safe_substitute({
            'genomeDir' : os.path.abspath(args.genomeDir),
            'runThreadN' : args.runThreadN,
            'outFilterMultimapScoreRange' : args.outFilterMultimapScoreRange,
            'outFilterMultimapNmax' : args.outFilterMultimapNmax,
            'outFilterMismatchNmax' : args.outFilterMismatchNmax,
            'fastq_left' : ','.join([os.path.join(x[0], x[1]) for x in align_sets[1]]),
            'alignIntronMax' : args.alignIntronMax,
            'alignMatesGapMax': args.alignMatesGapMax,
            'sjdbScore': args.sjdbScore,
            'alignSJDBoverhangMin' : args.alignSJDBoverhangMin,
            'genomeLoad' : args.genomeLoad,
            'readFilesCommand' : align_sets[0],
            'outFilterMatchNminOverLread' : args.outFilterMatchNminOverLread,
            'outFilterScoreMinOverLread' : args.outFilterScoreMinOverLread,
            'sjdbOverhang' : args.sjdbOverhang,
            'outSAMstrandField' : args.outSAMstrandField
        })
        if align_sets[2] == 'PE':
            cmd = string.Template(cmd).substitute({
                'fastq_right' : ','.join([os.path.join(x[0], x[2]) for x in align_sets[1]])
            })

        ### take temp directory from environment variable
        if args.useTMP is not None:
            align_dir_1st = os.path.abspath( tempfile.mkdtemp(dir=os.environ[args.useTMP], prefix="star_aligndir_1st_") )
            genome_dir_1st = os.path.abspath( tempfile.mkdtemp(dir=os.environ[args.useTMP], prefix="star_genomedir_1st_") )
        else:
            align_dir_1st = os.path.abspath( tempfile.mkdtemp(dir=args.workDir, prefix="star_aligndir_1st_") )
            genome_dir_1st = os.path.abspath( tempfile.mkdtemp(dir=args.workDir, prefix="star_genomedir_1st_") )
        print "Running", cmd
        subprocess.check_call(cmd, shell=True, cwd=align_dir_1st)

        ### build index using provided genome fasta as well as junctions from first run
        cmd = """STAR --runMode genomeGenerate --genomeDir %s \
--genomeFastaFiles %s \
--sjdbOverhang %i \
--runThreadN %i \
--sjdbFileChrStartEnd %s""" % (genome_dir_1st, args.genomeFastaFiles, args.sjdbOverhang, args.runThreadN, os.path.join(align_dir_1st, 'SJ.out.tab')) 
        print "Running", cmd
        subprocess.check_call(cmd, shell=True, cwd=align_dir_1st)

        ### replace index for the second run with the one currently built
        genome_dir = genome_dir_1st
    else:
        genome_dir = os.path.abspath(args.genomeDir)


    align_template_str = """STAR \
--genomeDir ${genomeDir} --readFilesIn %s \
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
--sjdbOverhang ${sjdbOverhang} \
--outSAMstrandField ${outSAMstrandField} \
--outSAMattributes ${outSAMattributes} \
--outSAMunmapped ${outSAMunmapped} \
--outSAMtype ${outSAMtype} \
--outSAMheaderHD ${outSAMheaderHD}""" % read_str

#--twopass1readsN ${twopass1readsN} \

    cmd = string.Template(align_template_str).safe_substitute({
        'genomeDir' : genome_dir,
        'runThreadN' : args.runThreadN,
        'fastq_left' : ','.join([os.path.join(x[0], x[1]) for x in align_sets[1]]), #os.path.abspath(pair[1]),
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
        'sjdbOverhang' : args.sjdbOverhang,
        'outSAMstrandField' : args.outSAMstrandField,
        'outSAMattributes' : " ".join(args.outSAMattributes),
        'outSAMunmapped' : args.outSAMunmapped, 
        'outSAMtype' : " ".join(args.outSAMtype),
        'outSAMheaderHD' : " ".join(args.outSAMheaderHD)
    })
#        'twopass1readsN' : args.twopass1readsN,
    if align_sets[2] == 'PE':
        cmd = string.Template(cmd).substitute({
            'fastq_right' : ','.join([os.path.join(x[0], x[2]) for x in align_sets[1]]) # os.path.abspath(pair[2]),
        })

    ### process read group information
    if args.metaDataTab is not None:
        RG_dict = spreadsheet2RGdict(args.metaDataTab, args.analysisID) 
    elif args.outSAMattrRGxml is not None:
        RG_dict = xml2RGdict(args.outSAMattrRGxml)
    elif args.outSAMattrRGline is not None:
        RG_dict = dict([(x.split(':', 1)[0], x.split(':', 1)[1]) for x in args.outSAMattrRGline.split()])
    elif args.outSAMattrRGfile is not None:
        _fh = open(args.outSAMattrRGfile, 'r')
        RG_dict = dict([(x.split(':', 1)[0], x.split(':', 1)[1]) for x in _fh.next().strip().split()])
        _fh.close()
    else:
        RG_dict = {'ID' : '', 'SM' : ''}

    ### perform sanity check on provided RG records
    if 'RG' in RG_dict:
        if args.weakRGcheck and len(align_sets[1]) == 1 and len(RG_dict['RG']) > 1:
            RG_dict['RG'] = [RG_dict['ID']]
            print >> sys.stderr, 'WARNING: generated generic RG ID: %s' % RG_dict['ID']
        else:
            assert (len(align_sets[1]) == len(RG_dict['RG'])), 'Number of input file (pairs) does not match read groups in given RUN-xml' 

    ### post-process RG-dict to comply with STAR conventions
    for key in RG_dict:
        if key == 'RG':
            continue
        sl = RG_dict[key].split(' ')
        if len(sl) > 1:
            RG_dict[key] = '"%s"' % RG_dict[key]

    ### convert RG_dict into formatted RG line
    RG_line = []
    for r, readgroup in enumerate(align_sets[1]):
        if 'RG' in RG_dict:
            tmp = 'ID:%s:%s' % (RG_dict['ID'], RG_dict['RG'][r])
        else:
            tmp = 'ID:%s:%s' % (RG_dict['ID'], readgroup[0])
        if len(RG_dict) > 1:
            tmp += '\t'
            tmp += '\t'.join(['%s:%s' % (key, RG_dict[key]) for key in RG_dict if key not in ['ID', 'RG', 'SI']])
        ### add read group label
        if 'RG' in RG_dict and 'CN' in RG_dict:
            tmp += '\tPU:%s:%s' % (RG_dict['CN'], RG_dict['RG'][r])
        RG_line.append('%s' % tmp)
    cmd += ' --outSAMattrRGline %s' % ' , '.join(RG_line)

    ### handle comment lines
    comment_file = None
    if 'SI' in RG_dict:
        if args.useTMP is not None:
            comment_file = os.path.abspath( tempfile.mkstemp(dir=os.environ[args.useTMP], prefix="star_comments_")[1] )
        else:
            comment_file = os.path.abspath( tempfile.mkstemp(dir=args.workDir, prefix="star_comments_")[1] )
        
        fd_com = open(comment_file, 'w')
        fd_com.write('@CO\tsubmitter_sample_id:%s\n' % RG_dict['SI'])

        fd_com.flush()
        fd_com.close()
        
        cmd += ' --outSAMheaderCommentFile %s' % comment_file


    ### take temp directory from environment variable
    if args.useTMP is not None:
        align_dir = os.path.abspath( tempfile.mkdtemp(dir=os.environ[args.useTMP], prefix="star_aligndir_") )
    else:
        align_dir = os.path.abspath( tempfile.mkdtemp(dir=args.workDir, prefix="star_aligndir_") )
    print "Running", cmd
    subprocess.check_call(cmd, shell=True, cwd=align_dir)

    ### move output file
    if 'BAM' in args.outSAMtype and 'SortedByCoordinate' in args.outSAMtype:
        shutil.move(os.path.join(align_dir, 'Aligned.sortedByCoord.out.bam'), args.out)
    elif 'BAM' in args.outSAMtype and 'Unsorted' in args.outSAMtype:
        shutil.move(os.path.join(align_dir, 'Aligned.out.bam'), args.out)
    else:
        raise Exception('STAR output file could not be determined') 

    ### move junctions if to be kept
    if args.keepJunctions:
        shutil.move(os.path.join(align_dir, 'SJ.out.tab'), args.out + '.junctions')

    ### clean up working directory
    shutil.rmtree(workdir)
    shutil.rmtree(align_dir)
    if args.twopass1readsN != 0:
        shutil.rmtree(align_dir_1st)
        shutil.rmtree(genome_dir_1st)
