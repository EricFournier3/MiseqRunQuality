from glob import glob
import subprocess
import sys, os
import logging
import argparse

logging.basicConfig(level=logging.INFO)

bowtie2_db_path = '/home/foueri01@inspq.qc.ca/InternetProgram/Bowtie2/DATABASE'

doc = """    

    usage : /usr/bin/python Qualimap2.py -o <outdir> -r <reference> -f <fastqdir> -m <bowtie2|bwa|snalt> -c <y|n> 

    Exemple : /usr/bin/python Qualimap2.py -o /home/foueri01@inspq.qc.ca/temp/TEMP2/SMALT_MAP -r  /home/foueri01@inspq.qc.ca/temp/TEMP2/GCF_000007645.1_ASM764v1_genomic.fna -f /home/foueri01@inspq.qc.ca/temp/TEMP2 -m smalt -c n
         
 """

parser = argparse.ArgumentParser(description="Compute mapping quality")
parser.add_argument("-o","--outdir",help="Output directory",required=True)
parser.add_argument("-r","--ref",help="Path vers la sequence de reference",required=True)
parser.add_argument("-f","--fastqdir",help="Path vers fastq directory",required=True)
parser.add_argument("-m","--mapper",help="Mapping program",nargs=1,type=str,choices=['bowtie2','bwa','smalt'],required=True)
parser.add_argument("-c","--concat",help="Fastq concatenation",nargs=1,type=str,choices=['y','n'],required=True)

args_commandline = parser.parse_args(sys.argv[1:])

args = args_commandline.__dict__
outdir =  args["outdir"]
path_ref = args["ref"]
path_fastq = args["fastqdir"]
mapper = args["mapper"][0]
concat_bool = args["concat"][0]

if not os.path.isdir(outdir):
    os.mkdir(outdir)

#fastq_list = [x for x in os.listdir(path_fastq) if x.endswith('fastq.gz')]
#print fastq_list

fastq_R1_list = sorted(glob(path_fastq + '/*_R1_*.fastq.gz'))
fastq_R2_list = sorted(glob(path_fastq + '/*_R2_*.fastq.gz'))

#print fastq_R1_list
#print fastq_R2_list

cat_fastq_R1 = os.path.join('TEMP2','TEMP_R1.fastq.gz')
cat_fastq_R2 = os.path.join('TEMP2','TEMP_R2.fastq.gz')
out_sam = os.path.join('TEMP2','map.sam')

fastq_paired_list = []

#repertoire temporaire
if not os.path.isdir('TEMP2'):
    os.mkdir('TEMP2')

#concatener les fastq
if concat_bool == 'y':
    logging.info("Concat reads")
    os.system("cat {0} > {1}".format(' '.join(fastq_R1_list), cat_fastq_R1))
    os.system("cat {0} > {1}".format(' '.join(fastq_R2_list), cat_fastq_R2))
else:
    fastq_paired_list = zip(fastq_R1_list,fastq_R2_list)

#print fastq_paired_list

def BowtieMap(*args):

    logging.info("Map with bowtie2")

    fq_r1 = None
    fq_r2 = None

    if len(args):
        fq_r1 = args[0]
        fq_r2 = args[1]
    else:
        fq_r1 = cat_fastq_R1
        fq_r2 = cat_fastq_R2

    logging.info("Map reads")
    os.system('bowtie2 -x {0} -1 {1} -2 {2} -S {3} --no-unal --no-discordant'.format('TEMP_DB',fq_r1, fq_r2, out_sam))

def SmaltMap(*args):

    logging.info("Map with smalt")

    fq_r1 = None
    fq_r2 = None

    if len(args):
        fq_r1 = args[0]
        fq_r2 = args[1]
    else:
        fq_r1 = cat_fastq_R1
        fq_r2 = cat_fastq_R2

    os.system('smalt map -l pe -i 500 -o {0} {1} {2} {3}'.format(out_sam,'ref_index',fq_r1,fq_r2))

def BwaMap(*args):

    logging.info("Map with Bwa")

    fq_r1 = None
    fq_r2 = None

    if len(args):
        fq_r1 = args[0]
        fq_r2 = args[1]
    else:
        fq_r1 = cat_fastq_R1
        fq_r2 = cat_fastq_R2

    os.system('bwa mem -t 20 {0} {1} {2} > {3}'.format(path_ref,fq_r1,fq_r2,out_sam))

def Samtools():
    os.system('{0} view -h -bS {1}.sam >  {1}.bam ; {0} sort {1}.bam {1}_sort ; {0} index {1}_sort.bam'.format('samtools',os.path.join('TEMP2','map')))

def Qualimap(*args):

    mapprog = args[0]

    out = os.path.join(outdir,mapprog) if len(args) == 1 else os.path.join(outdir,mapprog + '_' + args[1])

    print " ********************************************** OUT IS ", out

    os.system('{0} bamqc -bam {1} -outdir {2} -nt 10 --java-mem-size=4G -outformat PDF:HTML -c'.format('qualimap', os.path.join('TEMP2', 'map_sort.bam'),out))

def Clean():
    os.system('rm -r TEMP2/*')


if mapper == 'bowtie2':

    logging.info("Index reference")
    # indexer la reference. Les index sont places dans /home/foueri01@inspq.qc.ca/InternetProgram/Bowtie2/DATABASE car bowtie2 cible les db via la variable d environnement BOWTIE2_INDEXES
    os.system('bowtie2-build -q -f {0} {1}'.format(path_ref, 'TEMP_DB'))
    os.system('mv *.bt2 {0}'.format(bowtie2_db_path))


    if concat_bool == 'y':
        BowtieMap()
        Samtools()
        Qualimap('Bowtie')
        Clean()
    else:
        for r1,r2 in fastq_paired_list:
            BowtieMap(r1,r2)
            specname = str(os.path.basename(r1)).split('_')[0]
            Samtools()
            Qualimap('Bowtie',specname)
            Clean()

elif mapper == 'bwa':

    if concat_bool == 'y':
        BwaMap()
        Samtools()
        Qualimap('Bwa')
        Clean()
    else:
        for r1,r2 in fastq_paired_list:
            BwaMap(r1,r2)
            specname = str(os.path.basename(r1)).split('_')[0]
            Samtools()
            Qualimap('Bwa',specname)
            Clean()

elif mapper == 'smalt':

    os.system('smalt index -k 20 -s 13 ref_index {0}'.format(path_ref))

    if concat_bool == 'y':
        SmaltMap()
        Samtools()
        Qualimap('Smalt')
        Clean()
    else:
        for r1,r2 in fastq_paired_list:
            SmaltMap(r1,r2)
            specname = str(os.path.basename(r1)).split('_')[0]
            Samtools()
            Qualimap('Smalt',specname)
            Clean()
    try:
        os.remove('ref_index.sma')
        os.remove('ref_index.smi')
    except:
        pass
else:
    exit('Abort > no mapper')


#On conserve seulement les rapports pdf
for spec_map_dir in glob(outdir + "/*"):
    if os.path.isdir(spec_map_dir):
        #print spec_map_dir
        prefix = ''

        print

        if mapper == 'bowtie2':
            prefix = 'Bowtie_'
        elif mapper == 'bwa':
            prefix = 'Bwa'
        else:
            prefix = 'Smalt'

        report_name = os.path.basename(spec_map_dir) + '.pdf'
        #print report_name

        specname = os.path.basename(spec_map_dir).lstrip(prefix)
        #print specname

        map_dir = os.path.dirname(spec_map_dir)
        #print map_dir

        for k in glob(spec_map_dir + "/report.pdf"):
            #print k
            os.system('cp {0} {1}'.format(k,os.path.join(map_dir,report_name )))

        os.system('rm -r {0}'.format(spec_map_dir))




exit("Finish")


