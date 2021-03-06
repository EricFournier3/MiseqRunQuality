# coding=utf-8

from xml.dom.minidom import parse
import xml.dom.minidom
from decimal import *
import os
import subprocess
from subprocess import Popen, PIPE
import re
import logging
import argparse
import sys
import stat
import shutil
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary
import pandas as pd
import yaml
"""
Eric Fournier
2019-03-26

Programme permettant de calculer les metrics d'une run MiSeq

Example de commande sur slbio
/usr/bin/python MiSeqStat6.py --help
/usr/bin/python MiSeqStat6.py  --runno TESTINTEROP --bact salmonella --param path_to_yaml_file

Liste de modifications
- Modif_20190514: Eric Fournier 2019-05-14 => changer le repertoire FASTQ dans lequel se trouve les sequence brutes sur LSPQ_MiSeq. 2_SequencesBrutes au lieu de FASTQ 

- Modif_20190515: Eric Fournier 2019-05-15 => changer le repertoire dans lequel se trouve InterOp, RunInfo.xml, run Parameters.xml et le fichier de resultats pour /mnt/Partage/LSPQ_MiSeq/RunName/3_Analyse

-Modif_20190625: Eric FOurnier 2019-06-25 => ajouter l option --param pour lire le fichier de parametre contenant les path
"""



logging.basicConfig(level=logging.INFO)


##################################  Global Var #################################

#Parsing de la ligne de commande
parser = argparse.ArgumentParser(description="Calculateur des statistique de runs MiSeq")
parser.add_argument("-r","--runno",help="Nom de la run dans S/Partage/LSPQ_MiSeq",required=True)

#Taille du genome
parser.add_argument("-g","--gl",help="Taille de genome",required=True)

parser.add_argument("-p","--param",help="path vers le fichier de parametre",required=True)

args_commandline = parser.parse_args(sys.argv[1:])
args = args_commandline.__dict__
project_name =  args["runno"]
genome_length = int(args["gl"])
path_param_file = args["param"]

#exit(0)

snakemake_param_handle = open(path_param_file)
all_dict = yaml.load(snakemake_param_handle)
snakemake_param_handle.close()

#Modif_20190515
#Modif_20190625
'''
lspq_miseq_experimental_dir = "1_Experimental"
lspq_miseq_sequencebrute_dir = "2_SequencesBrutes"
lspq_miseq_analyse_dir = "3_Analyse"
'''
lspq_miseq_experimental_dir = all_dict["lspq_miseq_subdir"][0]
lspq_miseq_sequencebrute_dir = all_dict["lspq_miseq_subdir"][1]
lspq_miseq_analyse_dir = all_dict["lspq_miseq_subdir"][2]



#Repertoire local temporaire pour les calculs
#***********************************************************   INSPQ-9788 ***********************************************************
#temp_dir = "/home/ericf/TEMP_FASTQ"
#***********************************************************   slbio ***********************************************************
temp_dir = "/home/foueri01@inspq.qc.ca/temp/TEMP_FASTQ"

#Script R qui calcule les metrics
#***********************************************************    INSPQ-9788 ***********************************************************
#rScript = "/home/ericf/ProjetProgrammation/ProjetR/ComputeReadsStat2.R"
#***********************************************************    slbio ***********************************************************
rScript = "/home/foueri01@inspq.qc.ca/GitScript/MiSeqRunQuality/ComputeReadsStat2.R"


if not os.path.isdir(temp_dir):
    os.system("mkdir {0}".format(temp_dir))
else:
    os.system("rm -rf {0}".format(temp_dir))
    os.system("mkdir {0}".format(temp_dir))

#Repertoire de la run
#Modif_20190625
#basedir = os.path.join("/mnt/Partage/LSPQ_MiSeq/",project_name)
basedir = os.path.join(all_dict["path"][0],project_name)

#Quelques check-up
if not os.path.isdir(basedir):
    logging.error(basedir + " est inexistant")
    exit(0)

#Modif_20190514
#Repertoire contenant les fastq
#fastq_dir = os.path.join(basedir,"FASTQ")
fastq_dir = os.path.join(basedir,lspq_miseq_sequencebrute_dir)
if not os.path.isdir(fastq_dir):
    logging.error(fastq_dir + " est inexistant")
    exit(0)

#Modif_20190515
interop_dir = os.path.join(basedir,lspq_miseq_analyse_dir,"InterOp")
if not os.path.join(basedir,interop_dir):
    logging.error(interop_dir + " est inexistant")
    exit(0)

#On s assure qu il y a des fastq
if not os.listdir(fastq_dir):
    logging.error("Aucun fastq dans " + fastq_dir)
    exit(0)

#Le fichier  RunInfo.xml
#Modif_20190515
runinfo_file = os.path.join(basedir,lspq_miseq_analyse_dir,"RunInfo.xml")
print "run info ", runinfo_file
if not os.path.isfile(runinfo_file):
    logging.error("Le fichier {0} est absent".format(runinfo_file))
    exit(0)

#Le fichier runParameters.xml
#Modif_20190515
runparam_file = os.path.join(basedir,lspq_miseq_analyse_dir,"runParameters.xml")
if not os.path.isfile(runinfo_file):
    logging.error("Le fichier {0} est absent".format(runparam_file))
    exit(0)

#Fichier de resultats final
outfile = open(os.path.join(temp_dir,"MiSeqStat_" + project_name + "TEMP.txt"),'w')
outfile_append = open(os.path.join(temp_dir,"MiSeqStat_" + project_name + ".txt"),'a+')

#Metric de chacun des specimens
allfile_qc_dict = {}

#Key = nom du specimen  Value = [nombre total de nucleaotid, genome coverage]
allspec_cov_dict = {}


##################################  End Global Var #################################

##################################  Begin Function #################################

def ComputeGenomeCoverage(specname,nBnucleotid):
    """
    Calcul de la couverture du genome pour ce specimen
    :param specname:
    :param nBnucleotid:
    """
    #print specname, " ",nBnucleotid, " ",allspec_cov_dict[specname][0]
    nBnucleotid_r1_r2 = float(nBnucleotid) + float(allspec_cov_dict[specname][0])

    cov = round(nBnucleotid_r1_r2 / genome_length, 0)

    allspec_cov_dict[specname].append(cov)


##################################  End Function #################################


##################################  Begin Program #################################
logging.info("              Start Calculation")

logging.info("              Calcul des metrics de la run")

#Recuperation du Q30 pour la run
run_metrics = py_interop_run_metrics.run_metrics()
#Modif_20190515
run_folder = run_metrics.read(os.path.join(basedir,lspq_miseq_analyse_dir))
summary = py_interop_summary.run_summary()
py_interop_summary.summarize_run_metrics(run_metrics, summary)
summary.total_summary().yield_g()

columns = (('% Over Q30', 'percent_gt_q30'),)
rows = [('Total', summary.total_summary()),]

d = []
for label, func in columns:
    d.append((pd.Series([getattr(r[1], func)() for r in rows], index=[r[0] for r in rows])))

parse_d = re.search(r'Total\s{4}(\S+)',str(d[0]))
percent_gt_q30 = round(float(parse_d.group(1)),0)
#print percent_gt_q30

#Recuperation des valeurs de cluster pour la run
def format_value(val):
    if hasattr(val, 'mean'):
        return val.mean()
    else:
        return val

read = 0
columns = (('Density (K/mm2)', 'density'),('% Cluster PF','percent_pf'))
rows = [summary.at(read).at(lane) for lane in xrange(summary.lane_count())]
d2 = []
for label, func in columns:
    d2.append( (pd.Series([format_value(getattr(r, func)()) for r in rows])))

density = str(d2[0])
percent_pf = str(d2[1])

parse_density = re.search(r'\S+\s{4}(\S+)',density)
density = parse_density.group(1)
density = round(float(density) / 1000,0)

parse_percent_pf = re.search(r'\S+\s{4}(\S+)',percent_pf)
percent_pf = parse_percent_pf.group(1)
percent_pf = round(float(percent_pf),0)

#print density
#print percent_pf

#Calculs des Q30 pour les samples avec R
logging.info("              Calcul des metrics des samples dans R")
os.system("Rscript {0} {1} {2} ".format(rScript,fastq_dir,temp_dir))

#Contruction du dictionnaire de metrics
try:
    # Fichier de resultat metrics genere par le rScript
    metric_file_from_R = open(os.path.join(temp_dir,"fastqStat.txt"))

    for line in metric_file_from_R:

        #print "line is ", line

        line_parse = re.search(r'(\S+)\t(\S+)\t(\S+)\t(\S+)',line)
        fastq_file = line_parse.group(1)

        if(fastq_file.find("RUN") != -1):
            spec_name = "RUN"
        else:
            spec_name = fastq_file[:-3]

        min_q30_perc = line_parse.group(2)
        nb_read = line_parse.group(3)
        nb_nucleotid = line_parse.group(4)
        allfile_qc_dict[fastq_file] = [min_q30_perc,nb_read]

        if spec_name in allspec_cov_dict.keys():
            pass
            ComputeGenomeCoverage(spec_name,nb_nucleotid)
        elif spec_name.find("RUN") == -1:
            allspec_cov_dict[spec_name] = [nb_nucleotid]

    metric_file_from_R.close()

except:
    logging.error("Probleme de lecture du fichier fastqStat.txt")
    exit(0)

#Transfert des metrics dans le fichier de resultats final
logging.info("              Lecture des metrics")

outfile.write("ID\tNb_Reads\tCluster_Density_K_mm2\tCluster_Passing_Filter\tOver_Q30\n")

for fastqfile in allfile_qc_dict.keys():
    if(fastqfile != "RUN"):
        #print fastqfile
        min_q30_perc = allfile_qc_dict[fastqfile][0]
        nb_read = allfile_qc_dict[fastqfile][1]
        outfile.write("{0}{1}{2}{1}{3}{1}{4}{1}{5}\n".format(fastqfile, "\t", nb_read, 'NA', 'NA',min_q30_perc))
    else:
        min_q30_perc = allfile_qc_dict[fastqfile][0]
        nb_read = allfile_qc_dict[fastqfile][1]
        outfile.write("{0}{1}{2}{1}{3}{1}{4}{1}{5}\n".format(fastqfile, "\t", nb_read, density, percent_pf, percent_gt_q30))

outfile.close()

sortfile = os.path.join(temp_dir,"MiSeqStat_" + project_name + ".txt")
os.system("awk 'NR<2{print $0;next}{print $0 | \"sort -k1\" }' " + outfile.name + "> " + sortfile)

#On ajoute les valeurs de couverture
outfile_append.write("\nID\tCoverage\n")
for my_spec_name in allspec_cov_dict.keys():
    outfile_append.write("{0}\t{1}\n".format(my_spec_name,allspec_cov_dict[my_spec_name][1]))

outfile_append.close()

#print allspec_cov_dict
#Modif_20190515
os.system("sudo cp {0} {1}".format(sortfile,os.path.join(basedir,lspq_miseq_analyse_dir)))

logging.info("              End Calculation")

exit(0)

