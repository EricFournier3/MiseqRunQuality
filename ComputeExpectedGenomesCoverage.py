import sys
import os
import subprocess
from decimal import *
"""
Eric Fournier 2019-07-12
Calculs des couvertures de genome estimees selon la taille de genome de l organisme

"""

run_path, project_path, ref_file, lspq_miseq_run_path, genome_cov_file, slbio_fastqbrut_dir = sys.argv[1:8]

'''
run_path est le path slbio de la run
project_path est le path slbio du project
ref_file est le path slbio vers le fichier contenant les tailles des genomes de reference
lspq_miseq_run_path est le path LSPQ_MISEQ de la run
genome_cov_file est le fichier de sortie avec les couvetures
slbio_fastqbrut_dir est le sous repertoire 1_FASTQ_BRUT du projet

'''


def ComputeCoverage(count):
    """
    Calcul de la couverture selon le nombre total de nucleotides et la taille du genome
    :param count:
    :return:
    """
    try:
        cov = count[0]/count[1]
        return int(round(cov,0))
    except:
        return 'NA'

#le fichier de couverture en sortie
genome_cov_handle = open(genome_cov_file,'a')

#le fichier contenant les tailles de genomes de reference
ref_file_handle = open(ref_file)
ref_file_handle.close()

#la sample sheet de la run
sample_sheet_handle = open(os.path.join(project_path,os.path.basename(run_path[:-1]) + ".csv.temp3"))
sample_sheet_handle.readline()

#le fichier contenant les nombres de reads
read_count_file = os.path.join(project_path,slbio_fastqbrut_dir,'ReadCount.txt')


for line in sample_sheet_handle:
    line = line.split(',')[:-1]
    sample, organism = line[0],line[11]

    #calcul du nombre total de nucleotides
    cmd_1 = "sed -n /{0}/p {1}".format(sample,read_count_file) + " | awk 'NR==1{print $2*2*300}'"
    total_read_length = subprocess.check_output(cmd_1, shell=True)

    #taille du genome de reference de cet organisme
    try:
        cmd_2 = "sed -n '/^\"{0}[ \"]/Ip' ".format(organism) + ref_file + " | awk 'BEGIN{FS=\"\t\"}NR==1{print $2*1000000}'"
        genome_length = subprocess.check_output(cmd_2, shell=True)
        int_val = int(genome_length)
    except:
        genome_length = 0

    #la couverture estimee
    cov = ComputeCoverage(map(Decimal,[total_read_length,genome_length]))

    #on enregsitre la couverture dans le fichier de sortie
    genome_cov_handle.write('{0}\t{1}\t{2}\t{3}\n'.format(sample,organism,str(genome_length).strip(),str(cov)))
sample_sheet_handle.close()

genome_cov_handle.close()





