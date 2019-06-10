#inspire de https://gist.github.com/rchikhi/7281991

from glob import glob
import subprocess
import sys, os
import logging

logging.basicConfig(level=logging.INFO)

doc = """
Quickly estimates insert sizes of read datasets, given some sequence(s) they can be mapped to.
Author: Eric Fournier  2019-05-15    

short usage: <reference> <output file> <*.fastq.gz> 
         
         Exemple de commande:
/usr/bin/python ComputeInsertSize.py  /home/foueri01@inspq.qc.ca/temp/TEMP2/GCF_000007645.1_ASM764v1_genomic.fna /home/foueri01@inspq.qc.ca/temp/TEMP2/InsertSizeValues.txt /home/foueri01@inspq.qc.ca/temp/TEMP2/*.fastq.gz 
 """

if len(sys.argv) < 4:
    exit(doc)

try:
    subprocess.call(["bwa"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except:
    exit("Please make sure that the `bwa` binary is in your $PATH")

def Unzip(fastq_gz):
    logging.info("Unzipping fastq files " + os.path.basename(fastq_gz))
    os.system("gunzip {0}".format(fastq_gz))
    return str(fastq_gz).strip('.gz')

nb_threads = 5

#la reference
reference = sys.argv[1]

#fichier de resultats
outfile = sys.argv[2]
outfile_handle = open(outfile,'w')

#les fichiers reads
reads = sorted(sys.argv[3:])

#Unzip les fastq
reads = map(Unzip,reads)

#On verifie si tous les fichier y sont
for read in reads:
    if not os.path.isfile(read):
        exit("Error: %s does not exist" % read)

#ON index la reference
if not os.path.isfile(reference+".sa"):
    logging.info("Creating index file..")
    subprocess.call(["bwa", "index", reference], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def parse_list(line, nb_elts):
    #print "LINE IS ", line
    # specific to BWA-MEM stderr format
    return map(lambda x: int(float(x)), ' '.join(line.strip().replace(',','').split()[-nb_elts:])[1:-1].split())

stats = dict()

#Calcul de la taille moyenne des inserts pour chacun des paired end read
for read1, read2 in zip(reads[::2],reads[1::2]):
    #print read1, read2
    logging.info("Processing: %s with %s " % (os.path.basename(read1), os.path.basename(read2)))

    outfile_handle.write("Processing: %s with %s " % (os.path.basename(read1), os.path.basename(read2)) + "\n")

    cmd = ["bwa", "mem"] + (["-p"] if read2 == "" else []) + ["-t %d" % nb_threads, reference, read1, read2]
    print cmd
    break
    DEVNULL = open(os.devnull, 'wb')
    process = subprocess.Popen(cmd, stdout=DEVNULL, stderr=subprocess.PIPE)
    seen_candidate_line = False
    while True:
        line = process.stderr.readline()
        #print "line is ", line
        if line == '' and process.poll() != None:
            break
        if "worker" in line:
            break
        if "pestat" not in line:
            continue
        if "candidate unique pairs for" in line:
            if seen_candidate_line:
                break
            seen_candidate_line = True
            nb_pairs = parse_list(line,4)
            for i in xrange(4):
                stats[['FF', 'FR', 'RF', 'RR'][i]] = { 'nb_pairs' : nb_pairs[i] }
        if "orientation" in line:
            orientation = line.strip().split()[-1].replace('.','')
        if "mem_pestat] mean and std.dev:" in line:
            mean, stdev = parse_list(line,2)
            stats[orientation]['mean'] = mean
            stats[orientation]['stdev'] = stdev
            if orientation == 'RR':
                # stats are complete
                break
        #sys.stdout.write(line)
        #sys.stdout.flush()
    if process.poll() is None:
        process.terminate()

    results = sorted(stats.items(), key=lambda x: x[1]['nb_pairs'], reverse=True)
    most_likely = results[0]
    mean = most_likely[1]['mean']
    stdev = most_likely[1]['stdev']
    logging.info("Orientation " + most_likely[0] + " mean: " + str(mean) + " stdev: " + str(stdev))

    outfile_handle.write("Orientation " + most_likely[0] + " mean: " + str(mean) + " stdev: " + str(stdev) + '\n\n')

outfile_handle.close()



