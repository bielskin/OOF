import re
import os
import argparse
from multiprocessing import Pool, cpu_count

def get_origin_seqs(input_seq_file):
    origin_seqs=[]
    Infile=open(input_seq_file)
    for Line in Infile:
        Line = Line.strip()
        if Line[0]=='>':
            Name=Line[1:]  
            Name=Name.split()
            Name=Name[0]
            origin_seqs.append(Name)
        else:  
            continue
    Infile.close()
    return origin_seqs

def run_meme(input_seq_file, nmotifs):
    meme_command = "meme " + input_seq_file + " -nmotifs " + str(nmotifs) + " -dna"
    os.system(meme_command)

def get_genome_list():
    genome_list=[]
    for filename in os.listdir():
        if filename.startswith("."):
            continue
        else:
            genome_list.append(filename)
    return genome_list


def run_mast(genome_list, working_dir):
    with Pool(processes=cpu_count()) as pool:
  	    pool.starmap(os.system, [('mast ' + working_dir + "/meme_out/meme.txt " + working_dir + 
                        '/Genomes/'+ genome  + ' -o ' + working_dir + '/Genomes/'+ 'mastres_' + genome,) for genome in genome_list])

    


def get_hits(ref_folder_key, origin_seqs, working_dir):
    for mastresults in os.listdir():
        if mastresults.startswith("mastres_"):
            #Thresholding with original taxa
            #Retrieve back search results and determine the threshold
            if mastresults == ref_folder_key:
                os.chdir(ref_folder_key)
                file = open('mast.txt')
                pattern=re.compile('SEQUENCE NAME                      DESCRIPTION.                  E-VALUE  LENGTH')
                endpatt=re.compile('SECTION II: MOTIF DIAGRAMS')
                line_number=0
                mylines=[]
                for line in file:
                    line_number += 1
                    line=line.strip()
                    result=re.search(pattern, line)
                    end=re.search(endpatt, line)
                    mylines.append(line)
                    if result:
                        firstline=line_number+1
                    elif end:
                        lastline=line_number
                hits=mylines[firstline:lastline-7]    
                HitDict={}
                hitkeys=[]
                for line in hits:
                    line=str.split(line)
                    HitDict[line[0]]=line[(len(line)-2)]
                    hitkeys.append(line[0])
                seqcount=0
                dictcounter=0
                for key in hitkeys:
                    dictcounter=dictcounter+1
                    if key in origin_seqs:
                        seqcount=seqcount+1
                        if seqcount==len(origin_seqs):
                            seqid_next_it=dictcounter
                    else:
                        continue
                seqid_next=(hitkeys[seqid_next_it])   
                ethreshold=HitDict[seqid_next]
                chandir_command=working_dir+'/Genomes'
                os.chdir(chandir_command)
    return ethreshold

def get_sig_hits(ref_folder_key, ethreshold, working_dir):
    '''
    This function gets the significant hits and writes an error file if the HitDict is being overwritten
    '''
    sig_hits_master={}
    errur = None
    for mastresults in os.listdir():
        if mastresults.startswith("mastres_"):       
            if mastresults != ref_folder_key:
                current_folder=mastresults
                os.chdir(current_folder)
                file=open('mast.txt')
                pattern=re.compile('SEQUENCE NAME                      DESCRIPTION.                  E-VALUE  LENGTH')
                endpatt=re.compile('SECTION II: MOTIF DIAGRAMS')
                line_number=0
                mylines=[]
                sig_hits=[]


                for line in file:
                    line_number += 1
                    line=line.strip()
                    result=re.search(pattern, line)
                    end=re.search(endpatt, line)
                    mylines.append(line)
                    if result:
                        firstline=line_number+1
                    elif end:
                        lastline=line_number

                hits=mylines[firstline:lastline-7]       

                HitDict={}
                hitkeys=[]
                for line in hits:
                    line=str.split(line)
                    for key in HitDict: 
                        if key == line[0]: 
                            errur = [mastresults, key]
                    HitDict[line[0]]=line[(len(line)-2)]
                    hitkeys.append(line[0])
                for key in hitkeys:
                    if float(HitDict[key]) < float(ethreshold):
                        sig_hits.append(key)

                chandir_command=working_dir+'/Genomes'
                os.chdir(chandir_command)
                #Add all significant hits to dictionary
                sig_hits_master[mastresults]=sig_hits
                
                
    os.chdir(working_dir)
    if errur:
        Errorfile = 'Error_log.txt'
        WriteOutFile = True
        Errorfile = open(Errorfile, 'w')
        Errorfile.write('Your HitDict for the ' + str(errur[0]) + ' file was overwritten due to matching geneIDs')
        
    return sig_hits_master
    
    


def get_sig_hits_seqs(sig_hits_master):
    sig_hits_seqs={}
    for key in sig_hits_master:
        genome_file_name=key.lstrip("mastres_")
        Infile = open(genome_file_name,'r')
        RecordNum = -1 
        Sequences=[]
        SeqDict={}
        for Line in Infile:
            Line = Line.strip()
            if Line[0]=='>':
                Name=Line[1:]  
                Name=Name.split()
                Name=Name[0]
                Sequences.append([Name,''])
                RecordNum += 1
                SeqKey = Name
                SeqDict[SeqKey] = '' 
            else:  
                if RecordNum > -1:  
                    Sequences[RecordNum][1] += Line
                    SeqDict[SeqKey] += Line

        for hit in sig_hits_master[key]:
            actualseq=SeqDict[hit]
            sig_hits_seqs[hit]=actualseq
        
    return sig_hits_seqs



def write_files(sig_hits_master, sig_hits_seqs, ethreshold, working_dir):

    os.chdir(working_dir)

    #Write out significant hits summary
    Outfilesumname='Summary_of_ORFs.txt'
    Outfilesum=open(Outfilesumname, 'w')
    Outfilesum.write("Significant hits found in genomes searched, Genome:Number of hits" + '\n'+ '\n')
    Outfilesum.write("Evalue threshold determined:")
    Outfilesum.write(str(ethreshold))
    Outfilesum.write('\n'+ '\n')
    for key in sig_hits_master:
        numberofhits=len(sig_hits_master[key])
        Outfilesum.write(key.lstrip("mastres_") + ":")
        Outfilesum.write(str(numberofhits))
        Outfilesum.write('\n')


    #Write the significant hits in fasta format
    OutFileName='Results.fasta'
    OutFile=open(OutFileName, 'w')
    for hit in sig_hits_seqs:
        OutFile.write(">" + hit + '\n')
        OutFile.write(sig_hits_seqs[hit] +'\n')
        


def main():
    
    parser = argparse.ArgumentParser(description='Orthologous ORF Finder (OOF): Determines motifs for a set of gene family sequences and searches provided genomes for orthologs')
    parser.add_argument("-e", "--evalue_threshold",  help="Specify the e-value threshold for MAST hits")
    parser.add_argument("-n", "--nmotifs",  help="Specify the number of motifs to search for in MEME", default=10)
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-i", "--input_seq_file",  help="File of input gene family sequences", required=True)
    requiredNamed.add_argument("-o", "--origin_taxa", help="Name of original taxon genome annotation", required=True)
    args = parser.parse_args()

    input_seq_file=(vars(args)['input_seq_file'])
    origin_taxon=(vars(args)['origin_taxa'])
    nmotifs=(vars(args)['nmotifs'])
    evaluethreshold_spec=(vars(args)['evalue_threshold'])

    
    working_dir = os.getcwd()
    
    origin_seqs = get_origin_seqs(input_seq_file)
    run_meme(input_seq_file, nmotifs)


    chandir_command = working_dir +'/Genomes'
    os.chdir(chandir_command)

    genome_list = get_genome_list()

    run_mast(genome_list, working_dir)
        
    ref_folder_key= 'mastres_' + origin_taxon

    if evaluethreshold_spec==None:
        ethreshold = get_hits(ref_folder_key, origin_seqs, working_dir)
   
    else:
        ethreshold=evaluethreshold_spec



    sig_hits_master = get_sig_hits(ref_folder_key, ethreshold, working_dir)

    os.chdir(chandir_command)
    sig_hits_seqs = get_sig_hits_seqs(sig_hits_master) # returns sig_hit_seqs

    write_files(sig_hits_master, sig_hits_seqs, ethreshold, working_dir)
    

if __name__ == '__main__':
    main()

