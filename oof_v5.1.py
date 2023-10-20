import re
import os
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt

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
    print("Running MEME motif search")
    meme_command = "meme "  + input_seq_file + " -nmotifs " + str(nmotifs) + " -dna" + " -nostatus " + "-oc meme_out"
    os.system(meme_command)

def get_genome_list():
    genome_list=[]
    for filename in os.listdir():
        if filename.startswith("."):
            continue
        if filename.startswith("mastres_"):
            continue
        else:
            genome_list.append(filename)
    return genome_list

def sort_origin_seqs(genome_list, origin_seqs, genome_folder):
    os.chdir(genome_folder)
    origin_dict={}
    for genome in genome_list:
        genome_ids=[]
        Infile=open(genome)
        for Line in Infile:
            Line = Line.strip()
            if Line[0]=='>':
                Name=Line[1:]  
                Name=Name.split()
                Name=Name[0]
                genome_ids.append(Name)
            else:  
                continue
        genome_spec_ids=[]
        for original_sequence in origin_seqs:
            if original_sequence in genome_ids:
                genome_spec_ids.append(original_sequence)
        if genome_spec_ids==[]:
            continue
        else:
            origin_dict[genome]=genome_spec_ids
        Infile.close()
    return origin_dict

def run_mast(genome_list, working_dir, genome_folder, iterations):
	print("Running MAST searches for hits")
	with Pool(processes=cpu_count()) as pool:
            pool.starmap(os.system, [('mast ' + working_dir + "/meme_out/meme.txt " + genome_folder + '/' + genome  + ' -nostatus ' + '-oc ' + genome_folder + '/'+ 'mastres_'+'iteration_'+ str(iterations-1) + '_' + genome,) for genome in genome_list])
		
def ethresh_determine(working_dir, genome_folder):
    print("Determining E-Value Threshold")
    os.chdir(genome_folder)
    ethreshold_dict={}
    for mastresults in os.listdir():
        if mastresults.startswith("mastres_"):
            current_folder=mastresults
            os.chdir(current_folder)
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
            evals=[]
            ratios = []
            ratioratios = []
            subratios = []
            
            for line in hits:
                line=str.split(line)
                HitDict[line[0]]=line[(len(line)-2)]
                hitkeys.append(line[0])
                evals.append(line[(len(line)-2)])
            i=0
            ratios=[]
            for j in evals:
                if float(j)>0:
                    ratios.append(float(j)/float(evals[i-1]))
                i=i+1
            ratioratios=[]
            i2=0  
            for r in ratios:
                if float(r)>0:
                    ratioratios.append(r/ratios[i2-1])
            subratios=[]
            i3=-1
            for x in ratioratios:
                if i3>=0:
                    subratios.append(abs(x-ratioratios[i3]))
                i3=i3+1
            maximum_value = float(10000000000000000)
            maximum_index = 0
            for index, value in enumerate(subratios):
                if value >= maximum_value:
                    maximum_value = value
                    maximum_index = index

            os.chdir(genome_folder)
            if maximum_index==0:
                maximum_index=1
            
            ethreshold_dict[current_folder]=(evals[maximum_index-1])
            
    os.chdir(genome_folder)
    return ethreshold_dict

def get_sig_hits(ethreshold_dict, working_dir, genome_folder, iterations):
    '''
    This function gets the significant hits and writes an error file if the HitDict is being overwritten
    '''
    os.chdir(genome_folder)
    
    print("Thresholding hits")
    os.chdir(genome_folder)
    sig_hits_master={}
    errur = None
    for mastresults in os.listdir():
        if mastresults.startswith("mastres_iteration_" + str(iterations-1)):    
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
            ethreshold=ethreshold_dict[current_folder]
            for key in hitkeys:
                if float(ethreshold)<0.1:
                    if float(HitDict[key]) <= float(ethreshold):
                        sig_hits.append(key)
   
            os.chdir(genome_folder)
            #Add all significant hits to dictionary
            sig_hits_master[mastresults]=sig_hits
                
                
    os.chdir(working_dir)
    if errur:
        Errorfile = 'Error_log.txt'
        WriteOutFile = True
        Errorfile = open(Errorfile, 'w')
        Errorfile.write('Your HitDict for the ' + str(errur[0]) + ' file was overwritten due to matching geneIDs')
        
    return sig_hits_master
    
def get_sig_hits_seqs(sig_hits_master, genome_folder, iterations, iterative):
    sig_hits_seqs={}
    os.chdir(genome_folder)
    for key in sig_hits_master:
        if iterative==True:
            strip="mastres_iteration_"+ str(iterations-2) +'_'
        elif iterative==False:
            strip="mastres_iteration_"+ str(iterations-1) +'_'
        genome_file_name=key.replace(strip, "")
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

def write_files(sig_hits_master, sig_hits_seqs, ethreshold_dict, working_dir, iterations):

    os.chdir(working_dir)
    os.makedirs("OOF_Output", exist_ok=True)
	#Write out significant hits summary
    Outfilesumname=working_dir + '/OOF_Output/Summary_of_ORFs_Iteration_' + str(iterations-1) + '_.txt'
    with open(Outfilesumname, 'w') as Outfilesum:
        Outfilesum.write("Significant hits found in genomes searched" + '\n'+ '\n')
        Outfilesum.write("Total Number of Iterations:" +  str(iterations-1) + '\n'+ '\n')
        for key in sig_hits_master:
            numberofhits=len(sig_hits_master[key])
            Outfilesum.write(key.lstrip("mastres_") + ":" + '\n')
            Outfilesum.write('      ' + "E-value Threshold:" + ' '+ str(ethreshold_dict[key]) + '\n')
            Outfilesum.write('      ' + "Number of homologs detected:" + str(numberofhits) + '\n' + '\n')
            Outfilesum.write('\n')


    #Write the significant hits in fasta format
    OutFileName=working_dir + '/OOF_Output/Results_' + str(iterations-1) + '_.fasta'
    with open(OutFileName, 'w') as OutFile:
        for hit in sig_hits_seqs:
            OutFile.write(">" + hit + '\n')
            OutFile.write(sig_hits_seqs[hit] +'\n')

def create_barplots_from_summary_files(working_dir):
    # Get the list of summary files
    folder_path = working_dir + "/OOF_Output"
    file_list = sorted([file for file in os.listdir(folder_path) if file.startswith('Summary_of_ORFs_Iteration_')])

    # Extract the genome names from the summary files
    genome_names = set()
    for file in file_list:
        with open(os.path.join(folder_path, file), 'r') as f:
            contents = f.read()
            matches = re.findall(r"iteration_\d+_(.*?)\.cds\.fa", contents)
            genome_names.update(matches)

    # Determine the number of rows and columns for subplots
    num_genomes = len(genome_names)
    num_cols = 3 if num_genomes > 3 else num_genomes
    num_rows = (num_genomes - 1) // num_cols + 1

    # Create a figure and axes for the bar plots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 9))

    # Iterate over the genome names and create bar plots for each genome
    for i, genome in enumerate(genome_names):
        # Calculate the row and column index in the matrix
        row = i // num_cols
        col = i % num_cols

        # Create lists to store the iteration numbers and number of homologs
        iterations = []
        num_homologs = []

        # Iterate over the summary files and extract the data for the current genome
        for file in file_list:
            with open(os.path.join(folder_path, file), 'r') as f:
                contents = f.read()
                match = re.search(fr"iteration_\d+_{re.escape(genome)}\.cds\.fa:\s+E-value Threshold:.*?Number of homologs detected:(\d+)", contents, re.DOTALL)
                if match:
                    iteration = re.search(r"Iteration_(\d+)", file).group(1)
                    iterations.append(int(iteration))
                    num_homologs.append(int(match.group(1)))

        # Create a bar plot for the current genome
        if isinstance(axes, np.ndarray):  # If axes is a numpy array
            ax = axes[row, col]
        else:
            ax = axes  # If axes is not a numpy array (only one plot or 3 or fewer genomes)
        ax.bar(iterations, num_homologs)
        ax.set_xlabel('Iteration')
        ax.set_xticks(iterations) 
        ax.set_ylabel('Number of Homologs')
        ax.set_title(f'Genome: {genome}')

    # If there is only one genome or 3 or fewer genomes, remove any unused subplots
    if num_genomes <= 3:
        for i in range(num_genomes, num_rows * num_cols):
            row = i // num_cols
            col = i % num_cols
            fig.delaxes(axes[row, col])
        # If the number of subplots is not divisible by 3, remove the extra empty subplots
    if num_genomes % 3 != 0:
        num_empty_plots = 3 - (num_genomes % 3)
        for i in range(num_empty_plots):
            row = (num_genomes + i) // num_cols
            col = (num_genomes + i) % num_cols
            fig.delaxes(axes[row, col])
    # Adjust the spacing between subplots
    fig.tight_layout()

    # Save the figure as a PDF file
    output_pdf = working_dir + "/OOF_Output/genome_plots.pdf"
    fig.savefig(output_pdf, format='pdf')


def str_to_bool(s):
    if s.lower() in ('true', 't', 'yes', 'y', '1'):
        return True
    elif s.lower() in ('false', 'f', 'no', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected (true/false).')

def main():
    print("OOF: Orthologous ORF Finder")
    parser = argparse.ArgumentParser(description='Orthologous ORF Finder (OOF): Determines motifs for a set of gene family sequences and searches provided genomes for orthologs')
    parser.add_argument("-e", "--evalue_threshold",  help="Specify the e-value threshold for MAST hits")
    parser.add_argument("-n", "--nmotifs",  help="Specify the number of motifs to search for in MEME", default=10)
    parser.add_argument("-i", "--iterative",  type=str_to_bool, help="Specify if you want OOF to iteratively search", default=False)
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-s", "--input_seq_file",  help="File of input gene family sequences", required=True)
    requiredNamed.add_argument("-g", "--genome_folder_path", help="Path to folder of CDS genome annotations", required=True, type=str)
    args = parser.parse_args()
    
    input_seq_file=(vars(args)['input_seq_file'])
    genome_folder=os.path.abspath(args.genome_folder_path)
    nmotifs=(vars(args)['nmotifs'])
    evaluethreshold_spec=(vars(args)['evalue_threshold'])
    working_dir = os.getcwd()
    origin_seqs = get_origin_seqs(input_seq_file)
    iterative=args.iterative
    if args.iterative==True:
        iterations=1
        iteration_loop=True
        while iteration_loop:
            iteration_loop=False
            if iterations==1:
                run_meme(input_seq_file, nmotifs)
                iteration_loop=True
            else:
                new_seqs="OOF_Output/Results_" + str(iterations-1)+ "_.fasta"
                run_meme(new_seqs, nmotifs)
    
            os.chdir(genome_folder)
    
            genome_list = get_genome_list()
    
            #origin_dict = sort_origin_seqs(genome_list, origin_seqs, genome_folder)

            run_mast(genome_list, working_dir, genome_folder, iterations)
  

            if evaluethreshold_spec==None:
                ethreshold_dict = ethresh_determine(working_dir, genome_folder)
            else:
                ethreshold=evaluethreshold_spec
    
            sig_hits_master = get_sig_hits(ethreshold_dict, working_dir, genome_folder, iterations)
       
            if iterations > 1: 
                for key in sig_hits_master:
                    search = r"iteration_\d+_"
                    replace="iteration_" + str(iterations-2)+ "_"
                    previous_key = re.sub(search, replace, key)
                    if set(list(sig_hits_master[key]))!= set(list(previous_sig_hits_master[previous_key])):
                        iteration_loop=True
                        break
            elif iterations==11:
                break        

    
            previous_sig_hits_master = sig_hits_master
            iterations += 1


            sig_hits_seqs = get_sig_hits_seqs(sig_hits_master, genome_folder, iterations, iterative) # returns sig_hit_seqs
            write_files(sig_hits_master, sig_hits_seqs, ethreshold_dict, working_dir, iterations)
    
        create_barplots_from_summary_files(working_dir)
    else:
        iterations=2
        run_meme(input_seq_file, nmotifs)   
        os.chdir(genome_folder)
        genome_list = get_genome_list()
        run_mast(genome_list, working_dir, genome_folder, iterations)
        if evaluethreshold_spec==None:
            ethreshold_dict = ethresh_determine(working_dir, genome_folder)
        else:
            ethreshold_dict = {}
            for mastresults in os.listdir():
                if mastresults.startswith("mastres_iteration_"):
                    ethreshold_dict[mastresults] = evaluethreshold_spec
        sig_hits_master = get_sig_hits(ethreshold_dict, working_dir, genome_folder, iterations)
        sig_hits_seqs = get_sig_hits_seqs(sig_hits_master, genome_folder, iterations, iterative) # returns sig_hit_seqs
        write_files(sig_hits_master, sig_hits_seqs, ethreshold_dict, working_dir, iterations)
        create_barplots_from_summary_files(working_dir)

if __name__ == '__main__':
    main()
