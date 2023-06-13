# Orthologous ORF Finder (OOF)
OOF is a program created to utilize existing bioinformatic and sequence analysis tools to easily identify homologs of a gene family across a large variety of species

## Installation
Installing OOF requires downloading the python program and its dependencies in a Linux or Mac operating system. To run OOF you will need the python script, python3 (or later versions), and the MEME suite available at (https://meme-suite.org/meme/)

## Usage
python oof_v3.3.py - i FASTA formatted file of known gene family sequences -o genome file of species that the sequences in input are from (-n) specify number of motifs for MEME search, default: 10 (-e) pre-determined e-value threshold, if known (-h) help menu

## Example Files
The provided Example Files folder contains the file set-up for a regular OOF search. All genomes to be searched must be CDS genome files in fasta format contained in a "Genomes" directory folder. The example files genome folder contains six sample CDS genomes from within Angiosperms. Within the Example Files folder where is also a fasta format file of 4 gene sequences of the LORELEI and LORELI-like gene family (LLG) from Arabidopsis thaliana (LLG_seqs.fasta)
To run an example OOF search, download the example files folder and in a linux based terminal navigate to the folder with the cd command. Run OOF using the run command: 'python oof_v3.3.py -i LLG_seqs.fasta -o Athaliana_167_TAIR10.cds.fa'

## Output
OOF will output two files for review. The first is the Summary of Results file which is a plain text file that reports the determined e-value threshold from the fun and a list of the genomes searched with the number of hits found. The second file is the Results.fasta file which is a fasta formatted file containing all the hits found from all genomes searched.

## Contributing
Inform Nick of any contributions

## License
Copyright (c) 2021 Nicholas Bielski

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (OOF), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
