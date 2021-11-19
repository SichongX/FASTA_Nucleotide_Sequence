#6410 Assignment 2- Sichong Xu

#Import required modules.
from sys import argv 
import re
    
#Create functions needed for downstream analysis.
#Create function to sort nucleotie suquence in the right format(10pb per group, 60bp per line(66= 60+ 6 invisible characters))
def format(string, space = 10, group = 66):           
    string1= " ".join([string[i:i+space] for i in range(0, len(string), space)])
    return "\n".join([string1[j:j+group]for j in range(0,len(string1), group)])

#Create a function to join lines of FASTA file into one line(string)
def lines_to_line(list1):  
    a=''
    for line in list1:
        a += line.strip()
        b = a.split()
    list = ''.join(b)
    return(list)

#Reads in all required file, and create a output file for report.
#(The unix command line: python /Users/sx/Desktop/6410Programming/A2/sx_a2.py fasta_seq.fasta enzymes.txt out.txt)

sequence_file= open(argv[1], "r") #Reads in FASTA file with nucleotide sequence
enzyme_file= open(argv[2], "r") #Reads in text file with restriction enzymes
output_file= open(argv[3], "w") #Write the output to new file

output_file.write(f"Restriction enzyme analysis of nucleotide sequence from file {argv[1]}.\nCutting with enzymes from file {argv[2]}.")
output_file.write("\n" + '-' * 50 + "\n")


first_line = sequence_file.readline()  #Read the first line of FASTA file
if first_line.startswith(">"):      #If the first line contains the header information">"
    sequence_name = first_line.strip('>')   #Use as sequence name
    sequence1 = sequence_file.readlines() #Take the rest lines as nucleotide sequence
    sequence = lines_to_line(sequence1) #Use lines_to_line function to combine lines into one string line
    
else:
    sequence_name = "N/A"         #If not, use N/A as sequence name
    sequence1 = [first_line] + sequence_file.readlines() #Take the first and rest lines as whole sequence
    sequence= lines_to_line(sequence1)

output_file.write(f"Sequence name: {sequence_name}\n")
output_file.write(f"Sequence length: {len(sequence)} bp\n")



for line in enzyme_file.readlines(): #Read in each restriction enzyme, and their sequence, and their cutting site
    enzyme_name, enzyme_pattern = line.strip().split(";") 
    before_cutting_site, after_cutting_site = enzyme_pattern.split("%")

    cutting_sites = []   #Crearta a new list
    
    for seq in re.finditer(before_cutting_site + after_cutting_site, sequence):  #Find all hits that match the enzyme pattern in the sequence
        match_site = seq.start()      #Get the index of start matching position
        cutting_sites.append(match_site+len(before_cutting_site))    #Record all cuttung sites in to list


    if len(cutting_sites) != 0:  #If there are mactches for enzyme in sequence, print the count of cutting sites and sequence gragments
        output_file.write('-' * 50 + f"\nThere are {len(cutting_sites)} cutting sites for {enzyme_name}, cutting at {enzyme_pattern}\n")
        output_file.write(f"There are {len(cutting_sites)+1} fragments:\n\n")

        #Get length of each fragment and index info
        cutting_sites.append(len(sequence)) 
        start_cutting_site = 0
        for site in cutting_sites:
            fragment = sequence[start_cutting_site:site]
            output_file.write(f"Length: {len(fragment)} range: {start_cutting_site+1}-{site}\n")
            output_file.write(format(fragment) + "\n") #Use formating function to put fragment sequemce into the right format
            start_cutting_site = site


    if len(cutting_sites) == 0:   #If no cutting site found, print out result at the end
        no_cut = enzyme_name
output_file.write('-' * 50 + f"\nThere are no cutting sites for {no_cut} in this nucleotide sequence.\n")

  
sequence_file.close()  #Close all files
enzyme_file.close()
output_file.close()


