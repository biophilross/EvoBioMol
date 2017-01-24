#!python

###############################################################
# NWV 18 Jan 2013                                              
# This script will read in a multiple sequence alignment  
# return a file with the number of alternative amino acids at
# each site. This outfile is used in conjunction with ModifyBfactor.py
# to color amino acids by how diverse they are.                     
###############################################################

import os, re
from Bio import AlignIO


# "chdir" = "change directory" Tell Python to change where it
# will look for files. For example, you'll probably be storing
# your evol-biol-mol files in a directory like
# C:/Users/Nicholas/Documents/EvolBiolMol/HW
# Change this to be your working directory.

base_directory = "/Users/tylernstarr/Dropbox/evol-biol-mol/homeworks/"

os.chdir(os.path.expanduser(base_directory))

# Define a variable for each of your input and output files.
# If you have your directory structured correctly, you should not
# have to change these lines.

alignment_filename = "hw2/data/SR214-plus-1GWRa.fasta"
results_filename = "hw2/results/b-factorList.txt"

# Open the results file so we can write ('w') to it. 
results_file = file(os.path.join(base_directory,results_filename),'w')

# Use Biopython to read in the entire alignment. 
alignment = AlignIO.read(alignment_filename, "fasta")

# Print a summary of the alignment to the screen.
print(alignment)

# Store the length of the alignment into a file
aligned_seq_length = len(alignment[1,:])

print "Length of aligned sequences (including gaps): {}".format(aligned_seq_length)

# Store the number of sequences in your alignment and print it to your
# screen.

num_seqs = len(alignment)
print "Number of sequences in the alignment: {}".format(len(alignment))


# Initialize an array with all zeros. See the Python tutorial on data structure
# to better understand what an array (vs. a dictionary for example) is.

histogram = [0]*(aligned_seq_length+1)
print histogram
print len(histogram)

# For each site...
for site_number in range(aligned_seq_length):
                 
    # pull out the site data, and give it a convenient name
    site_aas = alignment[:,site_number]

    # if less than 20% of the sequences have gaps...
    fraction_gaps = site_aas.count('-')/float(num_seqs)
    if fraction_gaps < 0.2:
        
        # put the amino acids in a set 

        number_of_unique_aas = len(set(site_aas))

        # ...and print this number for the site that we're looking at
        histogram[site_number+1] = number_of_unique_aas
    else:
        histogram[site_number+1] = 0
for ind in range(len(histogram)):
    results_file.write("{site}\t{num}\n".format(site=(ind+313), num=histogram[ind]))
# We're done with this file, so close it.
results_file.close()




                         
