import sys
import random


sam_file = open(sys.argv[1], 'r')
sam_file_out = open('diversity','w')
flags = {}

diversity ={}

for line in sam_file:
        if line[0] =="#":
            continue 
        line.strip()
        word = line.split("\t") 
        otu_id = word[1]+ word[2]         

        print word[0], otu_id
           
for otu in diversity:
	sam_file_out.write(otu.ljust(40) + str(diversity[otu])+'\n')

sam_file.close()
sam_file_out.close()
