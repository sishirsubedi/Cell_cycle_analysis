import sys
import numpy as np
import matplotlib.pyplot as plt


in_file = open(sys.argv[1], 'r')

sample = 'x'
species =[]
composition =[]
for line in in_file:
	words = line.split("\t") 
	if "#Sample" in words:
                count = 1
		sample = words[1]
		newline = in_file.next()
                while  len(newline.strip()) != 0 :
                        words = newline.split("\t")
			species.append(words[0])
			composition.append(int(words[1]))
                        newline = in_file.next()
	
	        plt.figure()
                y_pos = np.arange(len(species))
       		plt.barh(y_pos, composition, align='center', alpha=0.4)
       		plt.yticks(y_pos, species)
                filename = sample + str(count)
                pdf.savefig(filename)
                plt.close()
                plt.clf()
                count = count +1
                for i in species:
			del i
                for i in composition:
                        del i
                

          
in_file.close()
