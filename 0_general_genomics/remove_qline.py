import sys
import random


sam_file = open(sys.argv[1], 'r')
sam_file_out = open('temp_qline_removed','w')
flags = {}

i=0
for line in sam_file:
        num = random.randrange(100,200,10)
        if line[0] =="+":
            sam_file.next()        
        elif line[0] == "@":
            sam_file_out.write( '>L1S'+ str(num)+ '_'+ str(i) +' '+ line[1:] )
            i= i+1
        else:   
            sam_file_out.write( line )
sam_file.close()
sam_file_out.close()
