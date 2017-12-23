import sys
import random


infile1 = open(sys.argv[1], 'r')
infile2 = open(sys.argv[2], 'r')
outfile1 = open('clusterwithGO','w')


goinfo = {}
for line in infile1:
        words=line.strip().split("\t")
        rv = words[0].lower()
	rv=rv[0:6]
        goinfo[rv]=str(words[2:])       

for i in goinfo:
	print i,goinfo[i]
           
for line in infile2:
        words=line.strip().split("\t")
        rv = words[0].lower()
	rv=rv[0:6]
        outfile1.write(words[0] +"\t"+goinfo[rv]+'\n') 

infile1.close()
infile2.close()
outfile1.close()
