import sys

# input temp file from match_go.py, go_info

if1= open(sys.argv[1], 'r')
if2 = open(sys.argv[2],'r')

of1 = open('temp3_goinfo_added','w')
of2 = open ('temp4_go_composition','w')

goinfo ={}
for line in if2:
    words = line.strip().split("\t")
    gonum = words[0]
    godef = words[1] 
    goinfo[gonum] = godef 

'''
for i in goinfo:
	print i , goinfo[i]
'''
goinfo_this ={}

for line in if1:
    words = line.strip().split("\t")
    templist = []
    for w in words:
    	templist.append(w)
    of1.write('\n\n'+ line.strip())
    for gos in templist:
    	if gos in goinfo:
    		if goinfo[gos] in goinfo_this:
    			goinfo_this[goinfo[gos]] += 1
    		else:
    			goinfo_this[goinfo[gos]] = 1
    		of1.write( '\n' + '@' + gos + '\t' + goinfo[gos])
           
for i in goinfo_this:
	of2.write( str(i) + '\t' + str(goinfo_this[i]) + '\n')

if1.close()
if2.close()
of1.close()
of2.close()


