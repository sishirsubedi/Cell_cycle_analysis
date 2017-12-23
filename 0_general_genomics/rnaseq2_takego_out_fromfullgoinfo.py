import sys


# input temp file from match_go.py, temp, goinfo

if1= open(sys.argv[1], 'r')

of1 = open('temp2_goinfo_longerversion','w')


goinfo ={}
for line in if1:
    words = line.strip().split()
    if words[0] == "id:":
        goinfo[words[1]] = "test"
        temp = if1.next()
        goinfo[words[1]] = temp[5:]



for items in goinfo:
    of1.write( items + '\t' + goinfo[items] + '\n')


            
if1.close()
of1.close()
