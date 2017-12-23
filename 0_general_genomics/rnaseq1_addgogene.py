import sys

# python ../genomics/match_go_add_goinfo.py 
############profile31.txt gene_names H37Rv.prot_table H37Rv_goids 


if1= open(sys.argv[1], 'r')
if2 = open(sys.argv[2],'r')
if3 = open(sys.argv[3],'r')
if4 = open(sys.argv[4],'r')
of1 = open('temp_1_gene_prot_go_added','w')

rvlist =[]

for line in if1:
    rv = line.strip().split("\t")
    rv1 = rv[0]
    rvlist.append(rv1.lower())

genename=[]
# to get gene names 

for line2 in if2:
    words = line2.strip().split("\t")
    w1 = words[0]
    if  w1.lower() in rvlist :
    	genename.append(line2)

# to get info from prot table

protname =[]
# to get gene names from prot table
for line2 in if3:
    words = line2.strip().split("\t")
    for w in words:
    	if  w.lower() in rvlist :
        	protname.append( w.lower() + ' ' + words[0] +'_' +words[1])


goname = []

for line2 in if4:
    words = line2.strip().split()
    for w in words:
    	if  w.lower() in rvlist :
        	goname.append( line2)




'''
for i in rvlist:
	print i

for i in genename:
	print i

for i in protname:
	print i

for i in goname:
	print i 

'''





i = 0
j = 0
k = 0
l=0
while i < len(rvlist):
	rvcode = rvlist[i]
	try:
		words = genename[j].strip().split()
		generv = words[0].lower()
		#print generv
		if rvcode == generv :
			gene = genename[j].strip()
			j = j + 1
		else:
			gene = 'null'
	except IndexError:
		gene = 'null'

	try:
		words = protname[k].strip().split()
		protrv = words[0].lower()
		#print protrv
		if rvcode == protrv :
			protein = protname[k].strip()
			k = k + 1
		else:
			protein = 'null'
	except IndexError:
		protein = 'null'
    
        try:
		words = goname[l].strip().split()
		gorv = words[0].lower()
		#print generv
		if rvcode == gorv :
			go = goname[l].strip()
			l = l + 1
		else:
			go = 'null'
	except IndexError:
		go = 'null'
        of1.write( gene + '\t' + protein[7:] + '\t' + go[7:] + '\n')
        i = i + 1

#'''rvcode + '\n' + gene + '\n' +'''


if1.close()
if2.close()
if3.close()
if4.close()
of1.close()
