

import numpy,scipy.optimize
import sys,math
import numpy as np

def stats(vals):
  sum,ss = 0,0
  for x in vals: sum += x; ss += x*x
  N = float(len(vals))
  mean = sum/N
  var = ss/N-mean*mean
  if var<0: var = 0
  stdev = math.sqrt(var)
  return mean,stdev





data = []
skip = 1
for line in open("total_deseq_norm_skip0hr.txt"):
  if skip>0: skip -= 1; continue
  w = line.split()
  data.append(w)

PI = 3.1415927
TP = [3., 6.5, 9., 12., 18.5, 21., 27., 31., 33., 36., 39.5, 42., 45.5, 52., 55.]
N = len(data)

#def F(x,ampl,freq,phas,const):
#  return ampl*numpy.sin(2.*PI*freq*(x/55.)+phas)+const
def F(x,ampl,freq,phas,trend,const):
  return ampl*numpy.sin(2.*PI*freq*(x/55.)+phas)+x*trend+const

for i,w in enumerate(data):
  strain = w[0]
  #X = numpy.array(TP+TP)
  #Y = [float(x) for x in w[1:31]]
  X = numpy.array(TP)

  ## cos
  Y1 = [float(w[i]) for i in range(1,16)]
  Y2 = [float(w[i]) for i in range(16,31)]

  # ## rv
  #Y1 = [float(w[i]) for i in range(31,46)]
  #Y2 = [float(w[i]) for i in range(46,61)]


  #m = sum(Y1+Y2)/len(Y1+Y2)
  #Y1 = [y/float(m) for y in Y1]
  #Y2 = [y/float(m) for y in Y2]
  m,s = stats(Y1+Y2)
  Y1 = [(y-m)/s for y in Y1]
  Y2 = [(y-m)/s for y in Y2]
  Y = [(y1+y2)/2. for y1,y2 in zip(Y1,Y2)]
  delta = [0.25*(y1-y2)**2 for y1,y2 in zip(Y1,Y2)]
  #delta = [abs (y1 - y2) for y1, y2 in zip(Y1, Y2)]
  ss = sum(delta)-max(delta)
  #absdelta = [abs (y1 - y2) for y1, y2 in zip(Y1, Y2)]
  kurt = [(y1 - m) ** 4 for y1,m in zip(Y1,Y)]
  krt = sum(kurt)-max(kurt)

  #try: params,covar = scipy.optimize.curve_fit(F,X,Y,bounds=([0.,1.0,-PI,-2.],[5.,2.0,PI,2.]))

  try: params,covar = scipy.optimize.curve_fit(F,X,Y,bounds=([0.,1.0,-PI,-0.02,-2.],[5.,2.0,PI,0.02,2.]))
  except: sys.stderr.write("fitting failed for %s\n" % strain); print strain,"fitting-failed"; continue


  stderrs = numpy.sqrt(numpy.diag(covar))

  Yhat = F(X,*params)
  RSS = sum([(y-yhat)**2 for y,yhat in zip(Y,Yhat)])



  vals = [strain,"%0.1f" % m,"|"]
  for y in Y1: vals.append("%0.3f" % y)
  vals.append("|")
  for y in Y2: vals.append("%0.3f" % y)
  vals.append("|")
  for y in Yhat: vals.append("%0.3f" % y)
  vals.append("|")
  for j in range(len(params)):
    vals.append("%0.3f" % params[j])
    vals.append("%0.3f" % stderrs[j])
  vals.append("%0.3f" % RSS)
  vals.append("%0.3f" % ss)
  vals.append("%0.3f" % krt)

  print '\t'.join(vals)

