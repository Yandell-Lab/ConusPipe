/usr/bin/python -w
from __future__ import division
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.semi_supervised import LabelSpreading
from sklearn.linear_model import Perceptron
import argparse


    
# Load your classification data here
X = []
Xtest=[]
y = []
parser=argparse.ArgumentParser()
parser.add_argument("sampleTestFeatureFile", type=str)
parser.add_argument("sampleTrainingFeatureFile", type=str)
parser.add_argument("sampleTrainingOutcomeFile", type=str)
parser.add_argument("logitOutputFile", type=str)
parser.add_argument("semiSupOutputFile", type=str)
parser.add_argument("neuroNetOutputFile", type=str)
args = parser.parse_args()

fx=open(args.sampleTrainingFeatureFile,'r')
for line in fx:
        #X.append(line.strip().split(','))
        a=[]
        a.append(line.strip().split('\t'))
        i=0;
        while i<16:
            X.append(float(a[0][i]))
            i=i+1
       
        #X.append(tuple(float(a[0][0]),float(a[0][1])??not work

X=np.asarray(X)
XLength=len(X)

X=X.reshape(XLength//16,16)
fxt=open(args.sampleTestFeatureFile,'r')
for line in fxt:
        #X.append(line.strip().split(','))
        a=[]
        a.append(line.strip().split('\t'))
        i=0;
        while i<16:
            Xtest.append(float(a[0][i]))
            i=i+1
       
        #X.append(tuple(float(a[0][0]),float(a[0][1])??not work

Xtest=np.asarray(Xtest)
XtestLength=len(Xtest)

Xtest=Xtest.reshape(XtestLength//16,16)
                        
                        

fy=open(args.sampleTrainingOutcomeFile,'r')
for line in fy:
	y.append(float(line))
 
y=np.asarray(y)

#logistic regression
log_reg = LogisticRegression(C=10**10)
log_reg.fit(X, y)
y_hat = log_reg.predict(Xtest)
out = open(args.logitOutputFile,'w')
out.write(','.join(map(str,y_hat)))
out.close()
#semi-supervised learning LabelSpreading
#choose  K = 5 alpha= 0.8 
clf = LabelSpreading(kernel='knn',n_neighbors=3,alpha=1,n_jobs=15)
Xin=np.concatenate((X,Xtest),axis=0)
lenXt=len(Xtest)
labels=np.empty(lenXt)
labels.fill(-1)
yin=np.concatenate((y,labels),axis=0)
clf.fit(Xin,yin)
y_hat = clf.predict(Xtest)
out = open(args.semiSupOutputFile,'w')
out.write(','.join(map(str,y_hat)))
out.close()

#neuroNet work perceptron
#choose  n_iter = 60 
clf = Perceptron(n_iter=5,random_state=0,n_jobs=20)
clf.fit(X,y)
y_hat = clf.predict(Xtest)
out = open(args.neuroNetOutputFile,'w')
out.write(','.join(map(str,y_hat)))
out.close()
