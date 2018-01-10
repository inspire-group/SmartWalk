from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.naive_bayes import MultinomialNB
# from sklearn.feature_selection import VarianceThreshold
# from sklearn.feature_selection import SelectKBest, SelectPercentile, f_regression
from sklearn import metrics, neighbors

import numpy
import getopt
import time
import sys
import math
import random
from math import sqrt
import matplotlib.pyplot as plt

def read_labels_txt(myfile):
    labels = numpy.loadtxt(myfile, dtype=int)
    # labels= numpy.reshape(labels,(numofemails,-1))
    return numpy.ravel(labels)

def main(argv):
  graph = ''
  classifier = ''
  start_time = time.time()

  # read path and classifier from command line
  try:
    opts, args = getopt.getopt(argv,'g:c:n:',['graph=','classifier=','num_nodes='])
  except getopt.GetoptError:
    print "python predictMixingTime.py -g <graph> -c <classifier> -n <num_nodes>"
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print "predictMixingTime.py -g <graph> -c <classifier>"
      sys.exit()
    elif opt in ("-g", "--graph"):
      graph = arg
    elif opt in ("-c", "--classifier"):
      classifier = arg
    elif opt in ("-n", "--num_nodes"):
      num_features = int(arg)

  # flag of feature selection
  fs=1
  if classifier == "svm":
    clf=SVR(C=1e3, epsilon=0.2)
  elif classifier == "NB":
    clf=MultinomialNB()
  elif classifier == "KNN":
    n_neighbors=10
    clf=neighbors.KNeighborsRegressor(n_neighbors, weights='uniform', p=1)
  elif classifier == "tree":
    clf=tree.DecisionTreeClassifier()
  elif classifier == "RF":
    clf=RandomForestRegressor(random_state=0, n_estimators=20)
  elif classifier == "OLS":
    clf = linear_model.LinearRegression()
  else:
    print ("classifier must be one of [svm, NB, KNN, tree, RF]\n")
    sys.exit(2)

  print 'Graph : ', graph
  print "Classifier : ", classifier

   # get the training and test examples, training and test labels
  rmse=[]
  r2=[]
  K=10
  
  # get label vector
  labels=read_labels_txt('../results/labels/length_'+graph+'_5000.txt')
  num_samples=numpy.shape(labels)[0]
  # generate feature matrix
  features=numpy.zeros((num_samples,num_features))
  out_file=open('../results/pred_length_'+graph+'_5000.txt', 'w')
  lines1=[line.rstrip('\n') for line in open('../results/features/neigh/neigh_'+graph+'_5000.txt')]
  lines2=[line.rstrip('\n') for line in open('../results/features/prob/prob_'+graph+'_5000.txt')]
  avg_num_neighs=0
  for j in range(num_samples):
    neighs=lines1[j].split(' ')
    probs=lines2[j].split(' ')
    for p in range(1,len(neighs)):
      if (neighs[p]!='' and probs[p]!=''):
        features[j, int(neighs[p])]=float(probs[p])
        avg_num_neighs+=1
  avg_num_neighs/=num_samples
  print("Avg number of 3-hop neighbours : %d  (%f)" % (avg_num_neighs, avg_num_neighs*1.0/num_features))

  for i in range(10):
    random.seed()
    # choose training examples
    train_ind=random.sample(range(num_samples), K)
    test_ind=list(set(range(num_samples))-set(train_ind))

    # # feature selection
    # # sel = VarianceThreshold(threshold=0.1)
    # sel = SelectPercentile(f_regression, percentile=5)
    # print sel.fit_transform(features[train_ind,:],labels[train_ind]).shape[1]
    # sel_ind = sel.get_support()

    # # fit and prefict
    # sel_features = features[:,sel_ind]
    sel_features = features
    clf.fit(sel_features[train_ind,:],labels[train_ind])
    results=clf.predict(sel_features[test_ind,:])
    pred_labels=[]
    
    for result in results:
      pred_labels.append(int(round(result)))
      if i==0:
        out_file.write(str(int(round(result))))
        out_file.write('\n')
    test_labels=labels[test_ind]

    r2.append(metrics.r2_score(test_labels, pred_labels))
    rmse.append(sqrt(metrics.mean_squared_error(test_labels, pred_labels)))
    print("RMSE = %f4.7" % rmse[len(rmse)-1])
    print("r2 = %f4.7" % r2[len(r2)-1])

  print("Avg. RMSE = %f4.7" % numpy.array(rmse).mean())
  print("Avg. r2 = %f4.7" % numpy.array(r2).mean())

  print str(time.time() - start_time)


if __name__ == "__main__":
   main(sys.argv[1:])


# for p in range(10):
#   clf=RandomForestRegressor(random_state=0, n_estimators=20)
#   clf.fit(features[sort_ind[p+500,1:50],:],labels[sort_ind[p+500,1:50]])
#   result=clf.predict(features[p+500,:])
#   print('pred len = %d' % int(round(result)))
#   print('original len = %d', test_labels[p])
#   print('prev. pred len = %d', pred_labels[p])
#   new_pred_labels.append(result)



