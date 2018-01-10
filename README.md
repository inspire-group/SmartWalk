# SmartWalk
This repository contains code for SmartWalk: SmartWalk: Enhancing Social Network Security via Adaptive Random Walks. SmartWalk is a security enhancing system which incorporates adaptive random walks in social network security applications. We utilize a set of supervised machine learning techniques to predict the necessary random walk length based on the structural characteristics of a social graph. 

## Usage
1. computeMixingTime.cpp

This provides the functionality of computing mixing time from an input graph. 
```
$ ./computeMixingTime.out graph num_nodes
```

2. getFeatures.cpp

This offers the functionality of obtaining features (probability distributions of k-hop neighbours) from an input graph. 
```
$ ./getFeatures.outgraph num_nodes num_chosen_nodes num_hops
```

3. predictMixingTime.py

This trains a classifier based on training samples and predict the local mixing time of a graph.
```
$ python predictMixingTime.py -g <graph> -c <classifier> -n <num_nodes>
```
