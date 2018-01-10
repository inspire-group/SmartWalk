#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <assert.h>
#include <boost/timer.hpp>

using namespace std;

/* Global variables */ 
// The following should be the user input
string graph_file, node_file, feature_file;
int num_nodes, num_chosen_nodes, num_hops;		// number of initial nodes		

double *degree;	// node degree
set<int> *graph;	// maintains topology
int *nodelist;

void generate_social_topology(); // generates real topology
void get_k_hop_neighbours(int node, int k, set<int> &neigh); 

int main(int argc, char * argv[])
{
	if(argc==5){
		num_nodes=atol(argv[2]);
        printf("num_nodes = %d\n", num_nodes);
        num_chosen_nodes=atol(argv[3]);
        printf("num_chosen_nodes = %d\n", num_chosen_nodes);
        num_hops=atol(argv[4]);
        printf("num_hops = %d\n", num_hops);
        stringstream ss, ss2;
        ss << num_chosen_nodes;
        ss2 << num_hops;
        graph_file="../datasets/ungraph_"+string(argv[1])+".txt";
        node_file="../datasets/nodelist/nodelist_"+string(argv[1])+".txt";
        feature_file="../results/features/feature_"+string(argv[1])+"_"+ss.str()+"_"+ss2.str()+".txt";

    }
    else{
        printf("Usage: $ ./%s graph num_nodes num_chosen_nodes num_hops\n", argv[0]);
        exit(1);
    }

    boost::timer t;
	
	// allocate memory
    degree=new double[num_nodes];
    graph=new set<int>[num_nodes];
    nodelist=new int[num_chosen_nodes];
    double feature[num_chosen_nodes][num_hops*3];
    for(int i=0;i<num_nodes;i++){
        degree[i]=0;
        graph[i].clear();
    }
    for (int i=0;i<num_chosen_nodes;i++){
    	for (int len=0;len<num_hops*3;len++)
        	feature[i][len]=0;
    }

    ifstream in;
    in.open(node_file.c_str(),ifstream::in);
    for (int i=0; i<num_chosen_nodes; i++){
        in >> nodelist[i];
    }
    in.close();
    cout << "initialized." << endl;

    generate_social_topology();
	cout << "social network topology generated." << endl;

	set<int>::iterator iter, iter2, iter3;
	set<int> neigh;
	int total_edges;
	for(int i=0;i<num_chosen_nodes;i++){
		int node=nodelist[i];
		for (int k=1; k<=num_hops; k++){
			neigh.clear();
			neigh.insert(node);
			total_edges=0;
			get_k_hop_neighbours(node, k, neigh);
			feature[i][k*3-3]=neigh.size(); // num of k-hop neighbours
			for (iter=neigh.begin();iter!=neigh.end();iter++) {
				feature[i][k*3-2]+=degree[*iter]; // k-hop neighbour total degree
				// for (iter2=neigh.begin();iter2!=neigh.end();iter2++) {
				// 	iter3=graph[*iter].find(*iter2);
				// 	if (iter3!=graph[*iter].end()) total_edges++;
				// }
				// if (feature[i][k*3-3]==1) {
				// 	feature[i][k*3-1]=1;
				// }
				// else
				// 	feature[i][k*3-1]=total_edges/feature[i][k*3-3]/(feature[i][k*3-3]-1);
			}
			// feature[i][k*3-2]/=feature[i][k*3-3]; // k-hop neighbour clustering coefficient
		}
	}

	ofstream out;
    out.open(feature_file.c_str(),ofstream::out);
    if(out.fail()){
        printf("2. error opening file\n");
        exit(1);
    }
    for (int i=0; i<num_chosen_nodes; i++) {
    	out << nodelist[i] << "\t";
        for (int len=0; len<num_hops*3; len++) {
            out << feature[i][len] << "\t";
        }
        out<<endl;
    }
    out.close();

    std::cout << "Time elapsed: " << t.elapsed() << "sec" << std::endl;

	return 0;
}

void generate_social_topology(){

	// use the real social network topology

	ifstream in;
	in.open(graph_file.c_str(),ifstream::in);
	if(in.fail()){
		printf("error opening file\n");
		exit(1);
	}

	long int i,j;
	string line;

	getline(in, line);
    getline(in, line);
    getline(in, line);
    while (!in.eof()) {
        getline(in, line);
        stringstream ss(line);
        if (ss) {
            ss >> i;
            if (ss.fail()) break;
            ss >> j;
            if (ss.fail()) break;
            degree[i]++;
            degree[j]++;
            graph[i].insert(j);
            graph[j].insert(i);
        }
    }
    in.close();
}

void get_k_hop_neighbours(int node, int k, set<int> &neigh) {
	set<int>::iterator iter;
	// cout<<"Node: "<<node<<endl;
	if (k==1) {
		for (iter=graph[node].begin();iter!=graph[node].end();iter++){
			neigh.insert(*iter);
			// cout<<*iter<<"\t";
		}
		// cout<<endl;
	}
	else {
		for (iter=graph[node].begin();iter!=graph[node].end();iter++){
			get_k_hop_neighbours(*iter, k-1, neigh);
		}
	}

}



