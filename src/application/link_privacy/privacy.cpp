#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <boost/timer.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <set>
#include <omp.h>
#include <string>
#include <time.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

# define estimation_of_entries 5e+05
# define num_chosen_nodes 20

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
typedef pair<int, int> Edge;


string graph_file, graph_new_file, privacy_file, privacy_file2;
int num_nodes, walk_len;
double num_edges, num_edges2, num_edges_new;  // num of edges * 2
int num_chosen_edges=100;


void read_network(string file_name, vector<Edge> &edgelist);
void generate_network(vector<Edge> &edgelist, SpMat &mat, double* degree);	// uses real social network topology
int randomInt(int range); // random integer between [0, range]

int main(int argc, const char * argv[]) {

    if(argc==4){
        graph_file="./graphs/graph_"+string(argv[1])+".txt";
        graph_new_file="./graphs/newGraph_"+string(argv[1])+".txt";
        privacy_file="./results/privacy_"+string(argv[1])+".txt";
        privacy_file2="./results/privacy2_"+string(argv[1])+".txt";
        num_nodes=atol(argv[2]);
        walk_len=atol(argv[3]);
        printf("num_nodes = %d\nwalk_len = %d\n", num_nodes, walk_len);
    }
    else{
        printf("Usage: $ ./%s graph num_nodes walk_len\n", argv[0]);
        exit(1);
    }

    boost::timer t;
    double *degree, *degree2, *degree_new;     // node degree
    vector<Edge> edgelist, edgelist2, edgelist_new;
    SpMat mat(num_nodes,num_nodes), mat2(num_nodes,num_nodes), mat_new(num_nodes,num_nodes);    // transition matrx

    srand((unsigned)time(NULL));

    // allocate memory
    degree=new double[num_nodes];
    degree2=new double[num_nodes];
    degree_new=new double[num_nodes];

    read_network(graph_file, edgelist); // G
    read_network(graph_new_file, edgelist_new); // G'
    
    generate_network(edgelist, mat, degree); // G
    generate_network(edgelist_new, mat_new, degree_new); // G'
    
    num_edges=0;
    for(int i=0;i<num_nodes;i++){
        num_edges+=degree[i];
    }
    // printf("num_edges=%f\n", num_edges/2);
    printf("num_edges=%d\n", edgelist.size());
    num_edges_new=0;
    for(int i=0;i<num_nodes;i++){
        num_edges_new+=degree_new[i];
    }
    // printf("num_edges_new=%f\n", num_edges_new/2);
    printf("num_edges_new=%d\n", edgelist_new.size());

    ofstream out,out2;
    out.open(privacy_file.c_str(),ofstream::out);
    out2.open(privacy_file2.c_str(),ofstream::out);
    // for chosen (p,q) in G
    for (int i=0; i<num_chosen_edges; i++){
        cout<<"edge = "<<i<<".\t";
        int random=randomInt(edgelist.size());
        int p,q;
        p=edgelist_new[random].first;
        q=edgelist_new[random].second;
        cout<<"("<<p<<",\t"<<q<<")\n"<<endl;
        // cout<<"degree: "<<degree[p]<<",\t"<<degree[q]<<endl;

        if (degree[p]==1 || degree[q]==1) continue;

        // generate G-L;
        double prob, ratio;
        int min=p<q? p:q;
        int max=p>q? p:q;
        vector<int> neighbors_min, neighbors_max;
        int common=0;

        edgelist2.clear();
        neighbors_min.clear();
        neighbors_max.clear();
        for (vector<Edge>::iterator iter=edgelist.begin();iter!=edgelist.end();iter++){
            Edge edge=*iter;
            int a=edge.first, b=edge.second;
            if (a==min){
                if (b!=max) neighbors_min.push_back(b);
            }
            if (a==max){
                if (b!=min) neighbors_max.push_back(b);
                for (vector<int>::iterator iter2=neighbors_min.begin();iter2!=neighbors_min.end();iter2++){
                    if (b==*iter2) common++;
                }
            }
            if ((a==p && b==q) ||  (a==q && b==p)) continue; 
            edgelist2.push_back(edge);
        }
        // calculate the number of common neighbours
        double Jacc=(double)common/(degree[p]+degree[q]-2-common);
        if (Jacc<0.1) ratio=10;
        else if (Jacc<0.2) ratio=5;
        else if (Jacc<0.3) ratio=2;
        else if (Jacc<0.4) ratio=1;
        else if (Jacc<0.5) ratio=0.5;
        else if (Jacc<0.75) ratio=0.2;
        else ratio=0.1;
        ratio=(1-Jacc)/Jacc;

        generate_network(edgelist2, mat2, degree2); // G-L

        // for each (a,b) in G'
        double ratio2=1;
        int cnt=0;
        for (vector<Edge>::iterator iter=edgelist_new.begin();iter!=edgelist_new.end();iter++){
            Edge edge=*iter;
            int a=edge.first, b=edge.second, flag=0;
            for (vector<int>::iterator iter2=neighbors_min.begin();iter2!=neighbors_min.end();iter2++){
                    if (a==*iter2 || b==*iter2) flag=1;
            }
            for (vector<int>::iterator iter2=neighbors_max.begin();iter2!=neighbors_max.end();iter2++){
                    if (a==*iter2 || b==*iter2) flag=2;
            }
            if (a==p || a==q || b==p || b==q) flag=3;
            if (flag==0){
                double r=(double)rand()/(RAND_MAX+1.0);
                if (r>0.0005) continue;
                
            }
            // cout<<"\tflag:"<<flag<<endl;
            // cout<<"\t\t# of edge: "<<++cnt<<endl;

            double prob1, prob2;
            VectorXd state, state2;
            // P^t(G)
            state=VectorXd::Zero(num_nodes);
            state2=VectorXd::Zero(num_nodes);
            state(a)=1;
            state2(b)=1;
            for (int len=0; len<walk_len; len++){
                state=mat * state;
                state2=mat * state2;
            }
            prob1=state(b)+state2(a);
            // P^t(G-L)
            state=VectorXd::Zero(num_nodes);
            state2=VectorXd::Zero(num_nodes);
            state(a)=1;
            state2(b)=1;
            for (int len=0; len<walk_len; len++){
                state=mat2 * state;
                state2=mat2 * state2;
            }
            prob2=state(b)+state2(a);

            ratio2*=prob2/prob1;
            // if (prob2>prob1 || prob2<prob1){
            //     cout<<"\t\tprob2/prob1: "<<prob2/prob1<<endl;
            //     cout<<"\t\tnew ratio: "<<ratio<<endl;
            // }
        }
        prob=1/(1+ratio*ratio2);
        cout<<"prob: "<<prob<<endl;
        out<<prob<<endl;
        if (Jacc==0) Jacc=0.001;
        prob=1/(1+(1-Jacc)*ratio2/Jacc);
        out2<<ratio2<<endl;
    }

    printf("Computation finished.\n");
    out.close();
    out2.close();
    
    std::cout << "Time elapsed: " << t.elapsed() << "sec" << std::endl;
    
    return 0;
}

// read the edgelist from the file
void read_network(string file_name, vector<Edge> &edgelist) {
    edgelist.clear();
    ifstream in;
    in.open(file_name.c_str(),ifstream::in);
    if(in.fail()){
        printf("1. error opening file: %s\n", file_name.c_str());
        exit(1);
    }
    int i,j;
    string line;
    while (!in.eof()) {
        getline(in, line);
        stringstream ss(line);
        if (ss) {
            ss >> i;
            ss >> j;
            edgelist.push_back(make_pair(i-1,j-1)); // renumber the nodes from 0
        }
    }
    in.close();
}

// generate the transition matrix and degree from the edgelist
void generate_network(vector<Edge> &edgelist, SpMat &mat, double* degree){
    std::vector<T> tripletList;
    tripletList.reserve(estimation_of_entries);
    for (int i=0; i<num_nodes; i++){
        degree[i]=0;
    }
    for (vector<Edge>::iterator iter=edgelist.begin();iter!=edgelist.end();iter++){
        Edge edge=*iter;
        degree[edge.first]++; // each bidirectional link appears twice in the edgelist
    }
    for (vector<Edge>::iterator iter=edgelist.begin();iter!=edgelist.end();iter++){
        Edge edge=*iter;
        int i=edge.first, j=edge.second;
        tripletList.push_back(T(i,j,1.0/degree[j]));
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    
}

int randomInt(int range){  // random integer between [0,range)
    int random_no = (int) ((float)range*rand()/(RAND_MAX+1.0));
    return random_no;
}

