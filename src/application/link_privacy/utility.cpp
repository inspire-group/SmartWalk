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

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

# define estimation_of_entries 5e+05
# define num_chosen_nodes 20

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


string graph_file, graph_new_file, node_file, utility_file, utility2_file;
long num_nodes;
double num_edges, num_edges2;  // num of edges * 2
int max_walk_len=65;   // max random walk length
int nodes[num_chosen_nodes]; 

void generate_social_topology(string file_name, SpMat &mat, double* degree);	// uses real social network topology

int main(int argc, const char * argv[]) {

    if(argc==3){
        graph_file="./graphs/graph_"+string(argv[1])+".txt";
        graph_new_file="./graphs/newGraph_"+string(argv[1])+".txt";
        node_file="./graphs/nodelist_"+string(argv[1])+".txt";
        utility_file="./results/utility_mean_"+string(argv[1])+".txt";
        utility2_file="./results/utility_max_"+string(argv[1])+".txt";
        num_nodes=atol(argv[2]);
        printf("num_nodes = %d\n", num_nodes);

    }
    else{
        printf("Usage: $ ./%s graph num_nodes\n", argv[0]);
        exit(1);
    }

    boost::timer t;
    double *degree, *degree2;     // node degree
    SpMat mat(num_nodes,num_nodes), mat2(num_nodes,num_nodes);    // transition matrx
    MatrixXd result=MatrixXd::Zero(num_chosen_nodes,max_walk_len); 
    // allocate memory
    degree=new double[num_nodes];
    for(int i=0;i<num_nodes;i++){
        degree[i]=0;
    }
    degree2=new double[num_nodes];
    for(int i=0;i<num_nodes;i++){
        degree2[i]=0;
    }
    
    generate_social_topology(graph_file,mat,degree);
    generate_social_topology(graph_new_file,mat2,degree2);
    cout<<"prob = "<<mat.sum()<<" " <<mat2.sum()<<endl;
    
    num_edges=0;
    for(int i=0;i<num_nodes;i++){
        num_edges+=degree[i];
    }
    printf("num_edges=%f\n", num_edges/2);
    num_edges2=0;
    for(int i=0;i<num_nodes;i++){
        num_edges2+=degree2[i];
    }
    printf("num_edges2=%f\n", num_edges2/2);

    ifstream in;
    in.open(node_file.c_str(),ifstream::in);
    for (int i=0; i<num_chosen_nodes; i++){
        in >> nodes[i];
    }
    in.close();

    ofstream out, out2;
    out.open(utility_file.c_str(),ofstream::out);
    out2.open(utility2_file.c_str(),ofstream::out);

    for (int i=0; i<num_chosen_nodes; i++){
        cout<<"node = "<<i<<endl;
        VectorXd state=VectorXd::Zero(num_nodes);
        VectorXd state2=VectorXd::Zero(num_nodes);
        state(nodes[i])=1;
        state2(nodes[i])=1;
        for (int len=0; len<max_walk_len; len++){
            state=mat * state;
            state2=mat2 * state2;
            // cout<<"sum: "<<state.sum()<<"\t"<<state2.sum()<<endl;
            result(i,len)=(state2-state).cwiseAbs().sum()*0.5;
            // cout<<"length = "<< len+1 <<"\t difference = "<<result(i,len)<<endl;

        }
    }

    for (int len=0; len<max_walk_len; len++){
        double utility=0, utility2=0;
        cout<<"walk length = "<< len+1 <<endl;
        for (int i=0; i<num_chosen_nodes; i++){
            utility+=result(i,len);
            if (result(i,len)>utility2) utility2=result(i,len);
        }
        utility/=num_chosen_nodes;
        cout<<"mean utility = "<< utility <<endl;
        out<< utility <<endl;
        cout<<"max utility = "<< utility2 <<endl;
        out2<< utility2 <<endl;
    }
    printf("Computation finished.\n");
    out.close();
    out2.close();
    
    std::cout << "Time elapsed: " << t.elapsed() << "sec" << std::endl;
    
    return 0;
}

void generate_social_topology(string file_name, SpMat &mat, double *degree) {
    ifstream in;
    in.open(file_name.c_str(),ifstream::in);
    if(in.fail()){
        printf("1. error opening file: %s\n", file_name.c_str());
        exit(1);
    }
    
    long int i,j;
    string line;
    std::vector<T> tripletList;
    tripletList.reserve(estimation_of_entries);

    while (!in.eof()) {
        getline(in, line);
        stringstream ss(line);
        if (ss) {
            ss >> i;
            if (ss.fail()) break;
            degree[i-1]++;
            ss >> j;
        }
    }
    in.close();

    in.open(file_name.c_str(),ifstream::in);
    if(in.fail()){
        printf("1. error opening file: %s\n", file_name.c_str());
        exit(1);
    }

    while (!in.eof()) {
        getline(in, line);
        stringstream ss(line);
        if (ss) {
            ss >> i;
            if (ss.fail()) break;
            ss >> j;
            if (ss.fail()) break;
            tripletList.push_back(T(i-1,j-1,1.0/degree[j-1]));
            // printf("(%d, %d) = %f\n", i-1,j-1,1.0/degree[j-1]);
            // tripletList.push_back(T(j-1,i-1,1.0/degree[i-1]));
            //printf("(%d, %d) = %f\n", j,i,1.0/degree[i]);
        }
    }
    in.close();
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

}




