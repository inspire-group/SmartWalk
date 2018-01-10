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
# define num_chosen_nodes 29060

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


string graph_file, node_file, result_file, result_hard_file, mixing_time_file;
long num_nodes;
double num_edges;  // num of edges * 2
int max_walk_len=200;   // max random walk length
double *degree;     // node degree
long nodes[num_chosen_nodes]; // 1000 randomly selected nodes

void generate_social_topology(SpMat &mat);	// uses real social network topology
void output_result(MatrixXd &result, string file_name1, string file_name2);

int main(int argc, const char * argv[]) {

    if(argc==3){
        graph_file="../datasets/ungraph_"+string(argv[1])+".txt";
        node_file="../datasets/nodelist/nodelist_"+string(argv[1])+".txt";
        stringstream ss;
        ss << num_chosen_nodes;
        result_file="../results/result_"+string(argv[1])+"_"+ss.str()+".txt";
        result_hard_file="../results/result_hard_"+string(argv[1])+"_"+ss.str()+".txt";
        mixing_time_file="../results/mixing_time_"+string(argv[1])+"_"+ss.str()+".txt";
        num_nodes=atol(argv[2]);
        printf("num_nodes = %d\n", num_nodes);

    }
    else{
        printf("Usage: $ ./%s graph num_nodes\n", argv[0]);
        exit(1);
    }

    boost::timer t;
    VectorXd curr_state=VectorXd::Zero(num_nodes);    // current state
    VectorXd final_state=VectorXd::Zero(num_nodes);    // final state
    SpMat mat(num_nodes,num_nodes);    // transition matrx
    MatrixXd result=MatrixXd::Zero(num_chosen_nodes,max_walk_len);  // variance distance
    // allocate memory
    degree=new double[num_nodes];
    for(int i=0;i<num_nodes;i++){
        degree[i]=0;
    }
    
    generate_social_topology(mat);
    
    num_edges=0;
    for(int i=0;i<num_nodes;i++){
        num_edges+=degree[i];
    }
    printf("num_edges=%f\n", num_edges/2);
    for(int i=0;i<num_nodes;i++){
        final_state(i)=degree[i]/num_edges;
    }
    ifstream in;
    in.open(node_file.c_str(),ifstream::in);
    for (int i=0; i<num_chosen_nodes; i++){
        in >> nodes[i];
    }
    in.close();
    ofstream out;
    out.open(mixing_time_file.c_str(),ofstream::out);
    for (int i=0; i<num_chosen_nodes; i++){
        curr_state=VectorXd::Zero(num_nodes);
        curr_state(nodes[i])=1;
        result(i,0)=(curr_state-final_state).cwiseAbs().sum()*0.5;
        int len=0;
        while (result(i,len)>0.25 && (len<max_walk_len-1)) {
            len++;
            curr_state=mat * curr_state;
            result(i,len)=(curr_state-final_state).cwiseAbs().sum()*0.5;
        }
        if (result(i,len)<=0.25) {
            cout<<"node "<<i<<" = "<<len<<endl;
            out<<len<<endl;
        }
    }
    printf("Computation finished.\n");
    output_result(result, result_file, result_hard_file);
    
    std::cout << "Time elapsed: " << t.elapsed() << "sec" << std::endl;
    
    return 0;
}

void generate_social_topology(SpMat &mat) {
    ifstream in;
    in.open(graph_file.c_str(),ifstream::in);
    if(in.fail()){
        printf("1. error opening file\n");
        exit(1);
    }
    
    long int i,j;
    string line;
    std::vector<T> tripletList;
    tripletList.reserve(estimation_of_entries);
    
    getline(in, line);
    getline(in, line);
    getline(in, line);
    while (!in.eof()) {
        getline(in, line);
        stringstream ss(line);
        if (ss) {
            ss >> i;
            if (ss.fail()) break;
            degree[i]++;
            ss >> j;
            if (ss.fail()) break;
            degree[j]++;
        }
    }
    in.close();

    in.open(graph_file.c_str(),ifstream::in);
    if(in.fail()){
        printf("1. error opening file\n");
        exit(1);
    }
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
            tripletList.push_back(T(i,j,1.0/degree[j]));
            //printf("(%d, %d) = %f\n", i,j,1.0/degree[j]);
            tripletList.push_back(T(j,i,1.0/degree[i]));
            //printf("(%d, %d) = %f\n", j,i,1.0/degree[i]);
        }
    }
    in.close();
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

}

void output_result(MatrixXd &result, string file_name1, string file_name2) {
    ofstream out;
    out.open(file_name1.c_str(),ofstream::out);
    if(out.fail()){
        printf("2. error opening file\n");
        exit(1);
    }
    for (int i=0; i<num_chosen_nodes; i++) {
        out << i << " ";
        for (int len=0; len<max_walk_len; len++) {
            out << result(i,len) << " ";
        }
        out<<endl;
    }
//    for (int k=0; k<mat.outerSize(); ++k){
//        for (SpMat::InnerIterator it(mat,k); it; ++it)
//        {
//            out<<it.value()<<'\t'<<it.row()<<'\t'<<it.col()<<endl;
//        }
//    }

    out.close();
    
    out.open(file_name2.c_str(),ofstream::out);
    if(out.fail()){
        printf("2. error opening file\n");
        exit(1);
    }
    for (int i=0; i<num_chosen_nodes; i++) {
        if (result(i,max_walk_len-1)>0.25){
            out << i << " ";
        }
    }
    out.close();
    
}




