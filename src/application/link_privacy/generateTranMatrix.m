function node_degree=generateTranMatrix(graph)

source=graph(:,1);
dest=graph(:,2);

row=size(graph,1);

M=max([max(source) max(dest)]);
N=min([min(source) min(dest)])-1;

node_degree=zeros(1,M-N);


for i=1:row
    node_degree(source(i)-N)=node_degree(source(i)-N)+1;
%     node_degree(dest(i)-N)=node_degree(dest(i)-N)+1;
end
