clear all; clc;
fid=fopen('../datasets/graph.txt','r');
fgetl(fid);
fgetl(fid);
fgetl(fid);
graph=fscanf(fid,'%d %d',[2 Inf]);
graph=graph';
graph=graph+1';
row=size(graph,1);
label=1:max(max(graph));
graph=[graph;graph(:,[2,1])];
graph(:,1),I]=sort(graph(:,1));
graph(:,2)=graph(I,2);

nodedegree=generateTranMatrix(graph);

M=30;


load('../results/walk_length_thres.mat','len_all');
len=len_all;

indices=[75,50,25];
for kk=1:3
    ind=indices(kk);
tic
Graph = [];
total_len=[];
for  i = 1:length(nodedegree) 
    count = 1;
    pos_nei=find(graph(:,1)==label(i));
    for j = 1:nodedegree(i)
        loop = 0;
        edge_exist = 1; 
        
%         fprintf('%d\n',t);
        while (loop <= M && edge_exist)
            tmp_len=1;
            now = graph(pos_nei(j),2);
            t=len(now,kk);
            while(t>0)
                tmp_len=tmp_len+1;
                pos = find (graph(:,1) == now) ;
                degree=nodedegree(find(label==now));
                now = graph(pos(randi([1,degree],1)),2);
                t=min([t-1;len(now,kk)]);
            end
            loop = loop + 1;
            if isempty(Graph) && now == label(i)
                edge_exist = 1;
            else if isempty(Graph) && now ~= label(i)
                    edge_exist=0;
                else if now == label(i) ||  ~isempty(find (Graph(find( Graph(:,1) == label(i)),2) == now))
                        edge_exist=1;
                    else if isempty(find (Graph(find( Graph(:,1) == label(i)),2) == now)) && now ~= label(i)
                        edge_exist=0;
                        end
                    end
                end
            end
        end
        if (loop <= M)
            total_len=[total_len;tmp_len];
            if count == 1
                Graph = [Graph; label(i), now ; now,label(i)];
            else
                p = (0.5*nodedegree(i)-1)/(nodedegree(i)-1);
                if rand(1) <= p
                    Graph = [Graph; label(i), now; now,label(i)];
                end
            end
            count= count +1;
        end
    end
    
end
Graph_new=Graph;
fprintf('avg walk length = %f\n',sum(total_len)/length(total_len));

fname=sprintf('./perturbation/graph/newGraph_adaptive_path.%d.txt',ind);
fid = fopen(fname,'w');
for i=1:size(Graph_new,1);
    fprintf(fid, '%d\t%d\n', Graph_new(i,1),Graph_new(i,2));
end
toc
end