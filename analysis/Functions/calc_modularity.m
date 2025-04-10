function Q = calc_modularity(A, Ci)
% Part of modularity_und, BCT toolbox

if ~exist('gamma','var')
    gamma = 1;
end

N=length(A);                            %number of vertices
K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B=A-gamma*(K.'*K)/m;                    %modularity matrix

s=Ci(:,ones(1,N));                      %compute modularity
Q=~(s-s.').*B/m;
Q=sum(Q(:));

end