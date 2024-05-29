function [A,M,b,hd,FreeNodes,Dirichlet] = StiffMass(N,NT,pde)

% Dirichlet nodes for a unit square domain 
bd = zeros(max(size(N)),min(size(N)));
for i = 1:max(size(N))
    for j = 1:min(size(N))
        if N(i,j) <= 0.1e-8 || N(i,j) >= .99999
            bd(i,j) = i;
        else
            bd(i,j) = 0;
        end
    end
end
Dirichlet = zeros(size(N,1),1);
for i = 1:max(size(N))
    if bd(i,1) == bd(i,2)
        Dirichlet(i,1) = bd(i,1);
    else
        Dirichlet(i,1) = bd(i,1) + bd(i,2);
    end
end
Dirichlet(Dirichlet == 0) = [];

ndof = size(N,1);
nel  = size(NT,1);
fullnodes = 1:ndof;
FreeNodes = setdiff(fullnodes,Dirichlet);

A = sparse(ndof,ndof);
M = sparse(ndof,ndof);
b = zeros(ndof,1);
Di = zeros(nel,1);

for el_id = 1:nel
    vert_id = NT{el_id};
    X   = N(vert_id,:);
    [~,~,hd,G,B,D,H] = localGBDH(X);
    n_sides = length(X);

    h1_p = (G\eye(size(G)))*B;
   
    h1_st = (eye(n_sides)- D * h1_p)' * ...
            (eye(n_sides)- D * h1_p);

    Di(el_id,1) = Di(el_id,1) + hd;
    be = loadterm(X,pde.f);
    b(vert_id,1) = b(vert_id,1) + be'; 
    
    % Classical Stabilization
    G_tilda = G; G_tilda(1,:) = 0;
    AK = (h1_p'*G_tilda*h1_p) + h1_st;

    C = H * (G\eye(size(G))) * B;
    l2_p = (H\eye(size(H))) * C;
%     l2_st = (eye(n_sides)- D * l2_p)' * ...
%             (eye(n_sides)- D * l2_p);
    G_bar = H;
    MK = h1_p'*G_bar*h1_p;
    M(vert_id,vert_id) = M(vert_id,vert_id) + MK;
    A(vert_id,vert_id) = A(vert_id,vert_id) + AK;      
    
dbstop if warning
end

hd = max(Di);