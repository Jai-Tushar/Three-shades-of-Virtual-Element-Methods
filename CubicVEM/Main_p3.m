%%%%%%%%%%% Cubic Vem for Poisson problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
% PARAMETERS, DATA and MESH INFO
pde  = LEdata;
NR = 5;
% kk = 2;
l2err = zeros(NR,1);
h1err = zeros(NR,1);
h = zeros(NR,1);

tic
for k = 1:NR

load(['Cp3SqEl_' num2str(k) '.mat']);
load(['Cp3SqDirichlet_' num2str(k) '.mat']);
load(['SquareEl_' num2str(k) '.mat']);
load(['SquareNo_' num2str(k) '.mat']);
N = Node; NT = Element; p3NT = Pkn4e; p3Dirichlet = DirichletPk; 
clear Node Element Pkn4e Dirichlet

% load Testp3PentagonDirichlet.mat
% load Testp3PentagonEl.mat
% N = [0 0; 3 0; 3 2; 3/2 4; 0 4];
% NT = {[1 2 3 4 5]};
% p3NT = Pkn4e; 
% p3Dirichlet = DirichletPk;

% load Testp3SqMeshDirichlet.mat;
% load ATestp3SqMeshEl.mat;
% N = [0 0; 0.5 0; 1 0; 0 0.5; 0.5 0.5; 1 0.5];
% NT = {[1 2 5 4]; [2 3 6 5]};
% p3NT = Pkn4e;
% p3Dirichlet = DirichletPk;

ndof = max(horzcat(p3NT{:})');
nel  = size(p3NT,1);
p3fullnodes = 1:ndof;
p3FreeNodes = setdiff(p3fullnodes,p3Dirichlet);

%% Assebling local Stiffness matrices
[K,F,dia] = Stiffness_p3(N,NT,p3NT,pde);
h(k) = dia;

%% Dirichlet Boundary Condition
u = zeros(ndof,1);
u(p3Dirichlet) = 0;
F = F - K*u;

%% Computation of the solution
u(p3FreeNodes) = K(p3FreeNodes,p3FreeNodes)\F(p3FreeNodes);

%% Error
q1=0;q2=0;
for id = 1:nel
    v_id = NT{id};
    X= N(v_id, :);
    Gl_id = p3NT{id,:};
    [l2e,h1e] = l2h1errp3R1(3,X,u(Gl_id),pde);
    q1=q1+l2e;
    q2=q2+h1e;
end
l2err(k)=sqrt(q1);
h1err(k)=sqrt(q2);
toc
end

l2order = zeros(NR-1,1);
h1order = zeros(NR-1,1);
for j=1:NR-1
    l2order(j) = log(l2err(j)./l2err(j+1))/log(h(j)/h(j+1));
    h1order(j) = log(h1err(j)./h1err(j+1))/log(h(j)/h(j+1));
end
