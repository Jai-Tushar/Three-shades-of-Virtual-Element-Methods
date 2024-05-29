%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear VEM to solve 
%      -\nabla . (K \nabla u) = f in \Omega
%                          u  = 0 on \partial \Omega 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jai Tushar, BITS-Pilani.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
%% PRE-PROCESSING
% PARAMETERS, DATA and MESH INFO
pde  = LEdata;
NR = 3;
kk = 5;
l2err = zeros(NR,1);
h1err = zeros(NR,1);
h = zeros(NR,1);

tic
for k = 1:NR
    % CONCAVE
    if kk == 1
        load(['ConcaveEl_' num2str(k) '.mat']);
        load(['ConcaveNo_' num2str(k) '.mat']);
        N = Node; NT = Element;
    end
    % STRUCTURED VORONOI
    if kk == 2
        load(['SVEl_' num2str(k) '.mat']);
        load(['SVNo_' num2str(k) '.mat']);
        N = Node; NT = Element;
    end
    % RANDOM VORONOI
    if kk == 3
        load(['RVEl_' num2str(k) '.mat']);
        load(['RVNo_'    num2str(k) '.mat']);
        N = Node; NT = Element;
    end
    if kk == 4
       load VoronoiElement.mat
       load VoronoiNode.mat
       N = X; NT = Element;
    end
    if kk == 5
        load(['VInt_' num2str(k) '.mat']);
        N = mesh.p; NT = mesh.t;
    end
    %% PROCESSING
    % Assemble Matrices
    [A,M,b,dia,FreeNodes,dirichlet] = StiffMass(N,NT,pde);
    b = b(FreeNodes);
    h(k) = dia;
    ndof = size(N,1);
    nel  = size(NT,1);
    
    % Dirichlet boundary conditions
    uD = zeros(ndof,1);
    uD(dirichlet) = pde.exactu(N(dirichlet,:));
    b = b - A(FreeNodes,:)*uD;
    % Direct solve
    u = uD;
    Z = (A(FreeNodes,FreeNodes)+M(FreeNodes,FreeNodes));
    u(FreeNodes) = Z\b;
    %% POST PROCESSING
    % Visualise
    uex = pde.exactu(N);
    figure(2); clf
    set(gcf,'Units','normal');
    set(gcf,'Position',[0.25,0.25,0.6,0.25]);
    subplot(1,2,1)
    plotsol(NT,N,uex);
    title('Exact solution')
    xlabel('x'); ylabel('y'); zlabel('u');
    axis tight;
    subplot(1,2,2)
    plotsol(NT,N,u);
    title('Approximate solution')
    xlabel('x'); ylabel('y'); zlabel('u');
    axis tight;
    Reorder = [mesh.p(mesh.InterfaceIdx,2) mesh.InterfaceIdx];
    Reorder =sortrows(Reorder);
    x = Reorder(:,1);
    figure(3)
    plot(x,uex(Reorder(:,2)),'--r',x,u(Reorder(:,2)),'b','LineWidth',2);    

%     % Error
%     q1=0;q2=0;
%     for j = 1:nel
%         v_id = NT{j,:};
%         verts= N(v_id(:),:);
%         [l2e,h1e]=err(verts,u(v_id),pde);
%         q1=q1+l2e;
%         q2=q2+h1e;
%     end
%     l2err(k) = sqrt(q1);
%     h1err(k) = sqrt(q2);

u_uex = pde.exactu(N(FreeNodes,:)) - u(FreeNodes,:);
l2err(k) = sqrt((u_uex)'*M(FreeNodes,FreeNodes)*(u_uex));
end
toc

l2order = zeros(NR-1,1);
% h1order = zeros(NR-1,1);
for j=1:NR-1
    l2order(j) = log(l2err(j)./l2err(j+1))/log(h(j)/h(j+1));
%     h1order(j) = log(h1err(j)./h1err(j+1))/log(h(j)/h(j+1));
end

% fprintf('h = max(diameter)');
% h
% fprintf('L^2 error and H^1 error');
% [l2err h1err]
% fprintf('L^2 EOC and H^1 EOC');
% [l2order h1order]
