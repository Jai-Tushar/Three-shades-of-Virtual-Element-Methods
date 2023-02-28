function [ar,cen,hD,G,B,D,H] = localGBDHp3R1(~,lc4n)

% ------------Test meshes----------------
% lc4n = [0 0; 3 0; 3 2; 3/2 4; 0 4]; 
% ln4e = [1 2 3 4 5]; ln4e = num2cell(ln4e,2);
% lc4n = [0 0; 0.1250 0.1250; 0 0.1250];
% ln4e = [1 2 3];  ln4e = num2cell(ln4e,2);
% lc4n = [0 0; 0.5 0; 0.5 0.5; 0 0.5];
% ln4e = [1 2 3 4]; ln4e = num2cell(ln4e,2);

%-----------Element Geometry-------------
Ne = length(lc4n);
ar_components = lc4n(:,1) .* lc4n([2:Ne,1],2) - lc4n([2:Ne,1],1) .* lc4n(:,2);
ar = 0.5 * abs(sum(ar_components));                                         % area
cen = sum((lc4n + lc4n([2:Ne,1],:)) .* repmat(ar_components,1,2))/(6*ar);   % centroid
hD = 0;                                                                     % diameter
 for i = 1:(Ne-1)                                              
    for j = (i+1):Ne                                 
        hD = max(hD, norm(lc4n(i,:) - lc4n(j,:)));    
    end                                                    
 end

 clear i j
 
%-------------Dimensions----------------
nk = 10;
vertexDof = Ne;
edgeDof = Ne*2;
intDof = 3;
Ndof = vertexDof+edgeDof+intDof;

%-------------Monomial Basis----------------------
m1 = @(x)(1); m2 = @(x)((x(:,1)-cen(1))./hD); m3 = @(x)((x(:,2)-cen(2))./hD);
m4 = @(x)(((x(:,1)-cen(1))./hD).^2); m5 = @(x)(((x(:,1)-cen(1))./hD).*((x(:,2)-cen(2))./hD));
m6 = @(x)(((x(:,2)-cen(2))./hD).^2); m7 = @(x)(((x(:,1)-cen(1))./hD).^3);
m8 = @(x)((((x(:,1)-cen(1))./hD).^2).*((x(:,2)-cen(2))./hD));
m9 = @(x)(((x(:,1)-cen(1))./hD).*(((x(:,2)-cen(2))./hD).^2));
m10 = @(x)(((x(:,2)-cen(2))./hD).^3);

%----------Gradient of Monomial Basis -------------------
gm1 = @(x) [0 0];   gm2 = @(x) [1 0]*(1/hD);    gm3 = @(x) [0 1]*(1/hD);           
gm4 = @(x) [2*m2(x) 0]*(1/hD);  gm5 = @(x) [m3(x) m2(x)]*(1/hD);   
gm6 = @(x) [0 2*m3(x)]*(1/hD);  gm7 = @(x) [3*m4(x) 0]*(1/hD);
gm8 = @(x) [2*m5(x) m4(x)]*(1/hD);  gm9 = @(x) [m6(x) 2*m5(x)]*(1/hD);
gm10 = @(x) [0 3*m6(x)]*(1/hD);

%------------Normals----------------------
ln4e = 1:length(lc4n); ln4e = num2cell(ln4e,2);
temp = circshift(ln4e{:},length(ln4e{:})-1)';
edgeCo = [ln4e{:}' temp];
v1 = 1:Ne; v2 = [2:Ne,1];                                                   % loop index for vertices or edges
% me = subedgeRoutine(1,lc4n);                                              % Ordering edge Dof: equidistant 
me = subedgeRoutineGLob(3,lc4n,edgeCo);                                     % Ordering edge Dof: Gauss Lobatto Points 
edgeNo = lc4n(v2,:)-lc4n(v1,:);
edgeLe = sqrt(edgeNo(:,1).^2+edgeNo(:,2).^2); 
edgeUnv = [edgeNo(:,2) -edgeNo(:,1)].*(1./edgeLe); 

%------------Triangulation----------------------
nodeT = [lc4n(ln4e{:},:); cen];
elemT = [(Ne+1)*ones(Ne,1), (1:Ne)', [2:Ne,1]'];

%------------B Matrix-----------------
B = zeros(nk,Ndof);

B(1,Ndof-2) = 1;

% Internal Dof
B(4,Ndof-2) = -(2/hD)*(1/hD)*ar;
B(6,Ndof-2) = -(2/hD)*(1/hD)*ar;
B(7,Ndof-1) = -(6/hD)*(1/hD)*ar; 
B(8,Ndof)   = -(2/hD)*(1/hD)*ar;
B(9,Ndof-1) = -(2/hD)*(1/hD)*ar;
B(10,Ndof)  = -(6/hD)*(1/hD)*ar;

% Edge Dof
mod_wrap = @(x,a)mod(x-1,a) + 1;
meEdge = me(vertexDof+1:end,:);
edgeLe4Dof = repmat(edgeLe,2,1);
edgeUnvDof = repmat(edgeUnv,2,1);
phiEdge   = eye(edgeDof);

gm = @(x) [gm1(x); gm2(x); gm3(x); gm4(x); gm5(x); gm6(x); gm7(x); gm8(x); gm9(x); gm10(x)];

for i = 1:edgeDof
    ma = gm(meEdge(mod_wrap(i,Ne),:));
    mb = gm(meEdge(mod_wrap(i,Ne)+Ne,:));
    B(:,vertexDof+i) = B(:,vertexDof+i) + 0.5*edgeLe4Dof(i).*...            % Gauss Lobatto Quadrature
                       ((5/6)*(sum(ma.*edgeUnvDof(mod_wrap(i,Ne),:),2)).*phiEdge(mod_wrap(i,Ne),i) + ...
                        (5/6)*(sum(mb.*edgeUnvDof(mod_wrap(i,Ne)+Ne,:),2)).*phiEdge(mod_wrap(i,Ne)+Ne,i));
end

clear i ma mb

% Vertex Dof
meVertex = me(1:vertexDof,:);

for i = 1:vertexDof
    ma = gm(meVertex(i,:));
    mb = gm(meVertex(i,:));
    B(:,i) = B(:,i) + 0.5*edgeLe(mod_wrap(i-1,Ne)).*...                     % Gauss Lobatto Quadrature
                     ((1/6)*(sum(ma.*edgeUnv(mod_wrap(i-1,Ne),:),2))) + ...
                      0.5*edgeLe(i).*((1/6)*(sum(mb.*edgeUnv(i,:),2)));
end

%-----------H Matrix-----------------
intH = @(x) [m1(x)*m1(x) m1(x)*m2(x) m1(x)*m3(x) m1(x)*m4(x) m1(x)*m5(x) m1(x)*m6(x) m1(x)*m7(x) m1(x)*m8(x) m1(x)*m9(x) m1(x)*m10(x);
             m2(x)*m1(x) m2(x)*m2(x) m2(x)*m3(x) m2(x)*m4(x) m2(x)*m5(x) m2(x)*m6(x) m2(x)*m7(x) m2(x)*m8(x) m2(x)*m9(x) m2(x)*m10(x);
             m3(x)*m1(x) m3(x)*m2(x) m3(x)*m3(x) m3(x)*m4(x) m3(x)*m5(x) m3(x)*m6(x) m3(x)*m7(x) m3(x)*m8(x) m3(x)*m9(x) m3(x)*m10(x);
             m4(x)*m1(x) m4(x)*m2(x) m4(x)*m3(x) m4(x)*m4(x) m4(x)*m5(x) m4(x)*m6(x) m4(x)*m7(x) m4(x)*m8(x) m4(x)*m9(x) m4(x)*m10(x);
             m5(x)*m1(x) m5(x)*m2(x) m5(x)*m3(x) m5(x)*m4(x) m5(x)*m5(x) m5(x)*m6(x) m5(x)*m7(x) m5(x)*m8(x) m5(x)*m9(x) m5(x)*m10(x);
             m6(x)*m1(x) m6(x)*m2(x) m6(x)*m3(x) m6(x)*m4(x) m6(x)*m5(x) m6(x)*m6(x) m6(x)*m7(x) m6(x)*m8(x) m6(x)*m9(x) m6(x)*m10(x);
             m7(x)*m1(x) m7(x)*m2(x) m7(x)*m3(x) m7(x)*m4(x) m7(x)*m5(x) m7(x)*m6(x) m7(x)*m7(x) m7(x)*m8(x) m7(x)*m9(x) m7(x)*m10(x);
             m8(x)*m1(x) m8(x)*m2(x) m8(x)*m3(x) m8(x)*m4(x) m8(x)*m5(x) m8(x)*m6(x) m8(x)*m7(x) m8(x)*m8(x) m8(x)*m9(x) m8(x)*m10(x);
             m9(x)*m1(x) m9(x)*m2(x) m9(x)*m3(x) m9(x)*m4(x) m9(x)*m5(x) m9(x)*m6(x) m9(x)*m7(x) m9(x)*m8(x) m9(x)*m9(x) m9(x)*m10(x);
             m10(x)*m1(x) m10(x)*m2(x) m10(x)*m3(x) m10(x)*m4(x) m10(x)*m5(x) m10(x)*m6(x) m10(x)*m7(x) m10(x)*m8(x) m10(x)*m9(x) m10(x)*m10(x);];

% Quadrature
H = integralTri(intH,6,nodeT,elemT);

%-----------D Matrix-----------------
D = zeros(Ndof,nk);

m = @(x) [m1(x) m2(x) m3(x) m4(x) m5(x) m6(x) m7(x) m8(x) m9(x) m10(x)];

%Boundary Dof
for i = 1:vertexDof+edgeDof
    D(i,:) = D(i,:) + m(me(i,:));
end

clear i

% Internal Dof
D(vertexDof+edgeDof+1:end,:) = (1/ar)*H(1:3,:);


%------------G Matrix--------------------
% intG = @(x) [dot(gm1(x),gm1(x)) dot(gm1(x),gm2(x)) dot(gm1(x),gm3(x)) dot(gm1(x),gm4(x)) dot(gm1(x),gm5(x)) dot(gm1(x),gm6(x)) dot(gm1(x),gm7(x)) dot(gm1(x),gm8(x)) dot(gm1(x),gm9(x)) dot(gm1(x),gm10(x));
%              dot(gm2(x),gm1(x)) dot(gm2(x),gm2(x)) dot(gm2(x),gm3(x)) dot(gm2(x),gm4(x)) dot(gm2(x),gm5(x)) dot(gm2(x),gm6(x)) dot(gm2(x),gm7(x)) dot(gm2(x),gm8(x)) dot(gm2(x),gm9(x)) dot(gm2(x),gm10(x));
%              dot(gm3(x),gm1(x)) dot(gm3(x),gm2(x)) dot(gm3(x),gm3(x)) dot(gm3(x),gm4(x)) dot(gm3(x),gm5(x)) dot(gm3(x),gm6(x)) dot(gm3(x),gm7(x)) dot(gm3(x),gm8(x)) dot(gm3(x),gm9(x)) dot(gm3(x),gm10(x));
%              dot(gm4(x),gm1(x)) dot(gm4(x),gm2(x)) dot(gm4(x),gm3(x)) dot(gm4(x),gm4(x)) dot(gm4(x),gm5(x)) dot(gm4(x),gm6(x)) dot(gm4(x),gm7(x)) dot(gm4(x),gm8(x)) dot(gm4(x),gm9(x)) dot(gm4(x),gm10(x));
%              dot(gm5(x),gm1(x)) dot(gm5(x),gm2(x)) dot(gm5(x),gm3(x)) dot(gm5(x),gm4(x)) dot(gm5(x),gm5(x)) dot(gm5(x),gm6(x)) dot(gm5(x),gm7(x)) dot(gm5(x),gm8(x)) dot(gm5(x),gm9(x)) dot(gm5(x),gm10(x));
%              dot(gm6(x),gm1(x)) dot(gm6(x),gm2(x)) dot(gm6(x),gm3(x)) dot(gm6(x),gm4(x)) dot(gm6(x),gm5(x)) dot(gm6(x),gm6(x)) dot(gm6(x),gm7(x)) dot(gm6(x),gm8(x)) dot(gm6(x),gm9(x)) dot(gm6(x),gm10(x));
%              dot(gm7(x),gm1(x)) dot(gm7(x),gm2(x)) dot(gm7(x),gm3(x)) dot(gm7(x),gm4(x)) dot(gm7(x),gm5(x)) dot(gm7(x),gm6(x)) dot(gm7(x),gm7(x)) dot(gm7(x),gm8(x)) dot(gm7(x),gm9(x)) dot(gm7(x),gm10(x));
%              dot(gm8(x),gm1(x)) dot(gm8(x),gm2(x)) dot(gm8(x),gm3(x)) dot(gm8(x),gm4(x)) dot(gm8(x),gm5(x)) dot(gm8(x),gm6(x)) dot(gm8(x),gm7(x)) dot(gm8(x),gm8(x)) dot(gm8(x),gm9(x)) dot(gm8(x),gm10(x));
%              dot(gm9(x),gm1(x)) dot(gm9(x),gm2(x)) dot(gm9(x),gm3(x)) dot(gm9(x),gm4(x)) dot(gm9(x),gm5(x)) dot(gm9(x),gm6(x)) dot(gm9(x),gm7(x)) dot(gm9(x),gm8(x)) dot(gm9(x),gm9(x)) dot(gm9(x),gm10(x));
%              dot(gm10(x),gm1(x)) dot(gm10(x),gm2(x)) dot(gm10(x),gm3(x)) dot(gm10(x),gm4(x)) dot(gm10(x),gm5(x)) dot(gm10(x),gm6(x)) dot(gm10(x),gm7(x)) dot(gm10(x),gm8(x)) dot(gm10(x),gm9(x)) dot(gm10(x),gm10(x))];
% 
% m = @(x) [m1(x) m2(x) m3(x) m4(x) m5(x) m6(x) m7(x) m8(x) m9(x) m10(x)];
%          
% % Quadrature      
% G = integralTri(intG,6,nodeT,elemT);         
% G(1,:) = (1/ar)*integralTri(m,6,nodeT,elemT);
G = B*D;


