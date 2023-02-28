%% local G, B, D and H for k = 3
% Author : Jai Tushar

function [area,cen,d,G,B,D,H] = localGBDH_p3(X)
% load 'test_pentagon_node_p3.mat'   % test element
% X = [0 0; 0.5 0; 0.5 0.5; 0 0.5; 0.1382 0; 0.5 0.1382; ...
%      0.1382 0.5; 0 0.1382; 0.3618 0; 0.5 0.3618; 0.3618 0.5; ...
%      0 0.3618; 0.25 0.25; 0.2535 0.2535; 0.2571 0.2571];

%% area, centroid and diameter for the p3 element
dof = length(X);   % degrees of freedom
nv  = (dof - 3)/3; % vertices of the polygon
newX = X(1:nv,:);
n_sides = length(newX);

area_components = newX(:,1) .* newX([2:nv,1],2) - newX([2:nv,1],1) .* newX(:,2);
area = 0.5 * abs(sum(area_components)); 
cen = sum((newX + newX([2:nv,1],:)) .* repmat(area_components,1,2))/(6*area);
d = 0;                                            
 for i = 1:(nv-1)                                              
    for j = (i+1):nv                                 
        d = max(d, norm(newX(i,:) - newX(j,:)));    
    end                                                    
 end

%% Weights for Gauss quadrature over a triangle. 
xw=TriGaussPoints(8);

%% Local G Matrix
n_polys = 10;
G = zeros(n_polys,n_polys);
Co = zeros(nv,6);
loc_con = zeros(nv,6);
for i = 1:nv
    if i < nv
        loc_con(i,:) = [newX(i,:) newX(i+1,:) cen];
    else
        loc_con(i,:) = [newX(nv,:) newX(1,:) cen];
    end
end

z1 = 0; z2 = 0; z3 = 0; z4 = 0; z5 = 0;
z6 = 0; z7 = 0; z8 = 0; z9 = 0;
z10= 0; z11= 0; z12= 0; z13= 0; z14= 0;
for k = 1:nv
    Co(k,:) = loc_con(k,:);
    At = (1/2)*abs(Co(k,1)*(Co(k,4)-Co(k,6))+Co(k,3)*(Co(k,6)-Co(k,2))...
        + Co(k,5)*(Co(k,2)-Co(k,4)));
    
    int_1 = 0; int_2 = 0; int_3 = 0; int_4 = 0; int_5 = 0;
    int_6 = 0; int_7 = 0; int_8 = 0; int_9 = 0;
    int_10= 0; int_11= 0; int_12= 0; int_13= 0; int_14= 0; 
    for j = 1:length(xw(:,1))
        %transformed coordinates
        x = Co(k,1)*(1-xw(j,1)-xw(j,2))+Co(k,3)*xw(j,1)+Co(k,5)*xw(j,2);
        y = Co(k,2)*(1-xw(j,1)-xw(j,2))+Co(k,4)*xw(j,1)+Co(k,6)*xw(j,2);
        
        int_1 = int_1 + xw(j,3)*((x-cen(1))/d);
        int_2 = int_2 + xw(j,3)*((y-cen(2))/d);
        int_3 = int_3 + xw(j,3)*((x-cen(1))/d)*((x-cen(1))/d);
        int_4 = int_4 + xw(j,3)*((x-cen(1))/d)*((y-cen(2))/d);
        int_5 = int_5 + xw(j,3)*((y-cen(2))/d)*((y-cen(2))/d);
        int_6 = int_6 + xw(j,3)*((x-cen(1))/d)*((x-cen(1))/d)*((x-cen(1))/d);
        int_7 = int_7 + xw(j,3)*((x-cen(1))/d)*((x-cen(1))/d)*((y-cen(2))/d);
        int_8 = int_8 + xw(j,3)*((x-cen(1))/d)*((y-cen(2))/d)*((y-cen(2))/d);
        int_9 = int_9 + xw(j,3)*((y-cen(2))/d)*((y-cen(2))/d)*((y-cen(2))/d);
        int_10= int_10+ xw(j,3)*((x-cen(1))/d)*((x-cen(1))/d)*((x-cen(1))/d)*((x-cen(1))/d);    % a1
        int_11= int_11+ xw(j,3)*((x-cen(1))/d)*((x-cen(1))/d)*((x-cen(1))/d)*((y-cen(2))/d);    % a2
        int_12= int_12+ xw(j,3)*((x-cen(1))/d)*((x-cen(1))/d)*((y-cen(2))/d)*((y-cen(2))/d);    % a3
        int_13= int_13+ xw(j,3)*((x-cen(1))/d)*((y-cen(2))/d)*((y-cen(2))/d)*((y-cen(2))/d);    % a4
        int_14= int_14+ xw(j,3)*((y-cen(2))/d)*((y-cen(2))/d)*((y-cen(2))/d)*((y-cen(2))/d);    % a5
    end
    z1 = z1 + At*int_1;    % int of m2 over the polygon
    z2 = z2 + At*int_2;    % int of m3 over the polygon
    z3 = z3 + At*int_3;    % int of m4 over the polygon
    z4 = z4 + At*int_4;    % int of m5 over the polygon
    z5 = z5 + At*int_5;    % int of m6 over the polygon
    z6 = z6 + At*int_6;    % int of m7 over the polygon
    z7 = z7 + At*int_7;    % int of m8 over the polygon
    z8 = z8 + At*int_8;    % int of m9 over the polygon
    z9 = z9 + At*int_9;    % int of m10 over the polygon
    z10= z10+ At*int_10;   % int of a1 over the polygon
    z11= z11+ At*int_11;   % int of a2 over the polygon
    z12= z12+ At*int_12;   % int of a3 over the polygon
    z13= z13+ At*int_13;   % int of a4 over the polygon
    z14= z14+ At*int_14;   % int of a5 over the polygon
end
G(1,1) = 1;             G(1,6) = z5*(1/area);
G(1,2) = z1*(1/area);   G(1,7) = z6*(1/area);
G(1,3) = z2*(1/area);   G(1,8) = z7*(1/area);
G(1,4) = z3*(1/area);   G(1,9) = z8*(1/area);
G(1,5) = z4*(1/area);   G(1,10)= z9*(1/area);

G(2,2) = (1/d)*(1/d)*area;  G(2,7) = (3/d)*(1/d)*z3;
G(2,4) = 2*(1/d)*(1/d)*z1;  G(2,8) = (2/d)*(1/d)*z4;
G(2,5) = (1/d)*(1/d)*z2;    G(2,9) = (1/d)*(1/d)*z5;    

G(3,3) = G(2,2);            G(3,8) = (1/d)*(1/d)*z3; 
G(3,5) = (1/d)*(1/d)*z1;    G(3,9) = (2/d)*(1/d)*z4;
G(3,6) = 2*(1/d)*(1/d)*z2;  G(3,10)= (3/d)*(1/d)*z5;

G(4,2) = G(2,4);            G(4,7) = (6/d)*(1/d)*z6; 
G(4,4) = 4*(1/d)*(1/d)*z3;  G(4,8) = (4/d)*(1/d)*z7;
G(4,5) = 2*(1/d)*(1/d)*z4;  G(4,9) = (2/d)*(1/d)*z8;

G(5,2) = G(2,5);                  G(5,7) = (3/d)*(1/d)*z7;  
G(5,3) = G(3,5);                  G(5,8) = (2/d)*(1/d)*z8 + (1/d)*(1/d)*z6;  
G(5,4) = G(4,5);                  G(5,9) = (1/d)*(1/d)*z9 + (2/d)*(1/d)*z7;
G(5,5) = (1/d)*(1/d)*(z3 + z5);   G(5,10)= (3/d)*(1/d)*z8;  
G(5,6) = 2*(1/d)*(1/d)*z4; 

G(6,3) = G(3,6);            G(6,8) = (2/d)*(1/d)*z7;
G(6,5) = G(5,6);            G(6,9) = (4/d)*(1/d)*z8;    
G(6,6) = 4*(1/d)*(1/d)*z5;  G(6,10)= (6/d)*(1/d)*z9;

G(7,2) = G(2,7);    G(7,7) = (9/d)*(1/d)*z10;  
G(7,3) = G(3,7);    G(7,8) = (6/d)*(1/d)*z11;
G(7,4) = G(4,7);    G(7,9) = (3/d)*(1/d)*z12;    
G(7,5) = G(5,7);
G(7,6) = G(6,7);

G(8,2) = G(2,8);    G(8,7) = G(7,8);
G(8,3) = G(3,8);    G(8,8) = (4/d)*(1/d)*z12 + (1/d)*(1/d)*z10;
G(8,4) = G(4,8);    G(8,9) = (2/d)*(1/d)*z13 + (2/d)*(1/d)*z11;
G(8,5) = G(5,8);    G(8,10)= (3/d)*(1/d)*z12;
G(8,6) = G(6,8);

G(9,2) = G(2,9);    G(9,7) = G(7,9);
G(9,3) = G(3,9);    G(9,8) = G(8,9);
G(9,4) = G(4,9);    G(9,9) = (1/d)*(1/d)*z14 + (4/d)*(1/d)*z12;
G(9,5) = G(5,9);    G(9,10)= (6/d)*(1/d)*z13;
G(9,6) = G(6,9);

G(10,2) = G(2,10);  G(10,7) = G(7,10);
G(10,3) = G(3,10);  G(10,8) = G(8,10);
G(10,4) = G(4,10);  G(10,9) = G(9,10);
G(10,5) = G(5,10);  G(10,10)= (9/d)*(1/d)*z14;
G(10,6) = G(6,10);

%% Local B Matrix
B = zeros(n_polys,dof);
B(1,dof-2) = 1;
mod_wrap = @(x,a)mod(x-1,a) + 1;             % utility function for wrapping around a vector
edges = zeros(n_sides,2);
edge_length = zeros(n_sides,1);
for vertex_id = 1:n_sides
    current = newX(mod_wrap(vertex_id, n_sides), :);                             % Calculating matrices
    next = newX(mod_wrap(vertex_id + 1, n_sides), :);
    edges(vertex_id,:) = next - current;
    edge_length(vertex_id) = sqrt((edges(vertex_id,1))^2+(edges(vertex_id,2))^2);
end
normal_vector = [edges(:,2),-edges(:,1)]; 
unit_nv = (1./edge_length).*normal_vector;

for j = 1:n_sides   
            gradm2 = [(1/d) 0];
            B(2,j) = dot(0.5*((1/6).*gradm2),((edge_length(mod_wrap(j-1,n_sides)))*...
                (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                (edge_length(mod_wrap(j,n_sides)))*...
                (unit_nv(mod_wrap(j,n_sides),:)))');
            B(2,nv+j) = dot(0.5*((5/6).*gradm2),((edge_length(mod_wrap(j,n_sides)))*...
                (unit_nv(mod_wrap(j,n_sides),:)))');
            B(2,2*nv+j)= B(2,nv+j);
            
            gradm3 = [0 (1/d)];
            B(3,j) = dot(0.5*((1/6).*gradm3),((edge_length(mod_wrap(j-1,n_sides)))*...
                (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                (edge_length(mod_wrap(j,n_sides)))*...
                (unit_nv(mod_wrap(j,n_sides),:)))');
            B(3,nv+j) = dot(0.5*((5/6).*gradm3),((edge_length(mod_wrap(j,n_sides)))*...
                (unit_nv(mod_wrap(j,n_sides),:)))');
            B(3,2*nv+j)= B(3,nv+j);
                
            gradm4a = [(2/d).*((X(j,1)-cen(1))/d), 0];
            B(4,j) = dot(0.5*((1/6).*gradm4a)...
                 ,((edge_length(mod_wrap(j-1,n_sides)))*...
                 (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                 (edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm4b = [(2/d).*((X(nv+j,1)-cen(1))/d), 0];
             B(4,nv+j) = dot(0.5*((5/6).*gradm4b)...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm4c = [(2/d).*((X(2*nv+j,1)-cen(1))/d), 0];
             B(4,2*nv+j) = dot(0.5*((5/6).*gradm4c)...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');   
             
             gradm5a = [(1/d).*((X(j,2)-cen(2))/d), (1/d).*((X(j,1)-cen(1))/d)];
             B(5,j) = dot(0.5*((1/6).*gradm5a),((edge_length(mod_wrap(j-1,n_sides)))*...
                 (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                 (edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm5b = [(1/d).*((X(nv+j,2)-cen(2))/d), (1/d).*((X(nv+j,1)-cen(1))/d)];
             B(5,nv+j) = dot((0.5*(5/6).*gradm5b),((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm5c = [(1/d).*((X(2*nv+j,2)-cen(2))/d), (1/d).*((X(2*nv+j,1)-cen(1))/d)];
             B(5,2*nv+j)= dot((0.5*(5/6).*gradm5c),((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             
             gradm6a = [0, (2/d).*((X(j,2)-cen(2))/d)];   
             B(6,j) = dot((0.5*((1/6).*gradm6a))...
                 ,((edge_length(mod_wrap(j-1,n_sides)))*...
                 (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                 (edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))'); 
             gradm6b = [0, (2/d).*((X(nv+j,2)-cen(2))/d)];
             B(6,nv+j) = dot((0.5*(5/6).*(gradm6b))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm6c = [0, (2/d).*((X(2*nv+j,2)-cen(2))/d)];
             B(6,2*nv+j)= dot((0.5*(5/6).*(gradm6c))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
            
             gradm7a = [(3/d).*((X(j,1)-cen(1))/d).^2, 0];
             B(7,j) = dot((0.5*(1/6).*(gradm7a))...
                 ,((edge_length(mod_wrap(j-1,n_sides)))*...
                 (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                 (edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm7b = [(3/d).*((X(nv+j,1)-cen(1))/d).^2, 0];
             B(7,nv+j) = dot(0.5*(5/6).*(gradm7b)...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm7c = [(3/d).*((X(2*nv+j,1)-cen(1))/d).^2, 0];
             B(7,2*nv+j)= dot(0.5*(5/6).*(gradm7c)...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');

             gradm8a = [(2/d)*((X(j,1)-cen(1))/d).*((X(j,2)-cen(2))/d), ((X(j,1)-cen(1))/d).^2.*(1/d)];
             B(8,j) = dot((0.5*(1/6).*(gradm8a))...
                 ,((edge_length(mod_wrap(j-1,n_sides)))*...
                 (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                 (edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm8b = [(2/d)*((X(nv+j,1)-cen(1))/d).*((X(nv+j,2)-cen(2))/d), ((X(nv+j,1)-cen(1))/d).^2.*(1/d)];
             B(8,nv+j) = dot((0.5*(5/6).*(gradm8b))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm8c = [(2/d)*((X(2*nv+j,1)-cen(1))/d).*((X(2*nv+j,2)-cen(2))/d), ((X(2*nv+j,1)-cen(1))/d).^2.*(1/d)];
             B(8,2*nv+j)= dot((0.5*(5/6).*(gradm8c))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');

             gradm9a = [(1/d).*((X(j,2)-cen(2))/d).^2, (2/d)*((X(j,1)-cen(1))/d).*((X(j,2)-cen(2))/d)];
             B(9,j) = dot(0.5*(1/6).*(gradm9a)...
                 ,((edge_length(mod_wrap(j-1,n_sides)))*...
                 (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                 (edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm9b = [(1/d).*((X(nv+j,2)-cen(2))/d).^2, (2/d)*((X(nv+j,1)-cen(1))/d).*((X(nv+j,2)-cen(2))/d)];
             B(9,nv+j) = dot((0.5*(5/6).*(gradm9b))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm9c = [(1/d).*((X(2*nv+j,2)-cen(2))/d).^2, (2/d)*((X(2*nv+j,1)-cen(1))/d).*((X(2*nv+j,2)-cen(2))/d)];
             B(9,2*nv+j)= dot((0.5*(5/6).*(gradm9c))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             
             gradm10a = [0, (3/d).*((X(j,2)-cen(2))/d).^2];
             B(10,j) = dot((0.5*(1/6).*(gradm10a))...
                 ,((edge_length(mod_wrap(j-1,n_sides)))*...
                 (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                 (edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm10b = [0, (3/d).*((X(nv+j,2)-cen(2))/d).^2];
             B(10,nv+j) = dot((0.5*(5/6).*(gradm10b))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))');
             gradm10c = [0, (3/d).*((X(2*nv+j,2)-cen(2))/d).^2];
             B(10,2*nv+j)= dot((0.5*(5/6).*(gradm10c))...
                 ,((edge_length(mod_wrap(j,n_sides)))*...
                 (unit_nv(mod_wrap(j,n_sides),:)))'); 
end

% Internal
B(4,dof-2) = -(2/d)*(1/d)*area;
B(6,dof-2) = B(4,dof-2);
B(7,dof-1) = -(6/d)*(1/d)*area; 
B(8,dof)   = -(2/d)*(1/d)*area;
B(9,dof-1) = -(2/d)*(1/d)*area;
B(10,dof)  = -(6/d)*(1/d)*area;

%% Local D Matrix
D = zeros(dof,n_polys);
D(:,1) = 1;

% Boundary dof
for i = 1:dof-1
    Dm(2) = (X(i,1)-cen(1))/d;
    Dm(3) = (X(i,2)-cen(2))/d;
    Dm(4) = ((X(i,1)-cen(1))/d)*((X(i,1)-cen(1))/d);
    Dm(5) = ((X(i,1)-cen(1))/d)*((X(i,2)-cen(2))/d);
    Dm(6) = ((X(i,2)-cen(2))/d)*((X(i,2)-cen(2))/d);
    Dm(7) = ((X(i,1)-cen(1))/d)*((X(i,1)-cen(1))/d)*((X(i,1)-cen(1))/d);
    Dm(8) = ((X(i,1)-cen(1))/d)*((X(i,1)-cen(1))/d)*((X(i,2)-cen(2))/d);
    Dm(9) = ((X(i,1)-cen(1))/d)*((X(i,2)-cen(2))/d)*((X(i,2)-cen(2))/d);
    Dm(10)= ((X(i,2)-cen(2))/d)*((X(i,2)-cen(2))/d)*((X(i,2)-cen(2))/d);
    
    for j = 2:n_polys
        D(i,j) = Dm(j);
    end
end

% Internal dof
D(dof-2,1)=1;              D(dof-1,1)=(1/area)*z1;    D(dof,1)=(1/area)*z2;  
D(dof-2,2)=(1/area)*z1;    D(dof-1,2)=(1/area)*z3;    D(dof,2)=(1/area)*z4;  
D(dof-2,3)=(1/area)*z2;    D(dof-1,3)=(1/area)*z4;    D(dof,3)=(1/area)*z5;  
D(dof-2,4)=(1/area)*z3;    D(dof-1,4)=(1/area)*z6;    D(dof,4)=(1/area)*z7;  
D(dof-2,5)=(1/area)*z4;    D(dof-1,5)=(1/area)*z7;    D(dof,5)=(1/area)*z8;  
D(dof-2,6)=(1/area)*z5;    D(dof-1,6)=(1/area)*z8;    D(dof,6)=(1/area)*z9;  
D(dof-2,7)=(1/area)*z6;    D(dof-1,7)=(1/area)*z10;   D(dof,7)=(1/area)*z11;  
D(dof-2,8)=(1/area)*z7;    D(dof-1,8)=(1/area)*z11;   D(dof,8)=(1/area)*z12;
D(dof-2,9)=(1/area)*z8;    D(dof-1,9)=(1/area)*z12;   D(dof,9)=(1/area)*z13;  
D(dof-2,10)=(1/area)*z9;   D(dof-1,10)=(1/area)*z13;  D(dof,10)=(1/area)*z14;

%% Local H Matrix
H = zeros(n_polys,n_polys);
Co = zeros(nv,6);
loc_con = zeros(nv,6);

for i = 1:nv
    if i < nv
        loc_con(i,:) = [newX(i,:) newX(i+1,:) cen];
    else
        loc_con(i,:) = [newX(nv,:) newX(1,:) cen];
    end
end
    
for k = 1:nv
    Co(k,:) = loc_con(k,:);
    At = (1/2)*abs(Co(k,1)*(Co(k,4)-Co(k,6))+Co(k,3)*(Co(k,6)-Co(k,2))...
        + Co(k,5)*(Co(k,2)-Co(k,4)));
    
    for j = 1:length(xw(:,1))
        %transformed coordinates
        x = Co(k,1)*(1-xw(j,1)-xw(j,2))+Co(k,3)*xw(j,1)+Co(k,5)*xw(j,2);
        y = Co(k,2)*(1-xw(j,1)-xw(j,2))+Co(k,4)*xw(j,1)+Co(k,6)*xw(j,2);
        
        tm = [1; (x - cen(1))/d; (y - cen(2))/d; ...
              ((x - cen(1))/d).^2; ((x - cen(1))/d).*...
              ((y - cen(2))/d); ((y - cen(2))/d).^2; ((x - cen(1))/d).^3; ...
              ((x - cen(1))/d).^2.*((y - cen(2))/d); ((x - cen(1))/d).*((y - cen(2))/d).^2; ...
              ((y - cen(2))/d).^3];
          
        f = zeros(n_polys,n_polys);
        for r = 1:n_polys
            for s = 1:n_polys
                f(r,s) = tm(r)*tm(s);
            end
        end
        for p = 1:n_polys
            for q = 1:n_polys
                H(p,q) = H(p,q) + At*f(p,q)*xw(j,3);
            end
        end 
    end
end
% end

        
    
    
    
    
    
