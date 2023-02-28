%% k = 1
function [area,centroid,diameter,G,B,D,H] = localGBDH(X)
% %tests
% clear all
% X = [0 0; 3 0; 3 2; 3/2 4; 0 4];
% %X = [0 0; 1 0; 1 1; 0 1];
% X = [0 0; 0.1 0; 0.1 0.8];
%% area, centroid and diameter
n_sides = length(X);
n_polys= 3;

area_components = X(:,1) .* X([2:end,1],2) - X([2:end,1],1) .* X(:,2);
area = 0.5 * abs(sum(area_components)); 
centroid = sum((X + X([2:end,1],:)) .* repmat(area_components,1,2))/(6*area);
diameter = 0;                                            
 for i = 1:(n_sides-1)                                              
    for j = (i+1):n_sides                                 
        diameter = max(diameter, norm(X(i,:) - X(j,:)));    
    end                                                    
 end
%% Local G Matrix
G = zeros(n_polys,n_polys);
G(1,1) = 1;
G(1,2) = (1/n_sides)*((sum(X(:,1)) - n_sides*centroid(1))/diameter);
G(1,3) = (1/n_sides)*((sum(X(:,2)) - n_sides*centroid(2))/diameter);
for i = 2:n_polys
    for j = 2:n_polys
        if i == j
            G(i,j) = (1/diameter).^2*area;
        else
            G(i,j) = 0;
        end
    end
end
%% Local B Matrix
B = zeros(n_polys,n_sides);
B(1,:) = 1/n_sides;

mod_wrap = @(x,a)mod(x-1,a) + 1;             % utility function for wrapping around a vector
edges = zeros(n_sides,2);
edge_length = zeros(n_sides,1);
for vertex_id = 1:n_sides
    current = X(mod_wrap(vertex_id, n_sides), :);          % Calculating matrices
    next = X(mod_wrap(vertex_id + 1, n_sides), :);
    edges(vertex_id,:) = next - current;
    edge_length(vertex_id) = sqrt((edges(vertex_id,1))^2+(edges(vertex_id,2))^2);
end
normal_vector = [edges(:,2),-edges(:,1)]; 
unit_nv = (1./edge_length).*normal_vector;
for j = 1:n_sides
    for i = 2:n_polys
        if i == 2
            B(i,j) = ((1/(2*diameter))*[1,0])*((edge_length(mod_wrap(j-1,n_sides)))*...
                (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                (edge_length(mod_wrap(j,n_sides)))*...
                (unit_nv(mod_wrap(j,n_sides),:)))';
        else
            B(i,j) = ((1/(2*diameter))*[0,1])*((edge_length(mod_wrap(j-1,n_sides)))*...
                (unit_nv(mod_wrap(j-1,n_sides),:)) + ...
                (edge_length(mod_wrap(j,n_sides)))*...
                (unit_nv(mod_wrap(j,n_sides),:)))';
        end
    end
end

%% Local D Matrix
D = zeros(n_sides,n_polys);
D(:,1) = 1;
for i = 1:n_sides
    for j = 2:n_polys
        if j == 2
            D(i,j) = (X(i)-centroid(1))/diameter;
        else
            D(i,j) = (X(i,2)-centroid(2))/diameter;
        end
    end
end
%% Local H Matrix
H = zeros(n_polys,n_polys);
% Triangulation of the polygon
loc_con = zeros(n_sides,6);
for i = 1:length(X)
    if i < length(X)
        loc_con(i,:) = [X(i,:) X(i+1,:) centroid];
    else
        loc_con(i,:) = [X(end,:) X(1,:) centroid];
    end
end

xw=TriGaussPoints(5);    
    
Co = zeros(n_sides,6);
for k = 1:n_sides
    Co(k,:) = loc_con(k,:);
    At = (1/2)*abs(Co(k,1)*(Co(k,4)-Co(k,6))+Co(k,3)*(Co(k,6)-Co(k,2))...
        + Co(k,5)*(Co(k,2)-Co(k,4)));
    
    for j = 1:length(xw(:,1))
        %transformed coordinates
        x = Co(k,1)*(1-xw(j,1)-xw(j,2))+Co(k,3)*xw(j,1)+Co(k,5)*xw(j,2);
        y = Co(k,2)*(1-xw(j,1)-xw(j,2))+Co(k,4)*xw(j,1)+Co(k,6)*xw(j,2);
        
        tm = [1; (x - centroid(1))/diameter; (y - centroid(2))/diameter];
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


        
    
    
    
    
    
