%% Load term k = 2 
% Author - Jai Tushar, BITS-Pilani

function b_local = load_term_p2(X,pde)

dof = length(X);   % degrees of freedom
nv  = (dof - 1)/2; % vertices of the polygon
newX = X(1:nv,:);
n_sides = length(newX);

[~,centroid,diameter,G,B,~,H] = localGBDH_p2(X);
C = H * (G\eye(size(G))) * B;
l2_projector = (H\eye(size(H))) * C;

xw=TriGaussPoints(8);

% Keep for polygon but not necessary for triangular mesh%%%%%%%
loc_conn = zeros(nv,6);                                     
    for i = 1:nv                                        
        if i < nv                                       
            loc_conn(i,:) = [newX(i,:), newX(i+1,:), centroid];    
        else                                                     
            loc_conn(i,:) = [newX(nv,:), newX(1,:), centroid];    
        end
    end
    % gauss quad on triangular elements
    Co = zeros(nv,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    int_1 = 0; int_2 = 0; int_3 = 0;
    int_4 = 0; int_5 = 0; int_6 = 0;

    for k = 1:n_sides
        Co(k,:) = loc_conn(k,:);
        At = 0.5*abs(Co(k,1)*(Co(k,4)-Co(k,6))+Co(k,3)*(Co(k,6)-Co(k,2))...
            +Co(k,5)*(Co(k,2)-Co(k,4)));
        z1 = 0; z2 = 0; z3 = 0;
        z4 = 0; z5 = 0; z6 = 0;
        for j = 1:length(xw(:,1))
            x = Co(k,1)*(1-xw(j,1)-xw(j,2))+Co(k,3)*xw(j,1)+Co(k,5)*xw(j,2);
            y = Co(k,2)*(1-xw(j,1)-xw(j,2))+Co(k,4)*xw(j,1)+Co(k,6)*xw(j,2);
%             f = 2*(pi^2)*sin(pi*x)*sin(pi*y);
            XxYy = [x y];
            z1 = z1 + pde.f(XxYy)*xw(j,3);
            z2 = z2 + pde.f(XxYy)*((x-centroid(1))/diameter)*xw(j,3);
            z3 = z3 + pde.f(XxYy)*((y-centroid(2))/diameter)*xw(j,3);
            z4 = z4 + pde.f(XxYy)*(((x-centroid(1))/diameter).^2)*xw(j,3);
            z5 = z5 + pde.f(XxYy)*(((x-centroid(1))/diameter)...
                        *((y-centroid(2))/diameter))*xw(j,3);
            z6 = z6 + pde.f(XxYy)*(((y-centroid(2))/diameter).^2)*xw(j,3);         
        end
        int_1 = int_1 + At*z1;
        int_2 = int_2 + At*z2;
        int_3 = int_3 + At*z3;
        int_4 = int_4 + At*z4;
        int_5 = int_5 + At*z5;
        int_6 = int_6 + At*z6;
    end
    b_local = l2_projector(1,:)*int_1 + l2_projector(2,:)*int_2 + ...
              l2_projector(3,:)*int_3 + l2_projector(4,:)*int_4 + ...
              l2_projector(5,:)*int_5 + l2_projector(6,:)*int_6;
end