function [l2err,h1err] = l2h1err_p3(X,u,pde)

dof = length(X);  % degrees of freedom
nv = (dof - 3)/3; % vertices of the polygon
newX = X(1:nv,:);

[area,cen,diameter,G,B,~,H] = localGBDH_p3(X);
% G = round(G,5); B = round(B,5); H = round(H,5);
n = size(X,1);
C = H * (G\eye(size(G))) * B;
C(1,:) = 0;  C(1,3*nv+1) = area;
C(2,:) = 0;  C(2,3*nv+2) = area;
C(3,:) = 0;  C(3,3*nv+3) = area;
l2_proj = (H\eye(size(H)))*C;

xw=TriGaussPoints(8);

loc_conn = zeros(nv,6);                                     
    for i = 1:nv                                        
        if i < nv                                       
            loc_conn(i,:) = [newX(i,:), newX(i+1,:), cen];    
        else                                                     
            loc_conn(i,:) = [newX(nv,:), newX(1,:), cen];    
        end
    end
    % gauss quad on triangular elements
    Co = zeros(nv,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int_1 = 0; int_2 = 0;
    for k = 1:nv
        Co(k,:) = loc_conn(k,:);
        At = 0.5*abs(Co(k,1)*(Co(k,4)-Co(k,6))+Co(k,3)*(Co(k,6)-Co(k,2))...
            +Co(k,5)*(Co(k,2)-Co(k,4)));

        z1 = 0; z2 = 0;
        for j = 1:length(xw(:,1))
            x = Co(k,1)*(1-xw(j,1)-xw(j,2))+Co(k,3)*xw(j,1)+Co(k,5)*xw(j,2);
            y = Co(k,2)*(1-xw(j,1)-xw(j,2))+Co(k,4)*xw(j,1)+Co(k,6)*xw(j,2);
            X = [x y];
            
            u_ex=pde.exactu(X);
            dt1=pde.Dux(X);
            dt2=pde.Duy(X);

            [l2u,l2ux,l2uy]=l2proj_uh_p3(x,y,l2_proj,u,dof,cen(1),cen(2),diameter);
    
            erl2 = u_ex - l2u;
            grad_xc = dt1 - l2ux;
            grad_yc = dt2 - l2uy;
    
            z1 = z1 + erl2*erl2*xw(j,3);
            z2 = z2 + (grad_xc*grad_xc+grad_yc*grad_yc)*xw(j,3);
        end
        int_1 = int_1 + At*z1;
        int_2 = int_2 + At*z2;
    end
l2err = int_1;
h1err = int_2;
end