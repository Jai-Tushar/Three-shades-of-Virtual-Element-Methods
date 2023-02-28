% function b = loadterm_p3(X,f)
    
%     load 'test_pentagon_node_p3.mat'   % test element
%     X = Node;
X = [0 0; 0.5 0; 0.5 0.5; 0 0.5; 0.1382 0; 0.5 0.1382; ...
     0.1382 0.5; 0 0.1382; 0.3618 0; 0.5 0.3618; 0.3618 0.5; ...
     0 0.3618; 0.25 0.25; 0.2535 0.2535; 0.2571 0.2571];
f = @(x)(1); 

    dof = length(X);
    nv  = (dof-3)/3;
    [area,cen,hd,G,B,~,H] = localGBDH_p3(X);
    C = H * (G\eye(size(G))) * B;
    C(1,:) = 0;  
    C(1,3*nv+1) = area;
    C(2,:) = 0;
    C(2,3*nv+2) = area;
    C(3,:) = 0;
    C(3,3*nv+3) = area;
    l2_p = (H\eye(size(H)))*C;
    
    wt=TriGaussPoints(6);
    
    int_1 = 0; int_2 = 0; int_3 = 0; int_4 = 0; int_5 = 0;
    int_6 = 0; int_7 = 0; int_8 = 0; int_9 = 0; int_10= 0;
    for i = 1:dof
        Xx = [X(i,1) X(mod(i,dof)+1,1) cen(1)];
        Yy = [X(i,2) X(mod(i,dof)+1,2) cen(2)];
    
        z1 = 0; z2 = 0; z3 = 0; z4 = 0; z5 = 0;
        z6 = 0; z7 = 0; z8 = 0; z9 = 0; z10= 0;
        for j = 1:length(wt(:,1))
            xn = Xx(1)*(1-wt(j,1)-wt(j,2)) + Xx(2)*wt(j,1) + Xx(3)*wt(j,2);
            yn = Yy(1)*(1-wt(j,1)-wt(j,2)) + Yy(2)*wt(j,1) + Yy(3)*wt(j,2);
        
            v = [xn,yn];
            z1 = z1 + f(v)*wt(j,3);
            z2 = z2 + f(v)*((xn - cen(1))/hd)*wt(j,3);
            z3 = z3 + f(v)*((yn - cen(2))/hd)*wt(j,3);
            z4 = z4 + f(v)*(((xn - cen(1))/hd).^2)*wt(j,3);
            z5 = z5 + f(v)*(((xn - cen(1))/hd)*((yn - cen(2))/hd))*wt(j,3);
            z6 = z6 + f(v)*(((yn - cen(2))/hd).^2)*wt(j,3);
            z7 = z7 + f(v)*(((xn - cen(1))/hd).^3)*wt(j,3);
            z8 = z8 + f(v)*((((xn - cen(1))/hd).^2)*((yn - cen(2))/hd))*wt(j,3);
            z9 = z9 + f(v)*(((xn - cen(1))/hd).*((yn - cen(2))/hd).^2)*wt(j,3);
            z10= z10+ f(v)*(((yn - cen(2))/hd).^3)*wt(j,3);
        end
        int_1 = int_1 + polyarea(Xx,Yy)*z1;
        int_2 = int_2 + polyarea(Xx,Yy)*z2;
        int_3 = int_3 + polyarea(Xx,Yy)*z3;
        int_4 = int_4 + polyarea(Xx,Yy)*z4;
        int_5 = int_5 + polyarea(Xx,Yy)*z5;
        int_6 = int_6 + polyarea(Xx,Yy)*z6;
        int_7 = int_7 + polyarea(Xx,Yy)*z7;
        int_8 = int_8 + polyarea(Xx,Yy)*z8;
        int_9 = int_9 + polyarea(Xx,Yy)*z9;
        int_10= int_10 + polyarea(Xx,Yy)*z10;
    end
    b = l2_p(1,:)*int_1 + l2_p(2,:)*int_2 + l2_p(3,:)*int_3 + l2_p(4,:)*int_4 + ...
        l2_p(5,:)*int_5 + l2_p(6,:)*int_6 + l2_p(7,:)*int_7 + l2_p(8,:)*int_8 + ...
        l2_p(9,:)*int_9 + l2_p(10,:)*int_10;    