function [l2err,h1err]=l2h1error(x,y,u)

X = [x y];
[area,centroid,diameter,G,B,D,H] = localGBDH(X);
n = length(x);     
C = H*(G\eye(size(G)))*B;
l2_proj = (H\eye(size(H)))*C;
x_d = centroid(1);
y_d = centroid(2);
% xw=[0.33333333333333 0.33333333333333 0.22500000000000
%     0.47014206410511 0.47014206410511 0.13239415278851
%     0.47014206410511 0.05971587178977 0.13239415278851
%     0.05971587178977 0.47014206410511 0.13239415278851
%     0.10128650732346 0.10128650732346 0.12593918054483
%     0.10128650732346 0.79742698535309 0.12593918054483
%     0.79742698535309 0.10128650732346 0.12593918054483]; 
xw=[0.33333333333333 0.33333333333333 -0.56250000000000
0.20000000000000 0.20000000000000 0.52083333333333
0.20000000000000 0.60000000000000 0.52083333333333
0.60000000000000 0.20000000000000 0.52083333333333];

q1=0;q2=0;
for i=1:n
   if i~=n
       xt=[x(i) x(i+1) x_d];
       yt=[y(i) y(i+1) y_d];
   else 
       xt=[x(n) x(1) x_d];
       yt=[y(n) y(1) y_d];
   end
   
   z1=0; z2=0;
   for j=1:length(xw(:,1))
    xn=xt(1)*(1-xw(j,1)-xw(j,2))+xt(2)*xw(j,1)+xt(3)*xw(j,2);
    yn=yt(1)*(1-xw(j,1)-xw(j,2))+yt(2)*xw(j,1)+yt(3)*xw(j,2);

    u_ex=sin(xn*pi)*sin(yn*pi);
    dt1=pi*cos(xn*pi)*sin(yn*pi);
    dt2=pi*sin(xn*pi)*cos(yn*pi);

    
    [l2u,l2ux,l2uy]=l2proj_uh(xn,yn,l2_proj,u,n,x_d,y_d,diameter);
    
    erl2 = u_ex - l2u;
    grad_xc = dt1 - l2ux;
    grad_yc = dt2 - l2uy;
    
    z1 = z1 + erl2*erl2*xw(j,3);
    z2 = z2 + (grad_xc*grad_xc+grad_yc*grad_yc)*xw(j,3);
   end
   q1 = q1 + Area(xt,yt)*z1;
   q2 = q2 + Area(xt,yt)*z2;
 end
l2err=q1;
h1err=q2;
end