function [l2er,h1er]=err(X,u,pde)

[~,cen,hd,G,B,~,H] = localGBDH(X);
dof = length(X);     
C = H * (G\eye(size(G))) * B;
l2_proj = (H\eye(size(H))) * C;
wt = TriGaussPoints(5); 


int_1 = 0; int_2=0;
for i = 1:dof
   if i ~= dof
       Xx = [X(i,1) X(i+1,1) cen(1)];
       Yy = [X(i,2) X(i+1,2) cen(2)];
   else 
       Xx = [X(dof,1) X(1,1) cen(1)];
       Yy = [X(dof,2) X(1,2) cen(2)];
   end
   
   z1=0; z2=0;
   for j=1:length(wt(:,1))
    xt=Xx(1)*(1-wt(j,1)-wt(j,2))+Xx(2)*wt(j,1)+Xx(3)*wt(j,2);
    yt=Yy(1)*(1-wt(j,1)-wt(j,2))+Yy(2)*wt(j,1)+Yy(3)*wt(j,2);
    
    v=[xt,yt];
    u_ex = pde.exactu(v); 
    dt1 = pde.Dux(v);
    dt2 = pde.Duy(v);
       
    [l2u,l2ux,l2uy]=l2projuh(xt,yt,l2_proj,u,dof,cen(1),cen(2),hd);
    
    erl2 = u_ex - l2u;
    grad_xc = dt1 - l2ux;
    grad_yc = dt2 - l2uy;
    
    z1 = z1 + erl2*erl2*wt(j,3);
    z2 = z2 + (grad_xc*grad_xc+grad_yc*grad_yc)*wt(j,3);
  end 
   int_1 = int_1 + (polyarea(Xx,Yy))*z1;  
   int_2 = int_2 + (polyarea(Xx,Yy))*z2;
 end
l2er = int_1;
h1er = int_2;
end