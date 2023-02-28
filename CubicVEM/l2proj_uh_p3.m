
function [l2_uh,uh_xc,uh_yc] = l2proj_uh_p3(x,y,l2_p,u,dof,xd,yd,d)

l2_uh = 0; uh_xc = 0; uh_yc = 0;
m2 = (x-xd)/d;  
m3 = (y-yd)/d;  

for i=1:dof
  %L^2 norm
  q1 = l2_p(1,i)*1 + l2_p(2,i)*m2 + l2_p(3,i)*m3 + ...
       l2_p(4,i)*m2*m2+l2_p(5,i)*m2*m3+l2_p(6,i)*m3*m3+ ...
       l2_p(7,i)*m2*m2*m2+l2_p(8,i)*m2*m2*m3+l2_p(9,i)*m2*m3*m3+...
       l2_p(10,i)*m3*m3*m3;
  l2_uh = l2_uh + u(i)*(q1);     
  %H1 norm
  q2 = l2_p(2,i)*(1/d) + l2_p(4,i)*(2/d)*(m2) + l2_p(5,i)*(1/d)*m3 + ...
       l2_p(7,i)*3*(1/d)*m2*m2 + l2_p(8,i)*2*(1/d)*m2*m3 + l2_p(9,i)*(1/d)*m3*m3;
  uh_xc = uh_xc + u(i)*(q2); 
  q3 = l2_p(3,i)*(1/d) + l2_p(5,i)*(1/d)*(m2) + l2_p(6,i)*(2/d)*m3 + ...
       l2_p(8,i)*(1/d)*m2*m2 + l2_p(9,i)*2*(1/d)*m2*m3 + l2_p(10,i)*3*(1/d)*m3*m3;
  uh_yc = uh_yc + u(i)*(q3);
end