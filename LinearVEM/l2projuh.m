function [l2_uh,uh_xc,uh_yc] = l2projuh(x,y,l2_proj,u,dof,xd,yd,d)

l2_uh = 0; uh_xc = 0; uh_yc = 0;
m2 = (x-xd)/d;  
m3 = (y-yd)/d;  

for i=1:dof
  %L^2 norm
  q1 = l2_proj(1,i)*1 + l2_proj(2,i)*m2 + l2_proj(3,i)*m3; 
  l2_uh = l2_uh + u(i)*q1;     
  %H1 norm
  q2 = l2_proj(2,i)*(1/d);   
  uh_xc = uh_xc + u(i)*(q2); 
  q3 = l2_proj(3,i)*(1/d); 
  uh_yc = uh_yc + u(i)*(q3);
end