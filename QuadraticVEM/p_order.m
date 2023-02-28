%% order of convergence Elliptic
% % Elements in mesh 5,10,25,50,100 elements on each axis. 
% % 
% % Matrix has h(max), l2error, h1error columnwise

A = [
     0.141421356592613   0.020381695811775   0.772791256074176;
     0.056568542637045   0.003219672007576   0.309936577698053;
     0.028284271318523   0.000785058353824   0.155027956663150];
%
%  
%  
%  l2order = zeros(length(A)-1,1);
%  h1order = zeros(length(A)-1,1);
%  for i = 1:length(A)-1
%         l2order(i) = (log(A(i,2) ./ A(i+1,2)))/(log(A(i,1) ./ A(i+1,1)));
%         h1order(i) = (log(A(i,3) ./ A(i+1,3)))/(log(A(i,1) ./ A(i+1,1)));
%  end
%  [l2order,h1order]
%  
%  figure(1)
%  plot(log(A(:,1)),log(A(:,2)),'*-',log(A(:,1)),log((A(:,1)).^2))
%  title('l2-order');
%  figure(2)
%  plot(log(A(:,1)),log(A(:,3)),'+-',log(A(:,1)),log(A(:,1)))
%  title('h1-order')

% Elements in mesh 100,200,400,800,1600 row wise.
% Matrix has h(max), l2error, h1error columnwise
 
% A = [0.150785836829702   0.006925631476439   0.283608200472801;
%      0.112515980677038   0.003329756751618   0.200422507246046;
%      0.079032815787871   0.001626160093891   0.141192393970099;
%      0.055184813282302   0.000804681863072   0.1000872775560021;
%      0.040004561869555   0.000401917663169   0.070728960160473];
% 
%  
%  
%  l2order = zeros(length(A)-1,1);
%  h1order = zeros(length(A)-1,1);
%  for i = 1:length(A)-1
%         l2order(i) = (log(A(i,2) ./ A(i+1,2)))/(log(A(i,1) ./ A(i+1,1)));
%         h1order(i) = (log(A(i,3) ./ A(i+1,3)))/(log(A(i,1) ./ A(i+1,1)));
%  end
%  [l2order,h1order]
%  
%  figure(1)
%  plot(log(A(:,1)),log(A(:,2)),'*-',log(A(:,1)),log((A(:,1)).^2))
%  title('l2-order');
%  figure(2)
%  plot(log(A(:,1)),log(A(:,3)),'+-',log(A(:,1)),log(A(:,1)))
%  title('h1-order')


% elements 25, 100, 225, 625, 900
% A = [0.306466010163105   0.029104972909065   0.564535632290907;
%      0.149163783720733   0.006988946577233   0.282984891846768;
%      0.101999235846437   0.002936031579795   0.188187007078123;
%      0.061732471442536   0.001045028157393   0.113146047691643;
%      0.051235641337296   0.000718971145785   0.094373347998923];



 
 
 l2order = zeros(length(A)-1,1);
 h1order = zeros(length(A)-1,1);
 for i = 1:length(A)-1
        l2order(i) = (log(A(i,2) ./ A(i+1,2)))/(log(A(i,1) ./ A(i+1,1)));
        h1order(i) = (log(A(i,3) ./ A(i+1,3)))/(log(A(i,1) ./ A(i+1,1)));
 end
 [l2order,h1order]
 
 figure(1)
 plot(log(A(:,1)),log(A(:,2)),'*-',log(A(:,1)),log((A(:,1)).^3))
 title('l2-order');
 figure(2)
 plot(log(A(:,1)),log(A(:,3)),'+-',log(A(:,1)),log((A(:,1)).^2))
 title('h1-order')