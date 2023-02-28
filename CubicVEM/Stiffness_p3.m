function [K,b,hd] = Stiffness_p3(N,NT,p3NT,pde)

ndof = max(horzcat(p3NT{:})');
nel  = size(NT,1);

K = sparse(ndof,ndof);
b = zeros(ndof,1);
Di = zeros(nel,1);

for el_id = 1:nel
    v_id = NT{el_id};
    X    = N(v_id,:);
    Gl_id = p3NT{el_id,:};
    [~,~,hd,G,B,D] = localGBDHp3R1(3,X);
    h1_p = (G\eye(size(G)))*B;
    h1_st = (eye(length(Gl_id))- D * h1_p)' * (eye(length(Gl_id))- D * h1_p);
    Di(el_id,1) = Di(el_id,1) + hd;
    
    be = loadtermp3R1(X,pde.f);
    b(Gl_id,1) = b(Gl_id,1) + be;
    
    G_tilda = G; G_tilda(1,:) = 0;
    Ke = (h1_p'*G_tilda*h1_p) + h1_st;
    K(Gl_id,Gl_id) = K(Gl_id,Gl_id) + Ke;
end

hd = max(Di);
