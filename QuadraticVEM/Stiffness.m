%% Stiffness matrix k = 1
function [K,diameter] = Stiffness(X)

[~,~,diameter,G,B,D,H] = localGBDH_p2(X);

dof = length(X); 
h1_projector = (G\eye(size(G)))*B;
h1_stabilising_term = (eye(dof)- D * h1_projector)' * ...
                      (eye(dof)- D * h1_projector);
G_tilda = G; G_tilda(1,:) = 0;
                  
K = h1_projector' * G_tilda * h1_projector + h1_stabilising_term;