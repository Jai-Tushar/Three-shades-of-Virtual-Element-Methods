function bdDofMatrix = subedgeRoutineGLob(order,lc4n,edgeCo)

Ne = length(lc4n);
edgeDof =  Ne*(order-1);
edgeDofMatrix = zeros(edgeDof,2);

[xGL] = lglnodes(2*order-1);

for j = 1:Ne
    xi = lc4n(edgeCo(j,1),1); yi = lc4n(edgeCo(j,1),2);
    xj = lc4n(edgeCo(j,2),1); yj = lc4n(edgeCo(j,2),2);
    pt = zeros(2+(order-1),2);   
    for i = 1:size(pt)
        t = 0.5*xGL(i) + 0.5;
        pt(i,1) = xi + t*(xj-xi);
        pt(i,2) = yi + t*(yj-yi);
    end
    pt(1,:) = [];   pt(end,:) = [];
    for k = 1:order-1
        edgeDofMatrix((k-1)*Ne+j,:) = pt(k,:);
    end
end

bdDofMatrix = [lc4n; edgeDofMatrix];