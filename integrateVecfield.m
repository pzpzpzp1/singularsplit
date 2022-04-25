function integrated = integrateVecfield(data,vec)
% vec = u;
    idx = [1:data.numEdges]';
    edgevecs = data.vertices(data.edges(idx,1),1:2) - data.vertices(data.edges(idx,2),1:2);
    edgediffs = dot(vec(data.edges2triangles(idx,1),:), edgevecs,2);
    D = sparse(1:data.numEdges,data.edges(:,1),data.edges(:,1)*0+1,data.numEdges,data.numVertices) - sparse(1:data.numEdges,data.edges(:,2),data.edges(:,2)*0+1,data.numEdges,data.numVertices);
    integrated = D\edgediffs;
end