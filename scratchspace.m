ei = find(all(data.edges == [vv1 vv2],2) | all(data.edges == [vv2 vv1],2));
fdiffs(find(order==find(innerEdges== ei )))
    
ev = data.vertices(data.edges(:,2),:)-data.vertices(data.edges(:,1),:);
figure; hold all; axis equal;
quiver(data.vertices(data.edges(:,1),1),data.vertices(data.edges(:,1),2),ev(:,1),ev(:,2))

escore = fdiffs;
% escore = -adiffs/3;
% escore = -adiff_part1;
escore = -adiff_part2/2;
figure; hold all; % title('fdiff scores for low rows to mid rows.')
title('fdiff scores for all rows.')
maxj = 2;
for j=1:maxj
    vstart = 103 + j*100 - 100;
    vend = 199 + j*100 - 100;
    vs = vstart:vend-1;
    fds(j,numel(vs))=nan;
    for vi = 1:numel(vs)
        vv1 = vs(vi);
        vv2 = vv1+1; % horiz
%         vv2 = vv1+1+100; % diag
%         vv2 = vv1+100; % vert
        xxs(j,vi) = norm(data.vertices(vv2,:) - data.vertices(vstart,:));
        try;
            eind = find(all(data.edges == [vv1 vv2],2) | all(data.edges == [vv2 vv1],2));
            fds(j,vi) = escore(find(order==find(innerEdges==eind)));
        catch; end;
    end
%     plot(xxs(j,:), fds(j,:), '.-','color',[j j/5 0]/maxj);
    plot(fds(j,:), '.-','color',[j j/5 0]/maxj);
    pause(.1)
end; legend(split(num2str(1:maxj)))

scatter(data.vertices(vstart,1),data.vertices(vstart,2),'g','filled')
scatter(data.vertices(vend,1),data.vertices(vend,2),'r','filled')

scatter(data.vertices(vv1,1),data.vertices(vv1,2),'g','filled')
scatter(data.vertices(vv2,1),data.vertices(vv2,2),'r','filled')

figure; plot(xxs, fds, 'r.-')
figure; plot( fds, 'r.-')
title('fdiff scores from 303 to 399. second from bottom row of edges. 2 in from boundary.')