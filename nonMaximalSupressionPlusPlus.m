function resI = nonMaximalSupressionPlusPlus(edgesIndices,edgesScores,pixels,imSize,coorThresh)

    resI = zeros(imSize);
    selected = false(imSize);
    
    edgeCounter = 0;
    [m,n] = size(edgesIndices);
    
    for i=1:n
        curIndices = edgesIndices(:,i);
        curIndices(curIndices == 0 | isnan(curIndices)) = [];
        curPixels = pixels(:,curIndices);
        curPixels = curPixels(:);
        curPixels(curPixels == 0 | isnan(curPixels)) = [];
        curI = false(imSize);
        curI(curPixels) = true;
        curIdialate = curI | [zeros(1,imSize(2)) ; curI(1:end-1,1:end)]|[curI(2:end,1:end);zeros(1,imSize(2))]|[curI(1:end,2:end) zeros(imSize(1),1)]|[zeros(imSize(1),1) curI(1:end,1:end-1)];
        L = sum(curI(:));
        
        coor = curIdialate & selected;
        
        if sum(coor(:))/L < coorThresh
            %figure,imshow(curI);
            edgeCounter = edgeCounter+1;
            selected = selected | curIdialate;
            resI = max(resI,curI*edgesScores(i));
        end
    end
    
    fprintf('Final Num of Edges = %d\n', edgeCounter);
end