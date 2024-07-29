function resI = nonMaximalSupression(pixels,edgeValues,edgePriority,imSize,coorThresh,nmsHolder)
    resI = nmsHolder;
    selected = nmsHolder>0;
    
    [edgePriority,I] = sort(edgePriority);
    pixels = pixels(:,I);
    edgeValues = edgeValues(1,I);
    
    lastPriority = 1;
    edgePixels = [];
    [m,n] = size(pixels);
    edgeCounter = 0;
    for i=1:n
        curPriority = edgePriority(i);
        if curPriority > lastPriority
            curEdge = zeros(imSize);
            edgePixels(isnan(edgePixels) | edgePixels == 0) = [];
            curEdge(edgePixels) = true;
            L = sum(curEdge(:));
            %curEdge = imclose(curEdge, true(5));
            curEdgeWide = imdilate(curEdge, true(3));
            
            coorIm = curEdgeWide & selected;
            coor = sum(coorIm(:))/L;
            
            if coor<coorThresh && L>5
                selected = selected | curEdgeWide;
                %figure,imshow(curEdge);
                curEdge = curEdge*edgeValues(i-1);
                resI = max(resI,curEdge);
                edgeCounter = edgeCounter+1;
            end
            edgePixels = pixels(:,i);
            lastPriority = curPriority;
        else
            edgePixels = [edgePixels;pixels(:,i)];
        end
    end
    
    fprintf('Final Num of Edges = %d\n', edgeCounter);
end