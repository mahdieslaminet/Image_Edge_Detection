classdef Image < handle
    %% properties
    properties(GetAccess = public)
        I
        resI
        resIgray
        param
        noiseSigma
        leftTree
        rightTree
        minLength
        treeDataCell
        lengthDataCell
        lowLevel
        maxContrast
    end 
    
    properties(Constant)
        WIDTH = 1;
        IND_COL = 2;
    end 
    
    %% methods
    methods (Access = public)
        function obj = Image(I,param,noiseSigma)
            obj.I = im2double(I);
            obj.param = param;
            obj.noiseSigma = noiseSigma;
            obj = obj.clearRes();
        end
        
        function obj = buildTree(obj,triFlag)
            prm = getPrm(obj.param);
            obj.minLength = prm.minLength;
            [m,n]  = size(obj.I);
            N = numel(obj.I);
            level = floor(log2(N))-1;
            if triFlag
                obj.leftTree = TrianglesTree([1,1],[m,1],[m,n],prm.minLength,obj.I,obj.param,level,obj.noiseSigma);
                obj.rightTree = TrianglesTree([m,n],[1,n],[1,1],prm.minLength,obj.I,obj.param,level,obj.noiseSigma);
            else
                obj.leftTree = TrianglesTree([1,n],[1,1],[m,1],prm.minLength,obj.I,obj.param,level,obj.noiseSigma);
                obj.rightTree = TrianglesTree([m,1],[m,n],[1,n],prm.minLength,obj.I,obj.param,level,obj.noiseSigma);
            end
            
            obj.treeDataCell = cell(1,level);
            obj.lengthDataCell = cell(1,level);

            curTree = obj.leftTree;
            for i=level:-1:1
                obj.lengthDataCell{i}.length = curTree.length;
                obj.lengthDataCell{i}.lengthHypo = curTree.lengthHypo;
                if curTree.isLeaf()
                    break;
                end
                curTree = curTree.leftTree;
            end
            
            obj.lowLevel = i;
            obj.treeDataCell{i}.respMap = [obj.leftTree.respMap,obj.rightTree.respMap];
            obj.treeDataCell{i}.respPixels = [obj.leftTree.respPixels,obj.rightTree.respPixels];
            %obj = obj.updateFirstDataLevel(obj.leftTree);
            %obj = obj.updateFirstDataLevel(obj.rightTree);
            obj.maxContrast = zeros(1,2^level)-1;

            first = true;
            
            for i=i:level
                
                if ~first
                    blockSize = obj.lengthDataCell{i-1}.length^2+2*obj.lengthDataCell{i-1}.length*obj.lengthDataCell{i-1}.lengthHypo;
                    [obj.treeDataCell{i}.respMap obj.treeDataCell{i}.respPixels] = TrianglesTree.updateRespMapPlusPlus(obj.treeDataCell{i-1}.respMap,obj.param,obj.lengthDataCell{i}.length,obj.lengthDataCell{i}.lengthHypo,blockSize,i,N,obj.noiseSigma);
                    obj.maxContrast(i) = max(abs(obj.treeDataCell{i}.respMap(TrianglesTree.C,:)));
                end
                
                curT = zeros(1,2^level)-1;
                curC = abs(obj.treeDataCell{i}.respMap(TrianglesTree.C,:));
                curL = 0.5*obj.treeDataCell{i}.respMap(TrianglesTree.Ll,:)+0.5*obj.treeDataCell{i}.respMap(TrianglesTree.Lr,:);
                
                curL = round(curL);
                toRemove = (curL <= 0 | isnan(curL));
                curL(toRemove) = [];
                curC(toRemove) = [];
                
                [curC,I] = sort(curC);
                curL = curL(I);  
                
                curT(curL) = curC;
                obj.maxContrast = max(obj.maxContrast,curT);
                
                first = false;
            end
        end
        
        function obj = updateFirstDataLevel(obj,curTree)
            i = obj.lowLevel;
            if(curTree.isLeaf())
                if isempty(obj.treeDataCell{i})
                    obj.treeDataCell{i}.respMap = curTree.respMap;
                    obj.treeDataCell{i}.respPixels = curTree.respPixels;
                else
                    obj.treeDataCell{i}.respMap = [obj.treeDataCell{i}.respMap curTree.respMap];
                    obj.treeDataCell{i}.respPixels = [obj.treeDataCell{i}.respPixels curTree.respPixels];
                end
            else
                obj = obj.updateFirstDataLevel(curTree.leftTree);
                obj = obj.updateFirstDataLevel(curTree.rightTree);
            end
        end
      
        function obj = detectEdgesPlusPlus(obj)
            prm = getPrm(obj.param);
            len = length(obj.treeDataCell);
            edgesIndices = zeros(2,0);
            edgesScores = zeros(1,0);
            
            for i = len:-1:obj.lowLevel+1
                map = obj.treeDataCell{i}.respMap;
                if i<len
                    leftInd = map(TrianglesTree.bestLeft,:);
                    rightInd = map(TrianglesTree.bestRight,:);
                    
                    isNonIndex = isnan(edgesIndices) | edgesIndices  == 0;
                    
                    edgesIndices(isNonIndex) = 1;
                    
                    tempLeft = leftInd(edgesIndices);
                    tempRight = rightInd(edgesIndices);
                    
                    tempLeft(isNonIndex) = 0;
                    tempRight(isNonIndex) = 0;
                    
                    if size(tempLeft,1) == 1
                        tempLeft = tempLeft';
                        tempRight = tempRight';
                    end
                    
                    edgesIndices = [tempLeft;tempRight];
                end
                
                isEdge = map(TrianglesTree.S,:)>0; 
                map(:,~isEdge) = [];
                
                if prm.doMinContrastTest
                    posContrast = map(TrianglesTree.C,:) > 0;
                    posTest = map(TrianglesTree.minC,:) >= 0.5*map(TrianglesTree.C,:);
                    negContrast = map(TrianglesTree.C,:) < 0;
                    negTest = map(TrianglesTree.maxC,:) <= 0.5*map(TrianglesTree.C,:);
                    minConstTest = (posContrast & posTest) | (negContrast & negTest);  
                    map(:,~minConstTest) = [];
                end
                
                %[m,n] = size(map);
                %fprintf('Level = %d   Edges = %d\n',i,n);
                
                [Y,I1] = sort(map(TrianglesTree.NORM,:));
                I1 = I1(end:-1:1);
                [Y,I2] = sort(map(TrianglesTree.S,I1));
                I2 = I2(end:-1:1);
                map = map(:,I1(I2));
                
                edgesNewIndices = map(TrianglesTree.bestLeft:TrianglesTree.bestRight,:);
                edgesNewScores = map(TrianglesTree.S,:);
                
                padd = size(edgesIndices,1)-size(edgesNewIndices,1);
                edgesNewIndices = [edgesNewIndices;zeros(padd,size(edgesNewIndices,2))];
                edgesIndices = [edgesIndices edgesNewIndices];
        
                edgesScores = [edgesScores edgesNewScores];
            end
            pixels = obj.treeDataCell{obj.lowLevel}.respPixels;
            
            obj.resIgray = nonMaximalSupressionPlusPlus(edgesIndices,edgesScores,pixels,size(obj.resIgray),prm.nmsThres);
            
            if prm.normalizeScore
                obj.resIgray = Image.normalize(obj.resIgray);
            end
        end

           
        function obj = detectEdgesPlus(obj,nmsFlag)
            prm = getPrm(obj.param);
            len = length(obj.treeDataCell);
            indices = [];
            scores = [];
            norms = [];
            priorities = [];
            n0 = 0;
            lastPriority = 0;
            for i = len:-1:obj.lowLevel
                map = obj.treeDataCell{i}.respMap;
                map(TrianglesTree.Fl,:) = inf;
                map(TrianglesTree.NORM,indices) = norms;    
                map(TrianglesTree.Fl,indices) = priorities;
                map(TrianglesTree.S,indices) = scores;
                
                if i == obj.lowLevel
                    map(TrianglesTree.S,:) = -1;
                    map(TrianglesTree.S,indices) = scores;    
                    break;
                end
                
                isEdge = map(TrianglesTree.S,:)>0; 
                validEdges = false(size(isEdge));
                validEdges(indices) = true;
                map(:,~isEdge) = [];
                validEdges(~isEdge) = [];
                
                if prm.doMinContrastTest
                    posContrast = map(TrianglesTree.C,:) > 0;
                    posTest = map(TrianglesTree.minC,:) >= 0.5*map(TrianglesTree.C,:);
                    negContrast = map(TrianglesTree.C,:) < 0;
                    negTest = map(TrianglesTree.maxC,:) <= 0.5*map(TrianglesTree.C,:);
                    minConstTest = (posContrast & posTest) | (negContrast & negTest);
                    
                    map(:,~minConstTest & ~validEdges) = [];
                end
                
                [m,n] = size(map);
                fprintf('Level = %d   Edges = %d\n',i,n-n0);
                n0 = n;
                
                %[Y,I1] = sort(map(TrianglesTree.NORM,:));
                %I1 = I1(end:-1:1);
                [Y,I2] = sort(map(TrianglesTree.S,:));
                I2 = I2(end:-1:1);
                curLength = length(I2);
                %map = map(:,I1);
                map = map(:,I2);

                indices = [map(TrianglesTree.bestLeft,:),map(TrianglesTree.bestRight,:)];
                scores = [map(TrianglesTree.S,:),map(TrianglesTree.S,:)];
                norms = [map(TrianglesTree.NORM,:),map(TrianglesTree.NORM,:)];
                curPriorities = [map(TrianglesTree.Fl,:),map(TrianglesTree.Fl,:)];
                priorities = lastPriority+1:lastPriority+curLength;
                priorities = [priorities,priorities];
                priorities = min(priorities,curPriorities);
                
                [scores,ind] = sort(scores);
                indices = indices(ind);
                norms = norms(ind);
                priorities = priorities(ind);
                
                lastPriority = curLength+lastPriority;
            end
            isEdge = map(TrianglesTree.S,:)>0;
            validLength = map(TrianglesTree.NORM,:) > prm.minContrast*2;
            if prm.eliminateShortEdges
                valid = validLength & isEdge;
            else
                valid = isEdge;
            end
            pixels = obj.treeDataCell{obj.lowLevel}.respPixels;
            map(:,~valid) = [];
            pixels(:,~valid) = [];
            
            
            %map = map(:,validLength);
            %pixels = pixels(:,validLength);
            
            %[Y,I1] = sort(map(TrianglesTree.NORM,:));
            %I1 = I1(end:-1:1);
            %[Y,I2] = sort(map(TrianglesTree.S,I1));
            %I2 = I2(end:-1:1);
            %map = map(:,I1);
            %map = map(:,I2);
            edgeValues = map(TrianglesTree.S,:);
            edgePriority = map(TrianglesTree.Fl,:);
            %pixels = pixels(:,I1);
            %pixels = pixels(:,I2);
            
            if nmsFlag
                obj.resIgray = nonMaximalSupression(pixels,edgeValues,edgePriority,size(obj.resIgray),prm.nmsThres);
                obj.resIgray = Image.normalize(obj.resIgray);
                return;
            end

            val = zeros(size(pixels));
            val = bsxfun(@plus,val',edgeValues')';
            pixels = pixels(:);
            val = val(:);
            nans = isnan(pixels) | (pixels == 0);
            pixels(nans) = [];
            val(nans) = [];
            curVal = zeros(size(obj.resIgray));
            [Y,Ind] = sort(val);
            curVal(pixels(Ind)) = Y;
            %curVal(pixels) = val;
            obj.resIgray = curVal;
            obj.resIgray = Image.normalize(obj.resIgray);
        end
        
        function obj = detectEdges(obj,nmsFlag)
            len = length(obj.treeDataCell);
            for i = len:-1:obj.lowLevel+1
                map = obj.treeDataCell{i}.respMap;
                isEdge = map(TrianglesTree.S,:)>0;
                
                posContrast = map(TrianglesTree.C,:) > 0;
                posTest = map(TrianglesTree.minC,:)>=0.5*map(TrianglesTree.C,:);
                negContrast = map(TrianglesTree.C,:) < 0;
                negTest = map(TrianglesTree.maxC,:)<=0.5*map(TrianglesTree.C,:);
            
                minConstTest = (posContrast & posTest) | (negContrast & negTest);
                
                map(:,~isEdge | ~minConstTest) = [];
                %map(:,~isEdge) = [];
                
                %[Y,I1] = sort(map(TrianglesTree.NORM,:));
                %I1 = I1(end:-1:1);
                %[Y,I2] = sort(map(TrianglesTree.S,I1));
                %I2 = I2(end:-1:1);
                %map = map(:,I1);
                %map = map(:,I2);
                
                [m,n] = size(map);
                fprintf('Level = %d   Edges = %d\n',i,n);
                
                for j=1:n
                    obj = obj.addEdge(i,map(TrianglesTree.S,j),map(TrianglesTree.bestLeft,j),map(TrianglesTree.bestRight,j),nmsFlag);
                end  
            end
            obj.resIgray = Image.normalize(obj.resIgray);
        end
        
        function obj = markLeftCorners(obj)         
            cur = obj.leftTree;
            while ~isempty(cur)
                obj.resI(cur.p1(1),cur.p1(2),obj.IND_COL) = 1;
                obj.resI(cur.p2(1),cur.p2(2),obj.IND_COL) = 1;
                obj.resI(cur.p3(1),cur.p3(2),obj.IND_COL) = 1;
                cur = cur.leftTree;
            end            
        end
        
        function obj = markRightCorners(obj)         
            cur = obj.rightTree;
            while ~isempty(cur)
                obj.resI(cur.p1(1),cur.p1(2),obj.IND_COL) = 1;
                obj.resI(cur.p2(1),cur.p2(2),obj.IND_COL) = 1;
                obj.resI(cur.p3(1),cur.p3(2),obj.IND_COL) = 1;
                cur = cur.rightTree;
            end            
        end
  
        function obj = markAllCorners(obj)
            obj = obj.markAllCornersHelp(obj.leftTree);
            obj = obj.markAllCornersHelp(obj.rightTree);
        end
        
        function obj = markAllEdges(obj)
            obj = obj.markAllEdgesHelp(obj.leftTree);
            obj = obj.markAllEdgesHelp(obj.rightTree);
            obj.resIgray = obj.resIgray./max(obj.resIgray(:));
        end    
        
        function obj = showAllEdges(obj,depth)
        end
        
        function prob = onEdgeProb(obj,x)
            prm = getPrm(obj.params);
            s = prm.s;
            p = prm.p;
            alpha = prm.alpha;
            z = 2*s/p*gamma(1/p)-alpha*(1+exp(-(alpha/s)^p));
            prob = -abs(x./s).^p-log(z);
            prob(abs(x)<=alpha) = -inf;
        end       
        function prob = onEdgeGaussProb(obj,x,sigma)
            prob = log(normpdf(x,0,sigma));
        end       
        function prob = onEdgeUniformProb(obj,x,a,b)
            a = a/2;
            b = b/2;
            inRange = (abs(x)<=b) & (abs(x) >=a);
            prob = log(inRange./(2*(b-a)));
        end      
        function prob = offEdgeNoiseProb(obj,x,len)
            prm = getPrm(obj.params);
            w = sum(abs(prm.filter)>0);
            sigmaL = obj.noiseSigma^2/(w*len);
            prob = -0.5*log(2*pi*sigmaL)-x.^2/(2*sigmaL);
        end      
        function prob = offEdgeProb(obj,x,len)
            prm = getPrm(obj.params);
            s = prm.sBack/len;
            p = prm.pBack;
            z = 2*s/p*gamma(1/p);
            prob = -abs(x./s).^p-log(z);
        end  
        function obj = clearRes(obj)
            [m,n] = size(obj.I);
            obj.resI = zeros(m,n,3);
            obj.resIgray = zeros(m,n);
            obj.resI(:,:,1) = obj.I;
            obj.resI = ntsc2rgb(obj.resI);
        end
    end
    
    methods (Access = private)
    
        function obj = addEdge(obj,level,score,leftInd,rightInd,nmsFlag)
            obj = obj.addEdgeHelp(level-1,score,leftInd,nmsFlag);
            if rightInd~=leftInd
                obj = obj.addEdgeHelp(level-1,score,rightInd,nmsFlag);
            end
        end
        
        function obj = addEdgeHelp(obj,level,score,ind,nmsFlag)
            if level == obj.lowLevel
                curVal = zeros(size(obj.resIgray));
                pixels = obj.treeDataCell{level}.respPixels(:,ind);
                pixels(isnan(pixels) | pixels == 0) = [];
                curVal(pixels) = score;
                newPixels = curVal>0;
                L = sum(newPixels(:));
                newPixels = imdilate(newPixels>0,true(3));
                oldPixels = imdilate(obj.resIgray>0,true(3));
                coor = newPixels & oldPixels;
                coor = sum(coor(:));
                if coor/L<2 || ~nmsFlag
                    obj.resIgray = max(obj.resIgray,curVal);
                end
                return;
            end
            obj.treeDataCell{level}.respMap(TrianglesTree.S,ind) = -1;
            leftInd = obj.treeDataCell{level}.respMap(TrianglesTree.bestLeft,ind);
            rightInd = obj.treeDataCell{level}.respMap(TrianglesTree.bestRight,ind);

            obj = obj.addEdgeHelp(level-1,score,leftInd,nmsFlag);
            
            if rightInd~=leftInd
                obj = obj.addEdgeHelp(level-1,score,rightInd,nmsFlag);
            end
        end
        
        function obj = markAllCornersHelp(obj,cur)
            if isempty(cur)
                return;
            end
            
            obj.resI(cur.p1(1),cur.p1(2),obj.IND_COL) = 1;
            obj.resI(cur.p2(1),cur.p2(2),obj.IND_COL) = 1;
            obj.resI(cur.p3(1),cur.p3(2),obj.IND_COL) = 1;
            
            obj = obj.markAllCornersHelp(cur.leftTree);   
            obj = obj.markAllCornersHelp(cur.rightTree); 
        end
        
        function obj = markAllEdgesHelp(obj,cur)
            if ~cur.isLeaf()
                edgePixels = cur.edgePixels;
                edgeValues = cur.edgeMap(cur.S,:); 
                
                edgeLength = 0.5*cur.edgeMap(cur.Lr,:)+0.5*cur.edgeMap(cur.Ll,:); 
                validLength = edgeLength >= obj.minLength;
                isEdge = edgeValues > 0;
                edgePixels = edgePixels(:,isEdge);
                edgeValues = edgeValues(:,isEdge);
                val = zeros(size(edgePixels));
                val = bsxfun(@plus,val',edgeValues')';
                nans = isnan(edgePixels) | edgePixels == 0;
                edgePixels(nans) = [];
                val(nans) = [];
                curVal = zeros(size(obj.resIgray));
                [Y,I] = sort(val);
                curVal(edgePixels(I)) = Y;
                obj.resIgray = max(obj.resIgray,curVal);
                obj = obj.markAllEdgesHelp(cur.leftTree);   
                obj = obj.markAllEdgesHelp(cur.rightTree); 
            end
        end
        
        function obj = showAllEdgesHelp(obj,cur,str,depth)
            figure('name',str);
            [m,n] = size(cur.edgePixels);
            for i = 1:n
                pix = cur.edgePixels(:,i);
                pix(isnan(pix) | pix == 0) = [];
                obj.resIgray(pix) = 1;
                subplot(ceil(n/9),9,i);
                imshow(obj.resIgray);
                obj = obj.clearRes();
            end

            if ~cur.isLeaf() && depth>0
                obj = obj.showAllEdgesHelp(cur.leftTree, strcat(str,'l '),depth-1);   
                obj = obj.showAllEdgesHelp(cur.rightTree, strcat(str,'r '),depth-1); 
            end
        end
    end
    
    methods(Static)
        function I = normalize(I)
            I = im2double(I);
            I = I-min(I(:));
            I = I/max(I(:));
        end  
    end
end