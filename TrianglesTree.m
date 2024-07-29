classdef TrianglesTree < handle
    % TRIANGLESTREE Data Structure of image triangles tree tiling
    % |p2-p1| = |p2-p3| = |p3-p1|/sqrt(2)
    properties
        
        % triangle vertices 
        p1
        p2
        p3
        
        % 2 sub triangle
        leftTree
        rightTree
        
        % responses data, S rows, number of edges indices columuns
        respMap
        edgeMap
        respPixels
        edgePixels
        
        % length of a triangle edge in pixels
        length
        
        %length of a triangle hypotenus in pixels
        lengthHypo
        
        % the edges vectors
        edgeVecCell
        edgeVecStartCell
        
        % parameter index
        param
        
        % level index
        level
        
        % number of pixels in the original image
        N
        
        % noise standart deviation
        sigma
    end
    
    properties(Constant)
        Fr = 1;
        Lr = 2;
        Fl = 3;
        Ll = 4;
        C = 5;
        minC = 6;
        maxC = 7;
        NORM = 8;
        bestLeft = 9;
        bestRight = 10;
        thetaBegin = 11;
        thetaEnd = 12;
        turnSqr = 13;
        turnNum = 14;
        S = 15;
    end 
    
    methods (Access = public)
    
        function obj = TrianglesTree(p1,p2,p3,minLength,I,param,level,sigma) 
            
            obj.p1 = p1;
            obj.p2 = p2;
            obj.p3 = p3;
            obj.length = max(abs(p2(1)-p1(1)),abs(p2(2)-p1(2)))+1;
            obj.lengthHypo = max(abs(p3(1)-p1(1)),abs(p3(2)-p1(2)))+1;
            obj.respMap = zeros(obj.S,obj.getNumOfEdges());
            %obj.respPixels = cell(1,obj.getNumOfEdges());
            %obj.respPixels = zeros(2^(level-3)*minLength,obj.getNumOfEdges());
            obj.param = param;
            obj.level = level;
            obj.N = numel(I);
            obj.sigma = sigma;
            obj = obj.updateEdgeVectors();            
            
            if(obj.length <= minLength && (p1(1) == p2(1) || p1(2) == p2(2)))
                    obj.leftTree = [];
                    obj.rightTree = [];
                    obj = obj.updateRespMapFirstLevel(I);
            else
                midPoint = round((p1+p3)/2);
                obj.leftTree = TrianglesTree(p2,midPoint,p1,minLength,I,param,level-1,sigma);
                obj.rightTree = TrianglesTree(p3,midPoint,p2,minLength,I,param,level-1,sigma);
                
                fname = sprintf('Files/level%d_%d.mat',obj.length,obj.lengthHypo);      
                obj.checkIfExistsMatFile(fname);
                
                obj.respMap = [obj.leftTree.respMap, obj.rightTree.respMap];
                obj.respPixels = [obj.leftTree.respPixels, obj.rightTree.respPixels];
                %obj = obj.updateRespMap(false,fname);
                %obj = obj.updateRespMapPlus(fname);
            end
            
            %obj = obj.markEdges();
        end
        
        function res = isLeaf(obj)
            res  = isempty(obj.leftTree) && isempty(obj.rightTree);
        end
        
        function n = getNumOfEdges(obj)
            n = (obj.length)^2+2*obj.length*obj.lengthHypo;
        end
        
    end
    
    methods (Access = public)

        function checkIfExistsMatFile(obj,fname)
            if ~exist(fname, 'file')
                obj.createLevelMatFile(fname,true);
            end
        end
        
        
        function createLevelMatFile(obj,fname,mex) 
            if mex
                res = createLevelMatFileFunc_mex( obj.length,obj.leftTree.length,obj.rightTree.length,obj.lengthHypo,obj.leftTree.lengthHypo,obj.rightTree.lengthHypo,obj.getNumOfEdges() );
                save(fname,'res');
                return
            end
            curSize = obj.leftTree.length*2;
            res = zeros(curSize,obj.getNumOfEdges()); 
            
            % first edges pair
            edge1 = 0;
            edge2 = 1;
            sign1 = -1;
            el1 = 0;
            el2 = 2;
            sign2 = 1;
            er1 = 1;
            er2 = 2;
            len1 = obj.length;
            len2 = obj.length;
            subLenLeft = obj.leftTree.length-1;
            subLenRight = obj.rightTree.length-1;
            
            for v1 = 0:len1-1
                for v2 = 0:len2-1
                    vl = v1;
                    vr = v2;
                    ind = obj.getIndexFromEdge(edge1,v1,edge2,v2,obj.length,obj.lengthHypo);
                    arr = zeros(curSize,1);
                    for vmid = 0:subLenLeft
                        indl = obj.getIndexFromEdge(el1,subLenLeft-vmid,el2,vl,subLenLeft+1,obj.leftTree.lengthHypo);
                        indr = obj.getIndexFromEdge(er1,vmid,er2,vr,subLenRight+1,obj.rightTree.lengthHypo);
                        %res{ind} = [res{ind} sign1*indl sign2*indr];
                        arr(vmid*2+1:vmid*2+2,1) = [sign1*indl ;sign2*indr];
                    end
                    res(:,ind) = arr;
                end
            end
            
            % second edges pair
            edge1 = 0;
            edge2 = 2;
            sign1 = -1;
            el1 = 0;
            el2 = 2;
            sign2 = -1;
            er1 = 0;
            er2 = 1;
            len1 = obj.length;
            len2 = obj.rightTree.length;
            gap = floor(obj.lengthHypo/2);
            subLenLeft = obj.leftTree.length-1;
            subLenRight = obj.rightTree.length-1;
            
            for v1 = 0:len1-1
                for v2 = 0:len2-1
                    ind = obj.getIndexFromEdge(edge1,v1,edge2,v2,obj.length,obj.lengthHypo);
                    arr = zeros(curSize,1);
                    vl = v1;
                    vr = v2;                   
                    for vmid = 0:subLenLeft
                        indl = obj.getIndexFromEdge(el1,vmid,el2,vl,subLenLeft+1,obj.leftTree.lengthHypo);
                        indr = obj.getIndexFromEdge(er1,vr,er2,subLenRight-vmid,subLenLeft+1,obj.rightTree.lengthHypo);
                        %res{ind} = [res{ind} sign1*indl sign2*indr];
                        arr(vmid*2+1:vmid*2+2,1) = [sign1*indl; sign2*indr];
                    end
                    res(:,ind) = arr;
                end
                
                for v2 = len2:obj.lengthHypo-1
                    ind = obj.getIndexFromEdge(edge1,v1,edge2,v2,obj.length,obj.lengthHypo);
                    indl = obj.getIndexFromEdge(1,v2-gap,2,v1,subLenLeft+1,obj.leftTree.lengthHypo);
                    %res{ind} = [res{ind} sign1*indl nan];
                    res(1:2,ind) = [sign1*indl; nan];
                end
            end
            
            % third edges pair
            edge1 = 1;
            edge2 = 2;
            sign1 = 1;
            el1 = 0;
            el2 = 1;
            sign2 = -1;
            er1 = 1;
            er2 = 2;
            len1 = obj.length;
            len2 = obj.leftTree.length;
            gap = floor(obj.lengthHypo/2);
            subLenLeft = obj.leftTree.length-1;
            subLenRight = obj.rightTree.length-1;
            
            for v1 = 0:len1-1
                for v2 = gap:len2-1+gap
                    ind = obj.getIndexFromEdge(edge1,v1,edge2,v2,obj.length,obj.lengthHypo);
                    arr = zeros(curSize,1);
                    vr = v1;
                    vl = v2-gap;                   
                    for vmid = 0:subLenLeft
                        indl = obj.getIndexFromEdge(el1,vmid,el2,vl,subLenLeft+1,obj.leftTree.lengthHypo);
                        indr = obj.getIndexFromEdge(er1,subLenRight-vmid,er2,vr,subLenRight+1,obj.rightTree.lengthHypo);
                        %res{ind} = [res{ind} sign1*indl sign2*indr];
                        arr(vmid*2+1:vmid*2+2,1) = [sign1*indl ; sign2*indr];
                    end
                    res(:,ind) = arr;
                end
                
                for v2 = 0:gap-1
                    ind = obj.getIndexFromEdge(edge1,v1,edge2,v2,obj.length,obj.lengthHypo);
                    indr = obj.getIndexFromEdge(0,v2,2,v1,subLenRight+1,obj.rightTree.lengthHypo);
                    %res{ind} = [res{ind} nan sign2*indr];
                    res(1:2,ind) = [ nan; sign2*indr];
                end
            end
            
            % save indices data
            save(fname,'res');
        end
        
        function obj = updateRespMapFirstLevel(obj,I)
            prm = getPrm(obj.param);
            w = prm.filterWidth;
            m = mean(I(:));
            Ior = I;
            str = 'symmetric';
            %str = 'replicate';
            I = padarray(padarray(I',w,str)',w,str);
            %I = padarray(I,[w w],'symmetric');
            
            %if prm.padWithNan
            %    I(1:w,:) = nan;
            %    I(end-w+1:end,:) = nan;
            %    I(:,1:w) = nan;
            %    I(:,end-w+1:end) = nan;
            %end
            
            triangleOrientation = obj.getTriangleOrientation();
            load(sprintf('Files/l%d.mat',triangleOrientation));
            minCoords = min(min(obj.p1,obj.p2),obj.p3);
            maxCoords = max(max(obj.p1,obj.p2),obj.p3);
            gap = w;
            patch = I(minCoords(1)-w+gap:maxCoords(1)+w+gap,minCoords(2)-w+gap:maxCoords(2)+w+gap);
            patch = reshape(patch,numel(patch),1);
            
            maskTableRight = maskTable < 0;
            maskTableLeft = maskTable > 0;
            
            filtersRight = sum(bsxfun(@times, maskTableRight, patch));
            filtersLeft = sum(bsxfun(@times, maskTableLeft, patch));
            
            nRight = sum(maskTableRight);
            nLeft = sum(maskTableLeft);
            sqrRight = (sum(bsxfun(@times, maskTableRight, patch.^2))-sum(bsxfun(@times, maskTableRight, patch)).^2./nRight)./(nRight-1);
            sqrLeft = (sum(bsxfun(@times, maskTableLeft, patch.^2))-sum(bsxfun(@times, maskTableLeft, patch)).^2./nLeft)./(nLeft-1);
            
            
            obj.respMap(1:4,:) = [filtersRight; metaTable(1,:); filtersLeft; metaTable(2,:)];
            obj.respMap(obj.C,:) = 0.5*obj.respMap(obj.Fr,:)./obj.respMap(obj.Lr,:)-0.5*obj.respMap(obj.Fl,:)./obj.respMap(obj.Ll,:);
            obj.respMap(obj.NORM,:) = metaTable(3,:);
            obj.respMap(obj.minC,:) = obj.respMap(obj.C,:);
            obj.respMap(obj.maxC,:) = obj.respMap(obj.C,:);
            
            if prm.radian
                myPi = pi;
                metaTable(4,:) = metaTable(4,:)*myPi/180;
            else
                myPi = 180;
            end
            
            obj.respMap(obj.thetaBegin,:) = mod(metaTable(4,:),2*myPi);
            obj.respMap(obj.thetaEnd,:) = mod(metaTable(4,:),2*myPi);
            %[m,n] = size(metaTable);
            coordsx = edgeIndx;
            coordsy = edgeIndy;
            coordsx = coordsx+obj.p1(1);
            coordsy = coordsy+obj.p1(2);
            obj.respPixels = sub2ind(size(Ior),coordsx,coordsy);
            
            prm = getPrm(obj.param);
            w = prm.filterWidth;
            %L = obj.respMap(obj.NORM,:);
            L = 0.5*(obj.respMap(obj.Lr,:)+obj.respMap(obj.Ll,:))./w;
            NL = obj.searchSpace(obj.level,obj.N,L,obj.param);
            T = obj.threshold(obj.sigma,2*w,L,NL,obj.param,obj.N);
            obj.respMap(obj.S,:) = abs(obj.respMap(obj.C,:))-T;
            
            if prm.eliminateEpsilons
                cont = obj.respMap(obj.C,:);
                epsilon =  abs(cont) < obj.sigma*prm.eliminateEpsilonsFactor;
                obj.respMap(1:4,epsilon) = nan;
            end
            
            if prm.doVarTest
                alpha = prm.varTestAlpha;
                varTestRight = abs(sqrRight./obj.sigma^2-1) > sqrt(2/alpha)./sqrt(nRight-1) ;
                varTestLeft = abs(sqrLeft./obj.sigma^2-1) > sqrt(2/alpha)./sqrt(nLeft-1);
                
                obj.respMap(1:4,varTestRight & varTestLeft) = nan;
            end
        end
        
        function obj = markEdges(obj)
            prm = getPrm(obj.param);
            w = prm.filterWidth;
            L = min(obj.respMap(obj.Lr,:),obj.respMap(obj.Ll,:));
            NL = obj.searchSpace(obj.level,obj.N,L);
            T = obj.threshold(obj.sigma,2*w,L,NL,obj.N);
            obj.respMap(obj.S,:) = abs(obj.respMap(obj.C,:))-T;    
                        
            isEdge = obj.respMap(obj.S,:)>0;
            obj.edgeMap = obj.respMap(:,isEdge);
            %obj.edgePixels = obj.respPixels(:,isEdge);
            
            posContrast = obj.edgeMap(obj.C,:) > 0;
            posTest = obj.edgeMap(obj.minC,:)>=0.5*obj.edgeMap(obj.C,:);
            negContrast = obj.edgeMap(obj.C,:) < 0;
            negTest = obj.edgeMap(obj.maxC,:)<=0.5*obj.edgeMap(obj.C,:);
            
            minConstTest = (posContrast & posTest) | (negContrast & negTest);
            obj.edgeMap(:,~minConstTest) = [];
            %obj.edgePixels(:,~minConstTest) = [];
            
            [Y,I1] = sort(obj.edgeMap(obj.NORM,:));
            I1 = I1(end:-1:1);
            [Y,I2] = sort(obj.edgeMap(obj.S,I1));
            I2 = I2(end:-1:1);
            obj.edgeMap = obj.edgeMap(:,I1);
            obj.edgeMap = obj.edgeMap(:,I2);
            %obj.edgePixels = obj.edgePixels(:,I1);
            %obj.edgePixels = obj.edgePixels(:,I2);
        end
        
        function orientation  = getTriangleOrientation(obj)
            v = obj.edgeVecCell{1};
            if v(1) == 1 && v(2) == 0
                orientation = 0;
            elseif v(1) == 0 && v(2) == 1
                orientation = 1;
            elseif v(1) == -1 && v(2) == 0
                orientation  = 2;
            elseif v(1) == 0 && v(2) == -1
                orientation = 3;
            else
                orientation = -1;
            end
        end
        
        function obj = updateEdgeVectors(obj)
            obj.edgeVecCell = cell(1,3);
            obj.edgeVecStartCell = cell(1,3);
            vEdge0 = (obj.p2-obj.p1)/(obj.length-1);
            vEdge1 = (obj.p3-obj.p2)/(obj.length-1);
            vEdge2 = (obj.p1-obj.p3)/(obj.lengthHypo-1);
            obj.edgeVecCell{1} = vEdge0;
            obj.edgeVecCell{2} = vEdge1;
            obj.edgeVecCell{3} = vEdge2;
            obj.edgeVecStartCell{1} = obj.p1;
            obj.edgeVecStartCell{2} = obj.p2;
            obj.edgeVecStartCell{3} = obj.p3;
        end
        
    end
    
    methods (Static)
        function  [respMapNew respPixelsNew]= updateRespMapPlusPlus(respMap,param,length,lengthHypo,blockSize,level,N,sigma) 
            
            prm = getPrm(param);
            w = prm.filterWidth;
            fname = sprintf('Files/level%d_%d.mat',length,lengthHypo);
            
            load(fname);
            
            if prm.radian
                myPi = pi;
            else
                myPi = 180;
            end
            
            [m,n] = size(res);
            res(res == 0) = nan;
            curArr = res;            
            signArr = sign(curArr);
            absArr = abs(curArr);
            absArr(2:2:end,:) =  absArr(2:2:end,:)+blockSize;

            blocks = size(respMap,2)/blockSize/2;
            
            signArrNew = repmat(signArr,1,blocks);
            absArrNew = repmat(absArr,1,blocks);
            
            [rows,cols] = size(absArr);
            
            vec = 0:blocks-1;
            vec = ones(cols,1)*vec;
            vec = vec(:)';
            vec = ones(rows,1)*vec;
            addInd = vec*2*blockSize;
            
            absArrNew = absArrNew+addInd;
            
            resNew = absArrNew.*signArrNew; 
            
            
            [m,n] = size(resNew);
            curArr = resNew;            
            signArr = signArrNew;
            absArr = absArrNew;
            absArrNoNan = absArr;
            absArrNoNan(isnan(absArrNoNan)) = 1;
            
            %% single edges
            isSingleLeft = isnan(curArr(2,:));
            isSingleRight = isnan(curArr(1,:));
            isSingle = isSingleLeft | isSingleRight;
            isRegular = ~isSingle;
            
            bestIndLeft = absArrNoNan(1,:);
            dataSingleLeft = respMap(:,bestIndLeft);
            %pixelsSingleLeft = respPixels(:,absArrNoNan(1,:));
            
            notSingleLeft = ~isSingleLeft;
            dataSingleLeft(:,notSingleLeft) = 0;
            bestIndLeft(:,notSingleLeft) = 0;
            %pixelsSingleLeft(:,~isSingleLeft) = 0;
            
            bestIndRight = absArrNoNan(2,:);
            dataSingleRight = respMap(:,bestIndRight);
            %pixelsSingleRight = respPixels(:,absArrNoNan(2,:));
            
            notSingleRight = ~isSingleRight;
            dataSingleRight(:,notSingleRight) = 0;
            bestIndRight(:,notSingleRight) = 0;
            %pixelsSingleRight(:,~isSingleRight) = 0;
            
            dataSingle = dataSingleRight+dataSingleLeft;
            %pixelsSingle = pixelsSingleRight+pixelsSingleLeft;
            
            bestInd = bestIndLeft+bestIndRight;
            
            dataSingle = [dataSingle(TrianglesTree.Fl,:); dataSingle(TrianglesTree.Ll,:); dataSingle(TrianglesTree.Fr,:); dataSingle(TrianglesTree.Lr,:); -dataSingle(TrianglesTree.C,:); -dataSingle(TrianglesTree.maxC,:); -dataSingle(TrianglesTree.minC,:); dataSingle(TrianglesTree.NORM,:);bestInd ; bestInd; mod(myPi+dataSingle(TrianglesTree.thetaEnd,:),2*myPi); mod(myPi+dataSingle(TrianglesTree.thetaBegin,:),2*myPi); dataSingle(TrianglesTree.turnSqr,:); dataSingle(TrianglesTree.turnNum,:);dataSingle(TrianglesTree.S,:)];
            %pixelsSingle = [pixelsSingle ;zeros(size(pixelsSingle))];
            
            dataSingle(:,isRegular) = 0;
            
            % regular edge
            absLeft = absArr(1:2:end,:);
            signLeft = signArr(1:2:end,:);
            signLeft = (signLeft == 1);
            absRight = absArr(2:2:end,:);
            signRight = signArr(2:2:end,:);
            signRight = (signRight == 1);
            
            % reading data into map
            absLeft(:,isSingle) = 1;
            absRight(:,isSingle) = 1;
            fll = respMap(TrianglesTree.Fl,:);
            fllMap = fll(absLeft);
            flr = respMap(TrianglesTree.Fr,:);
            flrMap = flr(absLeft);
            %frl = respMap(TrianglesTree.Fl,:);
            frlMap = fll(absRight);
            %frr = respMap(TrianglesTree.Fr,:);
            frrMap = flr(absRight);
            fllMap(:,isSingle) = 0;
            flrMap(:,isSingle) = 0;
            frlMap(:,isSingle) = 0;
            frrMap(:,isSingle) = 0;
            
            Lll = respMap(TrianglesTree.Ll,:);
            LllMap = Lll(absLeft);
            Llr = respMap(TrianglesTree.Lr,:);
            LlrMap = Llr(absLeft);
            %Lrl = respMap(TrianglesTree.Ll,:);
            LrlMap = Lll(absRight);
            %Lrr = respMap(TrianglesTree.Lr,:);
            LrrMap = Llr(absRight);
            LllMap(:,isSingle) = 0;
            LlrMap(:,isSingle) = 0;
            LrlMap(:,isSingle) = 0;
            LrrMap(:,isSingle) = 0;
            
            Tlb = respMap(TrianglesTree.thetaBegin,:);
            TlbMap = Tlb(absLeft);
            Tle = respMap(TrianglesTree.thetaEnd,:);
            TleMap = Tle(absLeft);
            %Lrb = respMap(TrianglesTree.Ll,:);
            TrbMap = Tlb(absRight);
            %Lrr = respMap(TrianglesTree.Lr,:);
            TreMap = Tle(absRight);
            TlbMap(:,isSingle) = 0;
            TleMap(:,isSingle) = 0;
            TrbMap(:,isSingle) = 0;
            TreMap(:,isSingle) = 0;
            
            Tsqr = respMap(TrianglesTree.turnSqr,:);
            TlsqrMap = Tsqr(absLeft);
            Tnum = respMap(TrianglesTree.turnNum,:);
            TlnumMap = Tnum(absLeft);
            TrsqrMap = Tsqr(absRight);
            TrnumMap = Tnum(absRight);
            
            Tsqrold = TlsqrMap+TrsqrMap;
            Tnumold = TlnumMap+TrnumMap;
           
            Tsqrold(:,isSingle) = 0;
            Tnumold(:,isSingle) = 0;           
            
            % calculation new filters
            notSignLeft = ~signLeft;
            notSignRight = ~signRight;
            filterLeft = signLeft.*fllMap+notSignLeft.*flrMap+signRight.*frlMap+notSignRight.*frrMap;
            filterRight = signLeft.*flrMap+notSignLeft.*fllMap+signRight.*frrMap+notSignRight.*frlMap;
            
            lengthLeft = signLeft.*LllMap+notSignLeft.*LlrMap+signRight.*LrlMap+notSignRight.*LrrMap;
            lengthRight = signLeft.*LlrMap+notSignLeft.*LllMap+signRight.*LrrMap+notSignRight.*LrlMap;
            
            case1 = notSignLeft & signRight;
            case2 = notSignLeft & notSignRight;
            case3 = signLeft & notSignRight;
            
            thetaBeginAll = case1.*(myPi+TleMap)+case2.*(myPi+TleMap)+case3.*(myPi+TreMap);
            thetaSplit1 = case1.*(myPi+TlbMap)+case2.*(myPi+TlbMap)+case3.*(myPi+TrbMap);
            thetaEndAll = case1.*TreMap+case2.*(myPi+TrbMap)+case3.*TleMap;
            thetaSplit2 = case1.*TrbMap+case2.*(myPi+TreMap)+case3.*TlbMap;
            
            splitAng = mod(abs(thetaSplit1-thetaSplit2),2*myPi);
            isHighAng = splitAng>myPi;
            splitAng(isHighAng) = 2*myPi-splitAng(isHighAng);
            
            if prm.radian
                isHighAng = splitAng>(prm.highestTurnAngle*myPi/180);
            else
                isHighAng = splitAng>prm.highestTurnAngle;
            end
            %splitAng(isHighAng) = splitAng(isHighAng)-2*myPi;
            
            Tsqrnew = splitAng.^2+Tsqrold;
            Tnumnew = Tnumold+1;
            respShape = Tsqrnew./Tnumnew;
            
            if prm.radian
                sigmaShape = (prm.sigmaShape*myPi/180)^2;
                shapeWeight = prm.radianScoreWeight;
            else
                sigmaShape = prm.sigmaShape^2;
                shapeWeight = prm.shapeScoreWeight;
            end
           
            
            %% continute turns code from here
            
            edgeLength = 0.5*(lengthLeft+lengthRight)./w;
            %edgeLength = 2^(level-1);
            resp = 0.5*filterRight./lengthRight-0.5*filterLeft./lengthLeft;
            NL = TrianglesTree.searchSpace(level,N,edgeLength,param);
            T = TrianglesTree.threshold(sigma,2*w,edgeLength,NL,param,N);
            
            %Tshape = sigmaShape*(2*log(2*myPi)-log(2*pi*sigmaShape)-2*log(NL)./Tnumnew);
            Tshape = sigmaShape*(2*log(1*myPi)-log(2*pi*sigmaShape));
            contrastScore = abs(resp)-T;
            shapeScore = sqrt(Tshape)-sqrt(respShape);
            shapeScore(Tnumnew == 0) = 0;
            
            %score = shapeWeight*shapeScore+(1-shapeWeight)*contrastScore;
            score = (1-shapeWeight)*contrastScore+shapeWeight*shapeScore;
            score(isHighAng) = nan;
            [maxScore,maxInd] = max(score,[],1);
            nanScore = isnan(maxScore);
            %a(0) = 1;
            indices = sub2ind(size(absLeft),maxInd,1:n);
            
            bestContrastScore = contrastScore(indices);
            posContrastScore = bestContrastScore>0;
            
            bestThetaBegin = mod(thetaBeginAll(indices),2*myPi);
            bestThetaEnd = mod(thetaEndAll(indices),2*myPi);
            bestTsqr = Tsqrnew(indices);
            bestTnum = Tnumnew(indices);
            
            bestAbsLeft = absLeft(indices);
            bestSignLeft = signLeft(indices);
            bestAbsRight = absRight(indices);
            bestSignRight = signRight(indices);

            maxResp = resp(indices);
            maxFilterLeft = filterLeft(indices);
            maxFilterRight = filterRight(indices);
            maxLengthLeft = lengthLeft(indices);
            maxLengthRight = lengthRight(indices);

            bestMinLeft = respMap(TrianglesTree.minC,bestAbsLeft);
            bestMaxLeft = respMap(TrianglesTree.maxC,bestAbsLeft);
            bestMinRight = respMap(TrianglesTree.minC,bestAbsRight);
            bestMaxRight = respMap(TrianglesTree.maxC,bestAbsRight);
            bestNormLeft = respMap(TrianglesTree.NORM,bestAbsLeft);
            bestNormRight = respMap(TrianglesTree.NORM,bestAbsRight);
            
            
            notBestSignLeft = ~bestSignLeft;
            notBestSignRight = ~bestSignRight;
            minLeft = bestSignLeft.*bestMinLeft+notBestSignLeft.*(-bestMaxLeft);
            maxLeft = bestSignLeft.*bestMaxLeft+notBestSignLeft.*(-bestMinLeft);
            minRight = bestSignRight.*bestMinRight+notBestSignRight.*(-bestMaxRight);
            maxRight = bestSignRight.*bestMaxRight+notBestSignRight.*(-bestMinRight);
            
            norm = bestNormLeft+bestNormRight;
            
            highNorm = norm>prm.minContrast;
            
            notHighNorm = ~highNorm;
            minContrast = highNorm.*min(minRight,minLeft)+notHighNorm.*maxResp;
            maxContrast = highNorm.*max(maxRight,maxLeft)+notHighNorm.*maxResp;
            
            dataRegular = [maxFilterRight;maxLengthRight;maxFilterLeft;maxLengthLeft;maxResp;minContrast;maxContrast;norm;bestAbsLeft;bestAbsRight;bestThetaBegin;bestThetaEnd;bestTsqr;bestTnum;maxScore.*posContrastScore];
            %pixelsRegular = [bestPixelsLeft;bestPixelsRight];
            dataRegular(:,nanScore & isRegular) = nan;
            dataRegular(:,isSingle) = 0;
            %pixelsRegular(:,nanScore) = 0;
            respMapNew = dataRegular+dataSingle;
            respPixelsNew = 0;%pixelsRegular+pixelsSingle;    
        end
    
        function ind = getIndexFromEdge(edgeNumber1,vertexNumber1,edgeNumber2,vertexNumber2,edgeSize,hypoSize)
            subIndEdge = vertexNumber1*edgeSize+vertexNumber2;
            subIndHypo = vertexNumber1*hypoSize+vertexNumber2;
            if edgeNumber1 == 0 && edgeNumber2 == 1
                ind = subIndEdge;
            elseif edgeNumber1 == 0 && edgeNumber2 == 2
                ind = edgeSize^2+subIndHypo;
            elseif edgeNumber1 == 1 && edgeNumber2 == 2
                ind = (edgeSize^2)+edgeSize*hypoSize+subIndHypo;
            end
            ind = ind+1;
        end
        
        % todo support hypotenus size
        function [x,y] = getXYfromIndex(edgeNumber,vertexNumber,edgeSize)
            x = 0;
            y = 0;
            if edgeNumber == 0
                x = x+vertexNumber;
            elseif edgeNumber == 1
                x = x+edgeSize-1;
                y = y+vertexNumber;
            elseif edgeNumber == 2
                x = x+edgeSize-1-vertexNumber;
                y = y+edgeSize-1-vertexNumber;
            else
                x = -1;
                y = -1;
            end
        end
        
        function NL = searchSpace(level,N,L,param)
            prm = getPrm(param);
            j = level-1;
            if prm.newSearchSpace
                %NL = exp(L/3);
                %NL = prm.searchSpaceFact*N*(2^(0.5*j));
                %NL = prm.searchSpaceFact*N*(2^(0.5*j));
            else
                NL = prm.searchSpaceFact*N;
            end
        end
        
        function T = threshold(sigma,w,L,NL,param,N)
            prm = getPrm(param);
            if prm.newThreshold
                T = sigma.*sqrt(log(NL.*w.*L)./(w.*L));
                %T = sigma.*sqrt(2*log(6*N)./(w.*L)+2/(3*w));
                %T = sigma.*sqrt(log(6*N*w*L)./(w.*L)+2/(3*w));
            else
                T = sigma.*sqrt(2.*log(NL)./(w.*L)+prm.expConst);
            end
        end
    end
    
end