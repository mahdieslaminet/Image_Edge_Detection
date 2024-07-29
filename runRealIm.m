function res = runRealIm(I,param)
    prm = getPrm(param);
    close all;
    addShift = prm.addShift;
    
    tic;
    [m,n] = size(I);
    R = zeros(size(I));
    block = prm.block;
    gap = prm.gap;
    
    xGrid = getGrid(block,gap,m);
    yGrid = getGrid(block,gap,n);
    
    iter = length(xGrid)*length(yGrid);
    
    iter
    
    for x0 = xGrid
        for y0 = yGrid
            curI = I(x0:x0+block-1,y0:y0+block-1);
            
            %im = Image(curI,param,sigma);
            %im = im.buildTree(true);
            %im = im.detectEdgesPlusPlus();
            %curR = im.resIgray;
            [curR,im] = runIm(curI,addShift,param);
            
            curR(1:2,:) = 0;
            curR(end-1:end,:) = 0;
            curR(:,1:2) = 0;
            curR(:,end-1:end) = 0;
            
            R(x0:x0+block-1,y0:y0+block-1) = max(R(x0:x0+block-1,y0:y0+block-1),curR);
        end
    end
    
    res = R;
end

function grid = getGrid(block,gap,max)
    grid = zeros(1,max);
    grid(1) = 1;
    for i= 2:max
         grid(i)= grid(i-1)+block-gap;
         if grid(i)+block-1>max
            grid(i) = max-block+1;
            break;
         end
    end
    grid(grid == 0) = [];
    
    if length(grid)>=2 && grid(end) == grid(end-1)
        grid(end) = [];
    end
    
end
