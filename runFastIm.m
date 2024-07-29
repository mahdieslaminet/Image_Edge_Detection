function res = runFastIm(I,param)
    prm = getPrm(param);
    close all;
    addpath('/net/mraid11/export/data/yehonato/app2/SGE');
    
    while 1
        folderName = sprintf('tmp%s', datestr(clock));
        if ~exist(folderName,'dir')
            mkdir(folderName);
            break;
        end
    end
    
    sigma = prm.sigma;
    [m,n] = size(I);
    R = zeros(size(I));
    block = prm.block;
    gap = prm.gap;
    
    xGrid = getGrid(block,gap,m);
    yGrid = getGrid(block,gap,n);
    
    if prm.addShift
        xShift = xGrid+2;
        yShift = yGrid+2;
        
        xShift(end) = xShift(end)-4;
        yShift(end) = yShift(end)-4;
        
        xGrid = [xGrid xShift];
        yGrid = [yGrid yShift];
    end
    
    iter = length(xGrid)*length(yGrid);
    
    data = cell(1,iter);
    curInd = 0;
    
    for x0 = xGrid
        for y0 = yGrid
            curInd = curInd+1;
            curI = I(x0:x0+block-1,y0:y0+block-1);
            s.x = x0;
            s.y = y0;
            data{curInd} = s;
            
            img.I = curI;
            
            save(sprintf('%s/%d.mat',folderName,curInd),'-struct','img');
        end
    end    
    
    code = sprintf('cd /net/mraid11/export/data/yehonato/app2; runImFromMat(''%s'', index,%d)',folderName,param);
    
    code;
    
    run_parallel(code, 'index', num2cell(1:iter),'-cluster', 'mcluster01');

    for i=1:iter
        load(sprintf('%s/%dres.mat',folderName,i));
        x0 = data{i}.x;
        y0 = data{i}.y;
        R(x0:x0+block-1,y0:y0+block-1) = max(R(x0:x0+block-1,y0:y0+block-1),curR);
    end
    res = R;
    
    rmdir(folderName,'s');
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
