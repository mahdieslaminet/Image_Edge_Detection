function [res,im] = runIm(I,addShift,param)
    close all;
   
    prm = getPrm(param);
    sigma = prm.sigma;
    
    tic;
    im = Image(I,param,sigma);
    im = im.buildTree(true);
    if prm.doNMS
        im = im.detectEdgesPlusPlus();
    end
    R = im.resIgray;
    

    
    if addShift
        Inew = zeros(size(I));
        Is = I(3:end,3:end);
        Inew(1:end-2,1:end-2) = Is;
        Inew(end-1:end,1:end) = I(1:2,1:end);
        Inew(1:end,end-1:end) = I(1:end,1:2);
        
        %Rnew = zeros(size(I));
        %Rs = R(3:end,3:end);
        %Rnew(1:end-2,1:end-2) = Rs;
        %Rnew(end-1:end,1:end) = R(1:2,1:end);
        %Rnew(1:end,end-1:end) = R(1:end,1:2);
        
        
        if max(Inew(:))>1
            Inew = Inew./255;
        end
        
        %if max(Rnew(:))>1
        %    Rnew = Rnew./255;
        %end
        
        im = Image(Inew,param,sigma);
        im = im.buildTree(true);
        
        if prm.doNMS
            im = im.detectEdgesPlusPlus();
        end
        ShiftedRnew = im.resIgray;
        Rnew = zeros(size(I));
        Rnew(3:end,3:end) = ShiftedRnew(1:end-2,1:end-2);
        
        R = max(R,Rnew);
    end
    R = R./max(R(:));
    res = R;
    toc;
    %imshow(res);
end
