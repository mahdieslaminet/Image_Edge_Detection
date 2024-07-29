function [ respMap,respPixels ] = updateRespMapFunc(res,leftMap,leftPixels,rightMap,rightPixels,respMap,respPixels,sigma,N,level,w,Fl,Ll,Fr,Lr,C,minC,maxC,NORM,S)
    [m,n] = size(res);
    for ind = 1:n
        curArr = res(:,ind)';
        curArr(curArr == 0) = [];
        signArr = sign(curArr);
        absArr = abs(curArr);

        if numel(curArr) == 2
            if isnan(curArr(2))
                data =  leftMap(:,absArr(1));
                pixels =  leftPixels(:,absArr(1));
            else
                data =  rightMap(:,absArr(2));
                pixels =  rightPixels(:,absArr(2));
            end

            data = [data( Fl) data( Ll) data( Fr) data( Lr) -data( C) -data( maxC) -data( minC) data( NORM) data( S)]';

             respMap(:,ind) = data;
             respPixels(:,ind) = [pixels ;zeros(size(pixels))+NaN];
        else
            absLeft = absArr(1:2:end);
            signLeft = signArr(1:2:end);
            absRight = absArr(2:2:end);
            signRight = signArr(2:2:end);

            dataLeft =  leftMap(:,absLeft);
            pixelLeft =  leftPixels(:,absLeft);
            dataRight =  rightMap(:,absRight);
            pixelRight =  rightPixels(:,absRight);

            lenGap = 0;
            if signLeft(1)>0
                filterLeft = dataLeft( Fl,:);
                lengthLeft = dataLeft( Ll,:);
                filterRight = dataLeft( Fr,:);
                lengthRight = dataLeft( Lr,:);
            else
                filterLeft = dataLeft( Fr,:);
                lengthLeft = dataLeft( Lr,:);
                filterRight = dataLeft( Fl,:);
                lengthRight = dataLeft( Ll,:);
            end

            if signRight(1)>0
                filterLeft = filterLeft+dataRight( Fl,:);
                lengthLeft = lengthLeft+dataRight( Ll,:);
                filterRight = filterRight+dataRight( Fr,:);
                lengthRight = lengthRight+dataRight( Lr,:);
            else
                filterLeft = filterLeft+dataRight( Fr,:);
                lengthLeft = lengthLeft+dataRight( Lr,:);
                filterRight = filterRight+dataRight( Fl,:);
                lengthRight = lengthRight+dataRight( Ll,:);
            end  
            lengthLeft = lengthLeft-lenGap;
            lengthRight = lengthRight-lenGap;


            minLength = min(lengthLeft,lengthRight);
            resp = 0.5*filterRight./lengthRight-0.5*filterLeft./lengthLeft;

            NL = 8*N*2^(0.25*level^2+0.25*level-0.5);
            T = 2.*sigma^2.*log(NL)./(w.*minLength);
            T = sqrt(T);
            score = abs(resp)-T;
            maxScore = max(score);
            bestInd = find(score == maxScore,1);
            

            if isempty(bestInd)
                 respMap(:,ind) = zeros( S,1)+nan;
            else
                bestInd = bestInd(1);
                maxResp = resp(bestInd);
                
                respMap(1:5,ind) = [filterRight(bestInd);lengthRight(bestInd);filterLeft(bestInd);lengthLeft(bestInd);resp(bestInd)];

                if signLeft(1) > 0
                    minLeft = dataLeft( minC,bestInd);
                    maxLeft = dataLeft( maxC,bestInd);
                else
                    maxLeft = -dataLeft( minC,bestInd);
                    minLeft = -dataLeft( maxC,bestInd);
                end
                if signRight(1) > 0
                    minRight = dataRight( minC,bestInd);
                    maxRight = dataRight( maxC,bestInd);
                else
                    maxRight = -dataRight( minC,bestInd);
                    minRight = -dataRight( maxC,bestInd);
                end

                 respMap( NORM,ind) = dataLeft( NORM,bestInd)+dataRight( NORM,bestInd);

                if  respMap( NORM,ind) > 4
                     respMap( minC,ind) = min(minLeft,minRight);
                     respMap( maxC,ind) = max(maxLeft,maxRight);
                else 
                     respMap( minC,ind) = maxResp;
                     respMap( maxC,ind) = maxResp;
                end


                 respPixels(:,ind) = [pixelLeft(:,bestInd) ;pixelRight(:,bestInd)];
            end
        end
    end       
end