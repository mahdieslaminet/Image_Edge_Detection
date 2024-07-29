classdef Line
    %LINE 2D line class with applicative methods for images.
    properties(GetAccess = public)
        % cooridinates of the line, Matlab matrix axes
        x0
        x1
        y0
        y1
        width % odd number
        xPixels % x pixels of the whole line
        yPixels % y pixels
        xPoints % x exact points of the whole line
        yPoints % y points
    end   
    properties(Constant)
        WHITE = 1;
    end 
    properties(Dependent)
    end
    
    methods
        % Constructor
        function obj = Line(x0,y0,x1,y1,w)
            obj.x0 = x0;obj.x1 = x1;obj.y0 = y0;obj.y1 = y1;
            if(mod(w,2) == 0)
                w = w+1;
            end
            obj.width = w;
            obj = obj.calcPixels();
        end
        
        % calculates the line pixels in the image
        function obj = calcPixels(obj)
            x0 = obj.x0;x1 = obj.x1;y0 = obj.y0;y1 = obj.y1;
            dx = abs(x1-x0);
            dy = abs(y1-y0);
            L = max(dx,dy);
            index = 1;
            obj.xPixels = zeros(1,L+1);
            obj.yPixels = zeros(1,L+1);
    
            if x0 < x1 
                sx = 1;
            else
                sx = -1;
            end
    
            if y0 < y1
                sy = 1;
            else
                sy = -1;
            end
        
            err = dx-dy;

            while 1>0
                obj.xPixels(index) = x0;
                obj.yPixels(index) = y0;
                index = index+1;
        
                if x0 == x1 && y0 == y1
                    break;
                end
                
                e2 = 2*err;
                
                if e2 > -dy 
                    err = err - dy;
                    x0 = x0 + sx;
                end
                if e2 < dx
                    err = err + dx;
                    y0 = y0 + sy;
                end
            end
        end
        
        % calculates the line points, geometricaly
        function obj = calcPoints(obj)
            x0 = obj.x0;x1 = obj.x1;y0 = obj.y0;y1 = obj.y1;
            v = [x1-x0,y1-y0];
            L = norm(v); 
            v = v/L;
            obj.xPoints =zeros(1,L+1);
            obj.yPoints = zeros(1,L+1);
            
            for i=0:L
                pos = v*i;
                obj.xPoints(i+1) = x0+pos(1);
                obj.yPoints(i+1) = y0+pos(2);
            end
            
            obj.xPoints(end) = round(obj.xPoints(end));
            obj.yPoints(end) = round(obj.yPoints(end));
        end
        
        % returns a binary image of the line, MxN size. Assumes that
        % calcPoints was called before.
        function I = getLineImage(obj,M,N)
            I = zeros(M,N);
            %shift = obj.getShift();
            for i=1:obj.width
                %loc = (obj.width+1)/2-i;
                %loc = loc*shift;
                loc = [0,0];
                %if obj.inRange(M,N,obj.xPixels+loc(1),obj.yPixels+loc(2))
                    ind = sub2ind(size(I),obj.xPixels+loc(1),obj.yPixels+loc(2));
                    I(ind) = obj.WHITE;
                %end
            end
            I(obj.x0,obj.y0) = obj.WHITE;
            I(obj.x1,obj.y1) = obj.WHITE;
        end
        
        % returns a binary image of the edges of the line, 2 lines of width
        % 2 parallel the the original line, MxN size. Assumes that
        % calcPoints was called before.
        function E = getEdgeImage(obj,M,N)
            E = zeros(M,N);
            shift = obj.getShift();
            w = (obj.width-1)/2;
            edges = [-w,w];
            for i=1:length(edges)
                loc = edges(i);
                loc = loc*shift;
                ind = sub2ind(size(E),obj.xPixels(1:end)+loc(1),obj.yPixels(1:end)+loc(2));
                E(ind) = obj.WHITE;
                
            end
        end
        
        function E = getEdgeFilter(obj,M,N)
            E = zeros(M,N);
            %shift = obj.getShift();
            vec = [obj.x1-obj.x0,obj.y1-obj.y0];
            vec = [0 -1; 1 0]*vec';
            if abs(vec(1))>abs(vec(2))
                vec(2) = 0;
            else
                vec(1) = 0;
            end
            vec = vec./norm(vec);
            shift = round(vec);
            
            w = (obj.width-1)/2;
            edges = [-w:-1,1:w];
            for i=1:length(edges)
                loc = edges(i);
                loc = loc*shift;
                ind = sub2ind(size(E),obj.xPixels(2:end-1)+loc(1),obj.yPixels(2:end-1)+loc(2));
                ind2 = sub2ind(size(E),obj.xPixels(1)+loc(1),obj.yPixels(1)+loc(2));
                ind3 = sub2ind(size(E),obj.xPixels(end)+loc(1),obj.yPixels(end)+loc(2));
                E(ind) = 0.5*obj.WHITE*sign(edges(i));
                E(ind2) = 0.5*0.5*obj.WHITE*sign(edges(i));
                E(ind3) = 0.5*0.5*obj.WHITE*sign(edges(i));
            end
        end
        
        
        
        % samples the points along the line in image I, width w. Using
        % Interp2. Assumes that calcPoints() was called before.
        function S = samplePoints(obj,I,w)
            normal = obj.getNormal();
            S = zeros(2*w,length(obj.xPoints));
            index = 1;
            for i=-w:w
                if(i == 0)
                    continue;
                end
                shift = normal*i;
                S(index,:) = interp2(I',obj.xPoints+shift(1),obj.yPoints+shift(2));
                index = index+1;
            end
        end
        
        % returns the normalized normal off the line.
        function normal = getNormal(obj)
            lineDir = getLineDir(obj);
            if lineDir(2) == 0
                normal  = [0 1];
            else
                normal = [1 -lineDir(1)/lineDir(2)];
            end
            normal = normal/norm(normal);
        end
        
        % return the normalized direction of the line
        function lineDir = getLineDir(obj)
            x0 = obj.x0;x1 = obj.x1;y0 = obj.y0;y1 = obj.y1;
            lineDir = [x1-x0,y1-y0];
            lineDir = lineDir/norm(lineDir);
        end
        
        % sets the width of the line
        function obj = setWidth(obj,width)
            obj.width = width;
        end
        
        % return the shift in points in order to draw a line with width
        function shift = getShift(obj)
                normal = getNormal(obj);
                shift = abs(normal);
                if shift(1)>shift(2)
                    shift(2) = 0;
                else
                    shift(1) = 0;
                end
                shift = round(shift);
        end
        
        % return the length of the line
        function length = getLength(obj)
            v = [obj.x1-obj.x0,obj.y1-obj.y0];
            length = norm(v);
        end
    end
    
    methods (Static)
        % return true of the x,y coordinate are in Range of an m x n image.
        function res = inRange(m,n,x,y)
            res = all(x>1 & y>0 & x<=m & y<=n);
        end
    end
    
end

