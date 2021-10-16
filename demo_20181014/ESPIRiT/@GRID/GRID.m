function  res = GRID(k,w,imSize, kernel,ksize,os)
%res = GRID(k,w,imSize, kernel,ksize,os)
%	2D GRIDDING



    res.adjoint = 0;
    if length(imSize) ==1
    	res.imSize = [imSize, imSize]*os;
    else
    	res.imSize = floor(imSize*os);
    end
    res.dataSize = size(k);
    res.w = w(:);
    res.k = k(:); 
    res.ksize = ksize;
    
    res.GRID = grid_setup(k, kernel,ksize,res.imSize, os);

	res = class(res,'GRID');
	
