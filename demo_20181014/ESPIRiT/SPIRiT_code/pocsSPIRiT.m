function x = pocsSPIRiT(data, GOP, nIter, x0, wavWeight, show)
%
%
% res = pocsSPIRIT(y, GOP, nIter, x0, wavWeight, show)
%
% Implementation of the Cartesian, POCS l1-SPIRiT reconstruction
%
% Input:
%		y - Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		GOP -	the SPIRiT kernel operator obtained by calibration
%		nIter -	Maximum number of iterations
%		x0 - initial guess
%		wavWeight - wavlet threshodling parameter
% 		show - >1 to display progress (slower)
%
% Outputs:
%		res - Full k-space matrix
%
% (c) Michael Lustig 2007
%



% if no l1 penalt then skip wavelet thresholding.
if wavWeight==0

	mask = (data==0);

	x = x0;

	for n=1:nIter
		tmpx =(x + GOP*x).*(mask); % Apply (G-I)x + x
		x = tmpx + data; % fix the data
		if show
			X = ifft2c(x);
			Xsqr = sqrt(sum(abs(X).^2,3));
			figure(show), imshow(Xsqr,[],'InitialMagnification',400);, drawnow
		end
	end

else
    
    % find the closest diadic size for the images
	[sx,sy,nc] = size(data);
	ssx = 2^ceil(log2(sx)); 
	ssy = 2^ceil(log2(sy));
	ss = max(ssx, ssy);
	W = Wavelet('Daubechies',4,4);
	%W = Wavelet('Haar',2,3);
	

	mask = (data==0);
	x = x0;
	x_old = x0;	
	for n=1:nIter
		x = (x + GOP*x ).*(mask) + data; % Apply (G-I)*x + x
        
        % apply wavelet thresholding
        X = ifft2c(x); % goto image domain
		X= zpad(X,ss,ss,nc); % zpad to the closest diadic 
		X = W*(X); % apply wavelet
		X = softThresh(X,wavWeight); % threshold ( joint sparsity)
		X = W'*(X); % get back the image
		X = crop(X,sx,sy,nc); % return to the original size
		xx = fft2c(X); % go back to k-space
		x = xx.*mask + data; % fix the data
		
		if show
			
			X = ifft2c(x);
			Xsqr = sqrt(sum(abs(X).^2,3));
			figure(show), imshow(Xsqr,[],'InitialMagnification',400);, drawnow
		end
	end

end

function x = softThresh(y,t)
% apply joint sparsity soft-thresholding 
absy = sqrt(sum(abs(y).^2,3));
unity = y./(repmat(absy,[1,1,size(y,3)])+eps);

res = absy-t;
res = (res + abs(res))/2;
x = unity.*repmat(res,[1,1,size(y,3)]);




