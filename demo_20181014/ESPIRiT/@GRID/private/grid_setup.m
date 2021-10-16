function mtx = grid_setup(k,kernel,ks,imSize, os)

k = k(:);

sk = length(k);
sx= imSize(1); sy = imSize(2);

ks = ks;

% cconvert k-space samples to matrix indices
nx = (sx/2+1) + sx*real(k) ;
ny = (sy/2+1) + sy*imag(k) ;

% make sparse matrix
mtx = sparse(sk,sx*sy);

% loop over kernel
for lx=-(ks-1)/2:(ks-1)/2
	for ly = -(ks-1)/2:(ks-1)/2
			
			% find nearest samples
			nxt = floor(nx+lx);
			nyt = floor(ny+ly);
			
			% find index of samples inside matrix
			idk = find(nxt<=sx & nxt >=1 & nyt <=sy & nyt>=1);

			% find index of neares samples
			idx = (nyt(idk)-1)*sx + nxt(idk);

			% compute distance
			distx = nx(idk)-nxt(idk);
			disty = ny(idk)-nyt(idk);
			
			% compute weights
			wx = KERNEL(distx, kernel,ks,os);

			wy = KERNEL(disty, kernel,ks,os);

			%d = d + wy.*wx.*img(idx);
			mtx = mtx + sparse(idk,idx, wx.*wy, sk, sx*sy); 
	end
end

%function w = KERNEL(dist,kernel,ks,os)
%w = sinc(real(dist)/os);

function w = KERNEL(dist, kernel,ks,os)

x = linspace(-(ks-1)/2,(ks-1)/2,length(kernel));
w = interp1(x,[kernel(:)]',dist,'linear');
w(find(isnan(w))) = 0;




