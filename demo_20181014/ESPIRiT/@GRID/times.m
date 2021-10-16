function ress = times(a,bb)
% performs only interpolation nufft

ress = zeros(a.imSize(1),a.imSize(2),size(bb,3));

for n=1:size(bb,3)
b = bb(:,:,n);
if a.adjoint
	b = b(:);
	res = full(a.GRID'*(b.*a.w(:)));
	res = reshape(res,a.imSize(1),a.imSize(2));
else
	b = b(:);
	res = full(a.GRID*b);
	res = reshape(res,a.dataSize(1),a.dataSize(2));
end
ress(:,:,n) = res;
end

