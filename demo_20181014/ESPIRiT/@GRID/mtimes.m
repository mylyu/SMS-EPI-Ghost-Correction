function res = mtimes(a,b)

res = times(a,b);
% 
% if a.adjoint
% 	b = b(:);
% 	res = a.GRID'*(b.*a.w(:));
% else
% 	b = b(:);
% 	res = a.GRID*b;
% end

