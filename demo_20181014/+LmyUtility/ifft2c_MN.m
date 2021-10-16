function res = ifft2c_MN(x,MROWS,NCOLS)

if nargin==1
    MROWS=size(x,1);
    NCOLS=size(x,2);
end

Previous_S = size(x);
S=Previous_S;
S(1:2)=[MROWS,NCOLS];
fctr = S(1)*S(2);

x = reshape(x,Previous_S(1),Previous_S(2),prod(Previous_S(3:end)));
new_x=zeros(MROWS,NCOLS,prod(Previous_S(3:end)));
row_start=1+floor(MROWS/2)-floor(Previous_S(1)/2);
row_range=row_start:row_start+Previous_S(1)-1;
col_start=1+floor(NCOLS/2)-floor(Previous_S(2)/2);
col_range=col_start:col_start+Previous_S(2)-1;
new_x(row_range,col_range,:)=x;

res = zeros(size(new_x));
for n=1:size(x,3)
res(:,:,n) = sqrt(fctr)*fftshift(ifft2(ifftshift(new_x(:,:,n))));
% res(:,:,n) = 1/sqrt(Previous_S(1)*Previous_S(2))*fctr*fftshift(ifft2(ifftshift(new_x(:,:,n))));% keep energy same
end


res = reshape(res,S);

