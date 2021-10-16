function img = addPhaseRamp(img1,equalvalent_shift_in_kspace)
%%% add phase ramp that is equivalent to
%%% shift image in k-space
% implemented by Victor and Mengye
n=equalvalent_shift_in_kspace(1);
if length(equalvalent_shift_in_kspace)==2
    m=equalvalent_shift_in_kspace(2);
else
    m=0;
end
if ismatrix(img1)
    [L,W]=size(img1);
    img = img1.*exp(1i*repmat((-L/2 : L/2-1)'/L*2*pi*n , [1 W]));
    img = img.*exp(1i*repmat((-W/2 : W/2-1)/W*2*pi*m , [L 1]));
else if ~ismatrix(img1)
        [L,W,H]=size(img1(:,:,:));
        img=img1*0;
        for j=1:H
            img(:,:,j) = img1(:,:,j).*exp(1i*repmat((-L/2 : L/2-1)'/L*2*pi*n , [1 W]));
            img(:,:,j) = img(:,:,j).*exp(1i*repmat((-W/2 : W/2-1)/W*2*pi*m , [L 1]));
        end
    end
end