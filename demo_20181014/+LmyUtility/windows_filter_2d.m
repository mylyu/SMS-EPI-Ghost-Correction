function output=windows_filter_2d(input,window_size,window_func1,window_func2)
input_size=size(input);

if nargin<2
    window_size=input_size(1:2);
end
if nargin<3
    window_func1=@hamming;
    window_func2=window_func1;
end
if nargin<4
    window_func2=window_func1;
end

if numel(window_size)==1
    window_size(2)=window_size(1);
end
w1D_column = window_func1(window_size(1)); % Some 1D window
w1D_row = window_func2(window_size(2)); % Some 1D window
w2D = w1D_column(:) * w1D_row(:).';
freq_domain_of_input=fft2c(input);
only_center=zeros(size(freq_domain_of_input));
kx_center=floor(size(freq_domain_of_input,1)/2)+1;
ky_center=floor(size(freq_domain_of_input,2)/2)+1;
kx_lower=kx_center-floor(window_size(1)/2);
kx_upper=kx_center+ceil(window_size(1)/2)-1;
ky_lower=ky_center-floor(window_size(2)/2);
ky_upper=ky_center+ceil(window_size(2)/2)-1;
only_center(kx_lower:kx_upper,ky_lower:ky_upper,:)=repmat(w2D,[1 1 prod(input_size(3:end))] ).*freq_domain_of_input(kx_lower:kx_upper,ky_lower:ky_upper,:);
output=reshape(ifft2c(only_center),input_size);
end