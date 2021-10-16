function MY_montage(varargin)
varargin{1}=reshape(varargin{1},size(varargin{1},1),size(varargin{1},2),1,[]);
if numel(varargin)==1
varargin{2}='DisplayRange';
varargin{3}=[];
end
orig_state = warning('off','images:initSize:adjustingMag');
montage(varargin{:});
warning(orig_state)