function ndthingout=do2d(func,ndthingin)
old_size=size(ndthingin);
ndthingin=reshape(ndthingin,old_size(1),old_size(2),[]);
first_slice_out=func(squeeze(ndthingin(:,:,1)));
ndthingout=zeros([size(first_slice_out,1),size(first_slice_out,2),size(ndthingin,3)]);
ndthingout(:,:,1)=first_slice_out;
for iSlice=2:size(ndthingin,3)
% ndthingout(:,:,iSlice)=func(squeeze(ndthingin(:,:,iSlice)));
ndthingout(:,:,iSlice)=func(ndthingin(:,:,iSlice));
end
output_size=size(ndthingout);
ndthingout=reshape(ndthingout,[output_size(1:2) old_size(3:end)]);