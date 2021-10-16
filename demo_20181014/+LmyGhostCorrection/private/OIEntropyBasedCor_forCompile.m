
% correct based on overall entropy
% correct multi-shot EPI by entropy based method
% step1: calculte the entropy on the phase grid, choose the minimum one as start point
% step2: interactively search the optimal phase correction factors
% step3: error detection
% reference : A fast and robust minimum entropy based non-interactive Nyquist ghost correction algorithm, Skare, ISMRM, 2006
% the same as SEntropyBasedCor.m, the difference is the phase correction is applied to each shot and entropy is calculated from reconstructed interleaved image.
% these two functions can be combined, but haven't been implemented.

function [AfterCor,PhasePara] = OIEntropyBasedCor_forCompile(xKy, Nshot,StartPoint)
if nargin < 3 || isempty(StartPoint)
    StartPoint = IGetStartPoint(xKy, [-pi/3,pi/3], [-0.1  0.1], [50 30],Nshot); % default  
end

[PhasePara] = fminsearch(@(x)OGetEntropy(xKy,x,Nshot),StartPoint);
PhasePara(1) = mod(PhasePara(1)+pi/2,pi)-pi/2;
AfterCor = OIPhaseCor(xKy,PhasePara,Nshot);

% error detection --- not robust enough 20130613
Cor = ifft(fftshift(AfterCor,2),[],2);
Corimage = fftshift(fft2(squeeze(Cor),256,256));
Corimage = sqrt(sum(abs(Corimage).^2 , 3));
Profile = sum(abs(Corimage),2);
DIM2 = size(Corimage,1);
% PhasePara
if  sum(Profile(DIM2/4+1:DIM2*3/4))< sum(Profile(1:DIM2/4))+sum(Profile(DIM2*3/4+1:end))
    if Nshot == 1
        if mod(PhasePara(1)+pi/2,pi)>=pi/2
            addphase=-pi/2;
        else
            addphase=pi/2;
        end
        PhasePara(1) = PhasePara(1)+addphase;
        for ii = 1:size(AfterCor,3)
            AfterCor(:,:,ii) = SPhaseCor(AfterCor(:,:,ii), [addphase,0]);
        end
%         display('----change Nshot 1');
    end
end
end

function Location = IGetStartPoint(xKy, ConRange, LinRange, Nstep, Nshot)
ConPhase = 0:Nstep(1);
LinPhase = 0:Nstep(2);
ConPhase = ConRange(1)+ConPhase.*(ConRange(2)-ConRange(1))/Nstep(1);
LinPhase = LinRange(1)+LinPhase.*(LinRange(2)-LinRange(1))/Nstep(2);

E=zeros(Nstep(1)+1,Nstep(2)+1);
for i = 1:Nstep(1)+1
    for j = 1:Nstep(2)+1
        E(i,j) = OGetEntropy(xKy,[ConPhase(i),LinPhase(j)],Nshot);
    end
end
[r, c]=find(E==min(min(E)));
Location = [ConPhase(r(1)) LinPhase(c(1))];
end

%calculte the entropy of a x-Ky data with certain phase correction
function entropy = OGetEntropy(xKy,PhaPara, Nshot)
if  nargin > 1
    x = PhaPara;
else
    x = [ 0 0 ];
end
% DIM = size(xKy);
AfterCor = OIPhaseCor(xKy, x, Nshot);
% image = fftshift(fft(AfterCor(end/4+1:end/4*3,end/4+1:end/4*3,:),[],1),1); %dimension of image = 3;
image = fft(AfterCor(ceil(end/4)+1:ceil(end/4)*3,ceil(end/4)+1:ceil(end/4)*3,:),[],1);
% image = fft(AfterCor(:,:,:),[],1);
%dimension of image = 3; image resolution was reduced by choosing data in central kspace.
image = sum(abs(image).^2 , 3);
%%
SquareSum = sum(image(:));
B = sqrt(image)./SquareSum;
PointwiseEntropy = B./log(B);
entropy = -sum(PointwiseEntropy(:));

end



