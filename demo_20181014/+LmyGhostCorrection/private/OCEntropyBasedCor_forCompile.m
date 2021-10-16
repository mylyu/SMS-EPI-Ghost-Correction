% correct based on overall entropy %20130613
%  funtion: correct theta3( phase discrepancy between two shots)
% 1. theta 3 was applied to EPI 2/3/4
% 2. interatively search theta 3 by entropy method.

function [AfterCor,PhasePara] = OCEntropyBasedCor_forCompile(xKy,Nshot,StartPoint)
if nargin<3||isempty(StartPoint)
StartPoint = CGetStartPoint(xKy, [-pi/3,pi/3], 50,Nshot);
end
[PhasePara, ~,exist_flag] = fminsearch(@(x)CGetEntropy(xKy,Nshot,x),StartPoint);
% % % exist_flag 1
% % % fminsearch converged to a solution x.
% % % exist_flag 0
% % % Maximum number of function evaluations or iterations was reached.
% % % exist_flag -1
% % % Algorithm was terminated by the output function.
if exist_flag~=1
    PhasePara=StartPoint;
end
AfterCor = xKy;
for ii = 2:Nshot
    AfterCor(ii:Nshot:end,:,:) = xKy(ii:Nshot:end,:,:)*exp(1i*PhasePara(ii-1));
end
% Entropy = CGetEntropy(AfterCor)

%% error detection
Cor = ifft(fftshift(AfterCor,2),[],2);
Corimage = fftshift(fft2(squeeze(Cor),256,256));
Corimage = sqrt(sum(abs(Corimage).^2 , 3));
Profile = sum(Corimage,2);
DIM2 = size(Corimage,1);
if sum(Profile(DIM2/4+1:DIM2*3/4))< sum(Profile(1:DIM2/4))+sum(Profile(DIM2*3/4+1:end))
    AfterCor = SPhaseCor(AfterCor, [pi/2,0]);
    PhasePara = PhasePara+pi/2;
end


end


function Location = CGetStartPoint(xKy, Range, Nstep,Nshot)
Location=zeros(1,Nshot-1);
ConPhase = 0:Nstep;
ConPhase = Range(1)+ConPhase.*(Range(2)-Range(1))/Nstep;
ConPhase_try=zeros(1,Nshot-1);
for iShot=2:Nshot
    E=zeros(Nstep+1,1);
    for i = 1:Nstep+1
        ConPhase_try(iShot-1)=ConPhase(i);
        E(i) = CGetEntropy(xKy,Nshot,ConPhase_try);
    end
    [r]=find(E==min(E));
    Location(iShot-1) = [ConPhase(r(1))];
end
end


function entropy = CGetEntropy(xKy,Nshot,PhaPara)
if  nargin > 2
    x = PhaPara;
else
    x = zeros(1,Nshot-1);
end
% DIM = size(xKy);
AfterCor = xKy;
for ii=2:Nshot
    AfterCor(ii:Nshot:end,:,:) = xKy(ii:Nshot:end,:,:)*exp(1i*x(ii-1));
end
% image = fft(AfterCor(end/4+1:end/4*3,end/4+1:end/4*3,:),[],1);
image = fft(AfterCor(ceil(end/4)+1:ceil(end/4)*3,ceil(end/4)+1:ceil(end/4)*3,:),[],1);
% image = fft(AfterCor,[],1);
%dimension of image = 3; image resolution was reduced by choosing data in central kspace.
image = sum(abs(image).^2 , 3);
% image = sum(abs(image).^0.2 , 3);
% image=wiener2(image);
% SquareSum = sum(sum(image));
SquareSum = sum(image(:));
B = sqrt(image)/SquareSum;
% image = fftshift(fft(AfterCor,[],1),1);
% image = abs(image);
PointwiseEntropy = B./log(B);
entropy = -sum(sum(PointwiseEntropy));


end

