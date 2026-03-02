% test script for dextre.m

% WEST dataset
I = imread("westconcordaerial.png"); 
% PIMENTEL dataset
I = imread("pimentel.ppm");

I = im2double(I);
[r c n] = size(I);
sigmaT = []; muT = [];
%sigmaT = ones(n,1); muT = zeros(n,1);

perc = 0.001;

if isempty(sigmaT)
	I0 = decorrstretch(I,Mode="covariance");
else
	I0 = decorrstretch(I,Mode="covariance",TargetMean=muT,TargetSigma=sigmaT);
end

[I1,M1,lam1] = dextre(I,1,sigmaT,muT);
[I2,M2,lam2] = dextre(I,2,sigmaT,muT);
[I3,M3,lam3,vz3] = dextre(I,3,sigmaT,muT);
[I4,M4,lam4,vz4] = dextre(I,4,sigmaT,muT);
[I5,M5,lam5,vz5] = dextre(I,[3 perc],sigmaT,muT);

montage({I,I0,I2,I3,I4,I5},"Size",[2 3])
title('Original and processed images')

