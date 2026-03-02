function [X,M,lam,vz] = dextre( X, met, sigmaT, muT, useCorr)
%DEXTRE decorrelation stretch for a multispectral image
%   [X,M] = dextre(X,met,sigmaT,muT,useCorr)
%   where
%      X       - multispectral image
%      met     - method: 1-standard, 2-SVD, 3-QR-SVD, 4-MGS, [X,perc]-any of 
%                the previous method with randomized subsampling
%      sigmaT  - target variances
%      muT     - target means
%      useCorr - if 0 use covariance matrix, if 1 use correlation matrix
%   and
%      X      - processed image
%      M      - transformation matrix

%   E. Crabu, F. Pes, and G. Rodriguez
%   University of Cagliari, Italy

% Last revised September 10, 2025

if (nargin<5) || isempty(useCorr), useCorr = 0; end
if (nargin<4) || isempty(muT), muT = []; end
if (nargin<3) || isempty(sigmaT), sigmaT = []; end
if (nargin<2) || isempty(met), met = 1; end
muT = muT(:)';

if length(met)>1	% randomized method
	perc = met(2);
	met = met(1);
else
	perc = [];
end

% tolerance for zero
tol = sqrt(eps);

% vectorize image
X = im2double(X);
[r c n] = size(X);
p = r*c;
if p<=1, error('Input is not an image.'), end
A = reshape(X,p,n);

% compute mean and shifted data
mu = mean(A);
if isempty(muT), muT = mu; end
%Am = A-ones(p,1)*mu;
Am = A-mu;
vz = false(n,1);

if ~isempty(perc)	% construct subsampling
	samp = sort(randperm(p,round(p*perc)));
	%samp = sprand(p,1,perc)>0;
	Am = Am(samp,:);
	p = size(Am,1);
end

switch met
case 1		% standard algorithm
	% compute variance-covariance matrix
	V = (Am'*Am)/(p-1);
	sigma = diag(V);
	if any(abs(sigma)<tol), error('zero variance'); end
	Sigma = diag(sqrt(sigma));
	if isempty(sigmaT), SigmaT = Sigma; else, SigmaT = diag(sigmaT); end

	% spectral factorization
	if useCorr
		C = Sigma\V/Sigma;
		[Q,Lam] = eig(C);
	else
		[Q,Lam] = eig(V);
	end
	lam = sqrt(diag(Lam));
	if any(abs(lam)<tol), error('rank-deficient V-C matrix'); end

	% create linear transformation
	M = Q*diag(1./lam)*Q'*SigmaT;
	if useCorr
		M = Sigma\M;
	end

case 2		% SVD factorization
	if useCorr, error('Correlation still unimplemented.'), end

	[~,S,Q] = svd(Am,'econ');
	lam = diag(S)/sqrt(p-1); % variance correction
	if any(abs(lam)<tol), error('rank-deficient V-C matrix'); end

	% compute variances
	if isempty(sigmaT)
		SigmaT = diag(sqrt(var(Am)));
	else
		SigmaT = diag(sigmaT);
	end

	% create linear transformation
	M = Q*diag(1./lam)*Q'*SigmaT;

case 3		% Q-less QR + SVD
	if useCorr, error('Correlation still unimplemented.'), end

	R = qr(Am,'econ');
	vz = abs(diag(R))<tol;
	if any(vz), R = R(~vz,~vz); Am = Am(:,~vz); end
	[U,S,Q] = svd(R,'econ');
	lam = diag(S)/sqrt(p-1); % variance correction

	% compute variances
	if isempty(sigmaT)
		SigmaT = diag(sqrt(sum(R.^2)/(p-1)));
	else
		SigmaT = diag(sigmaT(~vz));
	end

	% create linear transformation
	M = Q*diag(1./lam)*Q'*SigmaT;

case 4		% MGS (Modified Gram-Schmidt)
	if useCorr, error('Correlation still unimplemented.'), end

	[R, vz] = mgs(Am);
	if any(vz), Am = Am(:,~vz); end
	[U,S,Q] = svd(R,'econ');
	lam = diag(S)/sqrt(p-1); % variance correction

	% compute variances
	if isempty(sigmaT)
		SigmaT = diag(sqrt(sum(R.^2)/(p-1)));
	else
		SigmaT = diag(sigmaT(~vz));
	end

	% create linear transformation
	M = Q*diag(1./lam)*Q'*SigmaT;

otherwise
	error('Unknown method.')
end

% transform data and reassemble image
%A = A*M+ones(p,1)*(muT-mu*M);
if any(vz)
	A(:,~vz) = (A(:,~vz)-mu(~vz))*M+muT(~vz);
	A(:,vz) = 0;
else
	A = (A-mu)*M+muT;
end
X = reshape(A,r,c,n);

