function [R,vz,A,Rc] = mgs(A)

tol = sqrt(eps);

[p,n] = size(A);
				nz = 0;
vz = logical(zeros(n,1));
R = zeros(n);

for k = 1:n
	R(k,k) = norm(A(:,k));
	if R(k,k)>tol
		A(:,k) = A(:,k)/R(k,k);
		R(k,k+1:n) = A(:,k)'*A(:,k+1:n);
		A(:,k+1:n) = A(:,k+1:n) - A(:,k)*R(k,k+1:n);
	else
		nz = nz+1;
		vz(k) = 1;
		A(:,k) = 0;
	end
end

if nargout>3, Rc = R; end

if nz
	A = A(:,~vz);
	R = R(~vz,~vz);
end

