function [Phi, Lambda, U, S, V, Atilde] = DMDext(Data,r)
H1 = Data(:,1:end-1);
H2 = Data(:,2:end);

[U,S,V] = svd(H1,'econ');

if r ~= 0
    S = S(1:r,1:r);
    U = U(:,1:r);
    V = V(:,1:r);
end

%Sinv = S^(-1);
%change above for large matrices;
D = diag(S);
Dinv = 1./D;
Sinv = diag(Dinv);

Atilde = U'*H2*V*Sinv;
[Vtilde Lambda] = eig(Atilde);
Phi = H2*V*Sinv*Vtilde;


