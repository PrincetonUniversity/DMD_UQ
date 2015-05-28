function [ Admd, DMD_eigs, DMD_modes ] = dmd_tls(Data,r)
%%%%%%%
% Script to implement total least squares dynamic mode decomposition
% This implementation assumes a time series of data with uniform dt
%
% Modified by Scott Dawson from tlsdmd.m by Maziar Hemati
%
%%%%%%%
m = size(Data,1);
n = size(Data,2);
if nargin == 1;
    r = min(m,n);
end

%put in projection test?

[U,~,~] = svd(Data,'econ');
Ur = U(:,1:r);
% Project data onto first r POD modes
DataProjected = Ur'*Data;

Z1 = DataProjected(:,1:end-1);
Z2 = DataProjected(:,2:end);

[U,~,~] = svd([Z1;Z2],'econ');
U11 = U(1:r,1:r);
U21 = U((r+1):end,1:r);
if rank(U11)<r
    error('TLS solution does not exist')
end

%Admd = -(U12*U22^-1)'; %same as below, if full SVD is taken
Admd  = (U21*U11^-1);
[Evecs, DMD_eigs] = eig(Admd);
DMD_modes = Ur*Evecs;
end