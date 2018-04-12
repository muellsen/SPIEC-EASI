function [wList, thetaList,  pathList, lambdaList, errors] = QUICwrap(pop, lambdaList, verbose, penalDiag, tolThreshold, maxIter)
% Wrapper function for the QUIC, to make sure outputs (and most of the inputs)
% match that of Friedman's glasso matlab code, to reduce overhead/processing necessary for stars codesy
issym = @(m) isequal(tril(m), triu(m)');
% is pop the covariance matrix?
if ~issym(pop)
    S = cov(pop);
else
    S = pop;
end

if ~penalDiag
    SDs = diag(sqrt(diag(S)));
    S   = inv(SDs) * S * inv(SDs);
end

p     = size(pop, 2);
n     = size(pop, 1);
l     = length(lambdaList);
%Lmax  = max(lambdaList);
%lpath = lambdaList / Lmax;
[thetaList wList errors cputimeP iterP dGapP] = QUIC('path', S, 1.0, lambdaList, tolThreshold, verbose, maxIter);

mDDiag   =  repmat(diag(ones(p,1),0), [1, 1, l]);
pathList = sign(abs(thetaList)) - mDDiag;
