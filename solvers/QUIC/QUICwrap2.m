function [sigmaList, thetaList,  pathList, lambdaList, errors] = QUICwrap2(X, opts)
% Wrapper function for the QUIC, to make sure outputs (and most of the inputs)
% match that of Friedman's glasso matlab code, to reduce overhead/processing necessary for stars code

[n,p] = size(X);

lambdaList = sort(opts.lambdaList,'descend');

numLambda = length(lambdaList);

if strcmp(opts.matType,'corrMat')
    
    S = corr(X);
    
elseif strcmp(opts.matType,'covMat')
    
    S = cov(X);
    
end

if isfield(opts,'P')
    P = opts.P;
else
    P = 1;
end

% % is pop the covariance matrix?
% if ~issym(pop)
%     S = cov(pop);
% else
%     S = pop;
% end
%
% if ~penalDiag
%     SDs = diag(sqrt(diag(S)));
%     S   = inv(SDs) * S * inv(SDs);
% end

%Lmax  = max(lambdaList);
%lpath = lambdaList / Lmax;

if length(lambdaList)>1
    
    [thetaList sigmaList errors cputimeP iterP dGapP] = QUIC('path', S, P, lambdaList, opts.tolThreshold, opts.verboseMode, opts.maxIter);
else
    [thetaList sigmaList errors cputimeP iterP dGapP] = QUIC('default', S, lambdaList.*P, opts.tolThreshold, opts.verboseMode, opts.maxIter); 
end

mDDiag   =  repmat(diag(ones(p,1),0), [1, 1, numLambda]);
pathList = sign(abs(thetaList)) - mDDiag;


% Reverse order of estimate from small to large lambda
sigmaList = sigmaList(:,:,end:-1:1);
thetaList = thetaList(:,:,end:-1:1);
pathList = pathList(:,:,end:-1:1);
lambdaList = lambdaList(:,:,end:-1:1);
