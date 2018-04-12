function betaMat = lassoSchmidt(X,y,lambdaPath)

[n,p] = size(X);
nLams = length(lambdaPath);

betaMat = zeros(p,nLams);

X = [ones(n,1) X];

funObj = @(w)SquaredError(w,X,y); % Loss function that L1 regularization is applied to
beta_init = zeros(p+1,1);
lopts.verbose = 0;

for i=1:nLams
    temp = L1General2_PSSgb(funObj,beta_init,lambdaPath(i)*[0;ones(p,1)],lopts);
    beta_init = temp;
    betaMat(:,i) = temp(2:end);
end
