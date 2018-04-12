% Test SPIEC-EASI framework for Gaussian Data

% Parameters of the theta matrix

% Target matrix condition number
condT = 10; 

% Lower bounds (in absolute value)
theta_min=[2,-2];

% Upper bounds (in absolute value)
theta_max=[5,-5];

% Network types:(1) Scale-free graph
%               (2) Erdoes-Renyi graph
%               (3) Hub graph
%               (4) Cluster graph
%               (5) Band graph
%               (6) Block graph
%               (7) Geometric graph

% Different type of graphs
type = 4;

% Number of nodes
p = 100;

% Number of edges
e = 2*p;

% Create graph
[G,diffE] = createGraph(p,e,type);

% Generate Precision/Covariance matrix
[Theta,invTheta] = createPrecMat(G,theta_min,theta_max,condT);

n = 800;
disp(['Sample size: ',num2str(n)])

X = (chol(invTheta)'*randn(p,n));
C = cov(X');

figure;
imagesc(invTheta)
title('Population covariance matrix')
drawnow

figure;
imagesc(Theta)
title('Population precision matrix')
drawnow

%Correlation matrix
R=corr(X');

figure;
imagesc(R)
title('Sample correlation matrix')
drawnow

figure;
imagesc(C)
title('Sample covariance matrix')
drawnow

[p,n] = size(X);

% Test SPIEC-EASI

% General options
inopts.repNum = 20;
inopts.verboseMode = 0;
inopts.transType = 'none';
inopts.matType = 'corrMat';

% TREX-OR
% Specific options
inopts.method = 'trex_or';
inopts.nLams = 10; 

tic
[G_Trex,iCovMat_Trex,covMat_Trex,S_Trex,out_Trex] = spieceasi(X,inopts);
timeTREX_OR = toc

% TREX-AND
% Specific options
inopts.method = 'trex_and';
inopts.nLams = 10; 

tic
[G_TrexA,iCovMat_TrexA,covMat_TrexA,S_TrexA,out_TrexA] = spieceasi(X,inopts);
timeTREX_AND = toc

% QUIC
% Specific options
inopts.method = 'quic';
inopts.nLams = 20; 

% Optimal Prior
%P = 1*(G==0).*1+G.*0.1;
%P(1:p+1:p^2)=0;
%inopts.P = P;

tic
[G_Quic,iCovMat_Quic,covMat_Quic,S_Quic,out_Quic] = spieceasi(X,inopts);
timeQUIC = toc

% MB-OR
% Specific options
inopts.method = 'mb_or';
inopts.nLams = 20; 

tic
[G_MBor,iCovMat_MBor,covMat_MBor,S_MBor,out_MBor] = spieceasi(X,inopts);
timeMB_OR = toc

% MB-AND
% Specific options
inopts.method = 'mb_and';
inopts.nLams = 20; 

tic
[G_MBand,iCovMat_MBand,covMat_MBand,S_MBand,out_MBand] = spieceasi(X,inopts);
timeMB_AND = toc

figure;
plot(out_Quic.lambdaList./max(out_Quic.lambdaList),out_Quic.variability,'LineWidth',2)
hold on
plot(out_MBor.lambdaList./max(out_MBor.lambdaList),out_MBor.variability,'LineWidth',2)
plot(out_MBand.lambdaList./max(out_MBand.lambdaList),out_MBand.variability,'LineWidth',2)
plot(out_Trex.lambdaList./max(out_Trex.lambdaList),out_Trex.variability,'LineWidth',2)
plot(out_TrexA.lambdaList./max(out_TrexA.lambdaList),out_TrexA.variability,'LineWidth',2)
line([0,1],[out_Quic.opts.starsThresh, out_Quic.opts.starsThresh],'LineWidth',5)
grid on
box on
legend('QUIC','MB(or)','MB(and)','GTREX(or)','GTREX(and)')
ylabel('Variability')
xlabel('lambda path')
set(gca,'FontSize',20)

figure;
imagesc(G_MBor)
title(['Neighborhood selection graph (OR), Hamming: ',num2str(nnz(abs(G-G_MBor)))])

figure;
imagesc(G_MBand)
title(['Neighborhood selection graph (AND), Hamming: ',num2str(nnz(abs(G-G_MBand)))])

figure;
imagesc(G_Quic)
title(['Sparse Inverse Covariance graph, Hamming: ',num2str(nnz(abs(G-G_Quic)))])

figure;
imagesc(G_Trex)
title(['Graphical TREX (OR), Hamming: ',num2str(nnz(abs(G-G_Trex)))])

figure;
imagesc(G_TrexA)
title(['Graphical TREX (AND), Hamming: ',num2str(nnz(abs(G-G_TrexA)))])

threshICov = 5;
temp=abs(inv(C));
temp(1:p+1:p^2) = 0;

G_icov = (temp>threshICov);
G_icov(1:p+1:p^2) = 0;

figure;
imagesc(G_icov)
title(['Thresholded Inverse Covarinace Graph, Hamming: ',num2str(nnz(abs(G-G_icov)))])

figure;
imagesc(G)
title('True Graph')

%% Oracle runs without sub-sampling

% General options
inopts.repNum = 0;
inopts.verboseMode = 0;
inopts.transType = 'none';
inopts.matType = 'corrMat';
inopts.nLams = 200; 

% TREX-OR
% Specific options
inopts.method = 'trex_or';

tic
[G_TrexT,iCovMat_TrexT,covMat_TrexT,S_TrexT,out_TrexT] = spieceasi(X,inopts);
timeTREX_OR = toc
hammTREXor = hammDist(G_TrexT,G);

% TREX-AND
% Specific options
inopts.method = 'trex_and';

tic
[G_TrexAT,iCovMat_TrexAT,covMat_TrexAT,S_TrexAT,out_TrexAT] = spieceasi(X,inopts);
timeTREX_AND = toc

hammTREXand = hammDist(G_TrexAT,G);


% QUIC
% Specific options
inopts.method = 'quic';

% Optimal Prior
%P = 1*(G==0).*1+G.*0.1;
%P(1:p+1:p^2)=0;
%inopts.P = P;

tic
[G_QuicT,iCovMat_QuicT,covMat_QuicT,S_QuicT,out_QuicT] = spieceasi(X,inopts);
timeQUIC = toc

hammQUIC = hammDist(G_QuicT,G);

% MB-OR
% Specific options
inopts.method = 'mb_or';

tic
[G_MBorT,iCovMat_MBorT,covMat_MBorT,S_MBorT,out_MBorT] = spieceasi(X,inopts);
timeMB_OR = toc

hammMBor = hammDist(G_MBorT,G);

% MB-AND
% Specific options
inopts.method = 'mb_and';

tic
[G_MBandT,iCovMat_MBandT,covMat_MBandT,S_MBandT,out_MBandT] = spieceasi(X,inopts);
timeMB_AND = toc 

hammMBand = hammDist(G_MBandT,G);


figure;
plot(out_QuicT.lambdaList,hammQUIC,'LineWidth',2)
hold on
plot(out_MBorT.lambdaList,hammMBor,'LineWidth',2)
plot(out_MBandT.lambdaList,hammMBand,'LineWidth',2)
plot(out_TrexT.lambdaList,hammTREXor,'LineWidth',2)
plot(out_TrexAT.lambdaList,hammTREXand,'LineWidth',2)
grid on
box on
legend('QUIC','MB(or)','MB(and)','GTREX(or)','GTREX(and)')
ylabel('Hamming distance')
xlabel('regularization path')
set(gca,'FontSize',20)
set(gca,'YScale','log')






