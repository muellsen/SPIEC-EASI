function [G,iCovMat,covMat,S,starsOut] = spieceasi(X,inopts)
% SParse InversE Covariance for Ecological ASsociation Inferene (SPIEC-EASI) code
%
% Input:    X: pxn matrix of count data where p is number of OTUs, n is number of samples
%           inopts: option structure
% Output:   G:    Interaction graph (symmetrized, nz values of iCov)
%           iCov: penalized inverse covariance matrix (only in quic
%           setting)
%           cov:  penalized covariance/correlation matrix (only in quic
%           setting)
%           S:    sample covariance/correlation matrix
%           out:  Output of the StARS selection procedure
% Author: CL Müller, Simons Foundation
%
% When using this code please cite:
%
% Kurtz, Z. D., Müller, C. L., Miraldi, E. R., Littman, D. R., Blaser, M. J., & Bonneau, R. A. (2015).
% Sparse and Compositionally Robust Inference of Microbial Ecological Networks. PLoS Computational Biology, 11(5), e1004226.
% http://doi.org/10.1371/journal.pcbi.1004226
%

% Dimension of the input data
[p,n] = size(X);

% Standard options for SPIEC-EASI
defopts.method = 'quic';               % method to infer the graph: Options: 'mb_or','mb_and','quic','trex_or','trex_and'
defopts.matType = 'corrMat';            % Matrix type: 'corrMat' (correlation matrix) or 'covMat' (covariance matrix)
defopts.transType = 'clr';              % Transformation type: 'clr' (centered log-ratio) or 'none' (no transformation)

% Optimizer settings
defopts.verboseMode = 0;                % Write out progression of the algorithm (for debugging)
defopts.maxIter = 1e3;                  % Number of iterations for optimizer

% Numerical QUIC settings
defopts.tolThreshold = 1e-7;            % Standard threshold in QUIC setting for numerical zeros

% Prior modeling: The pxp matrix models edge weight priors in the 'quic'
% setting , e.g. P = ones(p,p) is the standard settting
defopts.P = 1;

% Regularization path settings
defopts.nLams = 50;                     % number of lambda values on path
defopts.lamMinFac = 0.01;               % Standard lower bound factor in percent of the upper bound
defopts.lamMax = [];

% Options for StARS model selection
defopts.repNum = 20;                    % Number of repetitions for StARS
defopts.starsThresh = 0.1;              % StARS beta threshold
defopts.uselamLB = 1;                   % Use the heuristic lower bound after two subsamples to cut off
% the lower part of the lambda path
if n > 144;
    defopts.starsSubsampleRatio = 10 * sqrt(n)/n;
else
    defopts.starsSubsampleRatio = 0.8;
end

% Check consistency of the input options
% Merge options inopts and defopts
if nargin < 2 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

if strcmp(opts.method,'mb_or')
    
    disp('Neighborhood selection (OR rule) used');
    
elseif strcmp(opts.method,'mb_and')
    
    disp('Neighborhood selection (AND rule) used');
    
elseif strcmp(opts.method,'trex_and')
    
    disp('Graphical TREX (AND rule) used');
    
elseif strcmp(opts.method,'trex_or')
    
    disp('Graphical TREX (OR rule) used');
    
elseif strcmp(opts.method,'quic')
    
    disp('Sparse Inverse Covariance (QUIC optimizer) used');
    
else
    
    error('No valid method provided')
    
end

% Check transformation type
if strcmp(opts.transType,'clr')
    disp('CLR transformation applied')
    Xin = clr(X);
elseif strcmp(opts.transType,'none')
    disp('No data transformation applied')
    Xin = X;
else
    disp('No valid data transformation given; using identity')
    Xin = X;
    
end

% Check whether (inverse) correlation or covariance should be estimated
if strcmp(opts.matType,'corrMat')
    
    disp('Correlation matrix used')
    
    % Sample correlation matrix
    S = corr(Xin');
    
elseif strcmp(opts.matType,'covMat')
    
    disp('Covariance matrix used')
    
    % Sample covariance matrix
    S = cov(Xin')
    
else
    error('Unknown matrix type')
    
end

if ~isempty(opts.lamMax)
    
    lamMax = opts.lamMax
    
else
    
    
    
    % Maximum value of the covariance matrix
    lamMax = 1/2*max(abs(S(:)));
    
    % Maximum lambda value for neighborhood selection
    if strcmp(opts.matType,'covMat') && ~(strcmp(opts.method,'quic'))
        lamMax = 0;
        for i=1:p
            ytemp = Xin(i,:);
            Xtemp = Xin([1:i-1,i+1:p],:);
            temp = max(abs(ytemp*Xtemp'));
            if temp>lamMax
                lamMax=temp;
            end
        end
        % Fix me
        lamMax = 1/2*lamMax;
        
    end
    
    % Maximum lambda value for TREX in neighborhood setting
    if strcmp(opts.method,'trex_or') || strcmp(opts.method,'trex_and')
        % Fix me
        lamMax = 1.5;
    end
    
end

% Minimum Lambda depends on maximum lambda
lamMin = opts.lamMinFac.*lamMax;

% Create actual lambda path
numLams = opts.nLams;
lamVec = linspace(lamMin,lamMax,numLams)
opts.lambdaList = lamVec;

disp('StARS options')

opts

% StARS model selection
starsOut = stars(Xin, opts);

covMat = starsOut.optSigma;
iCovMat = starsOut.optTheta;
temp = (iCovMat~=0).*1;

sizVec = size(temp);
if length(sizVec)==2
    % Remove diagonal elements
    temp(1:p+1:p^2)=0;
    G = (temp~=0).*1;
else
    for i=1:sizVec(3)
        temp2 = squeeze(temp(:,:,i));
        temp2(1:p+1:p^2)=0;
        temp(:,:,i) = temp2;
    end
    G = temp;
end


% ---------------------------------------------------------------
% FUNCTIONS BELOW ARE TAKEN FROM Niko Hansen's CMA-ES code
% ---------------------------------------------------------------
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%


if nargin < 2 || isempty(defopts) % no default options available
    opts=inopts;
    return;
elseif isempty(inopts) % empty inopts invoke default options
    opts = defopts;
    return;
elseif ~isstruct(defopts) % handle a single option value
    if isempty(inopts)
        opts = defopts;
    elseif ~isstruct(inopts)
        opts = inopts;
    else
        error('Input options are a struct, while default options are not');
    end
    return;
elseif ~isstruct(inopts) % no valid input options
    error('The options need to be a struct or empty');
end

opts = defopts; % start from defopts
% if necessary overwrite opts fields by inopts values
defnames = fieldnames(defopts);
idxmatched = []; % indices of defopts that already matched
for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
    if isoctave
        for i = 1:size(defnames, 1)
            idx(i) = strncmp(lower(defnames(i)), lower(name), length(name));
        end
    else
        idx = strncmp(lower(defnames), lower(name), length(name));
    end
    if sum(idx) > 1
        error(['option "' name '" is not an unambigous abbreviation. ' ...
            'Use opts=RMFIELD(opts, ''' name, ...
            ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
        defname  = defnames{find(idx)};
        if ismember(find(idx), idxmatched)
            error(['input options match more than ones with "' ...
                defname '". ' ...
                'Use opts=RMFIELD(opts, ''' name, ...
                ''') to remove the field from the struct.']);
        end
        idxmatched = [idxmatched find(idx)];
        val = getfield(inopts, name);
        % next line can replace previous line from MATLAB version 6.5.0 on and in octave
        % val = inopts.(name);
        if isstruct(val) % valid syntax only from version 6.5.0
            opts = setfield(opts, defname, ...
                getoptions(val, getfield(defopts, defname)));
        elseif isstruct(getfield(defopts, defname))
            % next three lines can replace previous three lines from MATLAB
            % version 6.5.0 on
            %   opts.(defname) = ...
            %      getoptions(val, defopts.(defname));
            % elseif isstruct(defopts.(defname))
            warning(['option "' name '" disregarded (must be struct)']);
        elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
            opts = setfield(opts, defnames{find(idx)}, val);
            % next line can replace previous line from MATLAB version 6.5.0 on
            % opts.(defname) = inopts.(name);
        end
    else
        warning(['option "' name '" disregarded (unknown field name)']);
    end
end

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function res = isoctave
% any hack to find out whether we are running octave
s = version;
res = 0;
if exist('fflush', 'builtin') && eval(s(1)) < 7
    res = 1;
end

