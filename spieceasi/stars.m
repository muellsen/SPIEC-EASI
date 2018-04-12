function out = stars(X,inopts)
% StARS model selection procedure (adapted from Liu, H., Roeder, K., & Wasserman, L. (2010).
% Stability approach to regularization selection (StARS) for high
% dimensional graphical models. NIPS)
%
% Input:    X: pxn matrix of count Xdata where p is number of OTUs, n is number of samples
%           inopts: option structure from spieceasi.m
% Output:   out:  Output structure of the StARS selection procedure (further processed in spieceasi)
%
% Author: CL Mueller, Simons Foundation
%         ZD Kurtz, NYU
% When using this code please also cite:
%
% Kurtz, Z. D., Müller, C. L., Miraldi, E. R., Littman, D. R., Blaser, M. J., & Bonneau, R. A. (2015).
% Sparse and Compositionally Robust Inference of Microbial Ecological Networks. PLoS Computational Biology, 11(5), e1004226.
% http://doi.org/10.1371/journal.pcbi.1004226
%

% Dimension of the input data
[p,n] = size(X);

% Necessary options for stars code
if ~isfield(inopts,'lambdaList')
    error('Regularization path needs to be provided')
else
    defopts.lambdaList = inopts.lambdaList;
end

% Options for StARS model selection
defopts.verboseMode = 0;                % Verbose during computation
defopts.method = 'quic';                % Standard graph inference method
defopts.repNum = 10;                    % Number of repetitions for StARS
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
    opts = getoptions(inopts,defopts);
end

% Augment stars options from SPIEC-EASI inopts
fNames = fieldnames(inopts);
for s = 1:length(fNames)
    if ~isfield(opts,fNames{s})
        opts = setfield(opts,fNames{s},getfield(inopts,fNames{s}));
    end
end

% Transpose of the input data X
Xdata = X';

% Dimension of Xdata
[n,p] = size(Xdata)

starsSubsampleRatio = opts.starsSubsampleRatio;

if floor(n * starsSubsampleRatio)>=n
    error('Sample size not sufficient for defined subsample ratio')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reversoptsr = '';

% Set opts fields to variables for potential multi-thread StARS
lambdaList = opts.lambdaList;
repNum = opts.repNum;
starsThresh = opts.starsThresh;

% Matrix that counts occurrence of edges across subsamples and the entire
% lambda path
stabMat = zeros(p, p, length(lambdaList));

gopts = [];
gopts.lambdaList = opts.lambdaList;
gopts.verboseMode = opts.verboseMode;
gopts.matType = opts.matType;

verbose = gopts.verboseMode;

for i = 1:repNum
    if verbose
        tic
    end
    indSample = randsample(1:n, floor(n * starsSubsampleRatio));
    %%% Which method to use
    if strcmp(opts.method,'mb_or')
        
        [~,~,G_union,~,~] = mblasso(Xdata(indSample,:), gopts);
        
        tempG = G_union;
        
    elseif strcmp(opts.method,'mb_and')
        [~,~,~,G_inter,~] = mblasso(Xdata(indSample,:), gopts);
        
        tempG = G_inter;
        
    elseif strcmp(opts.method,'quic')
        
        gopts.tolThreshold = opts.tolThreshold;
        gopts.maxIter = opts.maxIter;
        
        [~, ~, tempG,~, ~] = QUICwrap2(Xdata(indSample,:), gopts);
        
    elseif strcmp(opts.method,'trex_or')
        
        [~,~,G_union,~,~] = mbtrex(Xdata(indSample,:), gopts);
        
        tempG = G_union;
        
    elseif strcmp(opts.method,'trex_and')
        [~,~,~,G_inter,~] = mbtrex(Xdata(indSample,:), gopts);
        
        tempG = G_inter;
    else
        error('Method not supported')
    end
    
    if verbose
        msg = sprintf('Processed %s', int2str(floor(100*(i/repNum))));
        fprintf('\n');
        fprintf([reversoptsr, msg, '%%\n']);
        reversoptsr = repmat(sprintf('\b'), 1, length(msg)+1);
    end
    
    if verbose
        timePerSubsammple = toc
    end
    
    % Add edge counts to stability matrix
    stabMat = stabMat + tempG;
    
    % Compute relative edge frequency
    normStabMat       = stabMat / i;
    
    % Compute variabilty measure
    variability = 4 * sum(sum(normStabMat .* (1 - normStabMat))) / (p * (p - 1));
    variability = squeeze(variability);
    % Monotonize variability
    for j=1:length(variability)-1
        if variability(j)<variability(j+1)
            variability(j)=variability(j+1);
        end
    end
    
    % Identify optimal lambda using starsThresh value
    ind           = find(variability < starsThresh);
    [Vsel, VmaxI] = max(variability(ind)); % closest value to starsThresh from below
    optInd        = ind(VmaxI);
    optLam        = lambdaList(optInd);
    
    if verbose
        figure(42)
        plot(lambdaList,variability)
        grid on
        hold on
        xlabel('lambda path')
        ylabel('Variability')
        drawnow
    end
    
    % Use heuristic lower bound after two subsamples
    if i==2 && opts.uselamLB
        % Set small lambdas to optimal lambda in the two-sample case
        gopts.lambdaList(1:optInd) = optLam;
    end
end

if verbose
    fprintf('\n');
end

if repNum > 0
    % Monotonize variability
    for i=1:length(variability)-1
        if variability(i)<variability(i+1)
            variability(i)=variability(i+1);
        end
    end
    
    % Identify optimal lambda using starsThresh value
    ind           = find(variability < starsThresh);
    [Vsel, VmaxI] = max(variability(ind)); % closest value to starsThresh from below
    optInd        = ind(VmaxI);
    optLam        = lambdaList(optInd);
    
    if isempty(optLam)
        variability
        warning('No regularization parameter found that meets StARS selection criterion from below')
        ind           = find(variability > starsThresh);
        [Vsel, VminI] = min(variability(ind)); % closest value to starsThresh from below
        optInd        = ind(VminI);
        optLam        = lambdaList(optInd);
        
    end
end
% Final optimization with selected StARS lambda
fopts = gopts;

% When repNum==0, StARS returns the standard trace without subsampling
if repNum > 0
    fopts.lambdaList = optLam; % Run it with optimal lambda
end

%%% Final network estimate using the estimated lambda
if strcmp(opts.method,'mb_or')
    
    [w_union,~] = mblasso(Xdata, fopts);
    Theta = w_union;
    Sigma = [];
    
elseif strcmp(opts.method,'mb_and')
    
    [~,w_inter] = mblasso(Xdata, fopts);
    Theta = w_inter;
    Sigma = [];
    
elseif strcmp(opts.method,'trex_or')
    
    [w_union,~] =  mbtrex(Xdata, fopts);
    Theta = w_union;
    Sigma = [];
    
elseif strcmp(opts.method,'trex_and')
    
    [~,w_inter] = mbtrex(Xdata, fopts);
    Theta = w_inter;
    Sigma = [];
    
elseif strcmp(opts.method,'quic')
    
    fopts.tolThreshold = opts.tolThreshold;
    fopts.maxIter = opts.maxIter;
    [Sigma,Theta] = QUICwrap2(Xdata, fopts);
    
else
    error('Method not supported')
end

% Output all relevant computations and parameters
out.lambdaList = lambdaList;
out.opts = opts;
out.gopts = gopts;
out.fopts = fopts;
out.optSigma = Sigma;
out.optTheta = Theta;

% StARS output for subsampling
if repNum > 0
    out.stabMat    = stabMat;
    out.varSel   = Vsel;
    out.optInd   = optInd;
    out.optLam   = lambdaList(optInd);
    out.variability = variability;
    out.stabMat = stabMat;
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


