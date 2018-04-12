function [w_union,w_inter,G_union,G_inter,w_raw] = mbtrex(X,inopts)
% Neighborhood selection (Meinshausen-Buehlmann) approach with the TREX
% 
% Input:  Data X in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters (as set in e.g.
%         spieceasi)
% Output: all outputs are p x p x length of lasso path (numLams)
%         w_union: Sym. weights in the graph using the union 
%         w_inter: sym. weights in the graph using the intersection 
%         G_union: Graph using the union of p lasso problems 
%         G_inter: Graph using the intersection of p lasso problems 
%         w_raw:   Raw weights in the graph from p lasso problems 
%
% When using this code please cite:
%
% Lederer, J. and Müller, C. L. (2014). Topology Adaptive Graph Estimation in High Dimensions.
% arXiv preprint arXiv:1410.7279.


% Options for MBtrex
defopts.verboseMode = 0;
defopts.lambdaList = linspace(0.01,2,100);
defopts.matType = 'corrMat';

% Merge options inopts and defopts
if nargin < 2 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

% Dimension of the input data
[n,p] = size(X);

% Normalize data if necessary
if strcmp(opts.matType,'corrMat')
    
    X = (X-repmat(mean(X),n,1));
    
    normX = repmat(sqrt(sum(X.^2)),n,1);
    X = sqrt(n)*X./normX;
    
elseif strcmp(opts.matType,'covMat')
    
    X = (X-repmat(mean(X),n,1));
    
else
    error('Matrix type not specified')
end

% Set verbose flag
verbose = opts.verboseMode;

% lambda vector from max lambda to min lambda for warmstart
[lamVec,sortedLamInds] = sort(opts.lambdaList,'descend');

numLams = length(lamVec);

w_raw = zeros(p,p,numLams);
temp = zeros(p-1,p,numLams);
B_p = zeros(p,p,numLams);

for i=1:p
    
    if verbose
        disp('######################################')
        disp(['Regression on the ',num2str(i),'th variable'])
        disp('######################################')
    end
    
    currInds = setdiff(1:p,i);
    
    currX = X(:,currInds);
    Yvec = X(:,i);
        
    trexopts.beta0 = zeros(p-1,1);
    trexopts.boots= 0;
    trexopts.rep = 2;
    trexopts.cpath = lamVec;
    tempBeta = trexp(currX,Yvec,trexopts);
    temp(:,i,:) = tempBeta(:,sortedLamInds);

end

for i=1:p
    currInds = setdiff(1:p,i);
    B_p(currInds,i,:)=temp(:,i,:);
end

if verbose
    fprintf('\n')
    
    disp('######################################')
    disp(['Regression completed'])
    disp('######################################')
    fprintf('\n')
    
end


G_union = zeros(p,p,numLams);
G_inter = zeros(p,p,numLams);

w_union = zeros(p,p,numLams);
w_inter = zeros(p,p,numLams);

for j=1:numLams
    
    i=sortedLamInds(j);
    
    currB = squeeze(B_p(:,:,i));
    temp = abs(currB);
    signMat = sign(currB+currB');
    
    w_union(:,:,i) = signMat.*max(temp',temp);
    w_inter(:,:,i) = signMat.*min(temp',temp);
    G_union(:,:,i) = (w_union(:,:,i)~=0);
    G_inter(:,:,i) = (w_inter(:,:,i)~=0);
    
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


