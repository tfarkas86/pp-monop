%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metacommunity predator-prey model with evolving partners
% model aims to look at effects of predator evolution on coexistence of
% prey, and or predators
%
% writen by Tim Farkas, 2015/16 after Mark Urban model 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize script

clear all

%% initialize habitat conditions

nhabs = 1; % number of habitats
npredsites = 100; % maximum number of predators per habitat
maxprey = 9; % size of predator territory (max number of prey attacked)
npreysites = npredsites * maxprey; % number of microsites per habitat

%% initialize population parameters

npredsp = 1; % number of predator species
npreysp = 1; % number of prey species
maxint = npredsp * npreysp * npreysites; % max p-p interactions

mort = [0.2] .* ones(1, npredsp); % baseline predator mortality rates

% initialize predator-prey conversion efficiencies
ceff = [0.1] .* ones(npredsp, npreysp);

% attack rates, first as 2D matrix with prey cols and pred rows
% then as 3D array with sites as rows, prey cols, and pred z-dim
attkm = [0.9] .* ones(npredsp, npreysp); % m for matrix
attki = permute(attkm, [3 2 1]); % make 3D (i for intermediate)
attka = repelem(attki, npreysites, 1, 1); % individual attack rates (a for array)

birt = [.1] .* ones(1, npreysp); % prey birth rates

%% initialize matrices to follow occupancy in sites

preyocc = zeros(npreysites, npreysp); % occupancy of prey at their sites
predoccp = zeros(npredsites, npredsp); % occupancy of predators

nprey = ceil(npreysites / 2 / npreysp * ones(1, npreysp)); % starting number of prey
npred = ceil(npredsites / 10 / npredsp * ones(1, npredsp)); % starting number of predators

% set initial locations of prey and predators
preyloc = datasample(1:npreysites, sum(nprey), 'Replace', false); % prey
sampvec = repelem(1:npreysp, nprey);
lind = preyloc + (sampvec - 1) * size(preyocc, 1); % get linear indices
preyocc(lind) = 1; % fill in prey

predloc = datasample(1:npredsites, sum(npred), 'Replace', false); % predators
sampvec = repelem(1:npredsp, npred);
lind = predloc + (sampvec - 1) * size(predoccp, 1); % get linear indices
predoccp(lind) = 1; % fill in predators

% make prey-site sized matrix showing which sites are threatened by the
% presence of a predator

predocc = permute(predoccp, [1 3 2]); % z-dimension predator species
predocc3D = repelem(predocc, maxprey, 1, 1);
predocc = permute(predocc3D, [1, 3, 2]);
predoccrep = repmat(predocc3D, 1, npreysp, 1); % expand for all prey sp.
preyoccrep = repmat(preyocc, 1, 1, npredsp);
%% start time loop
time = 100;

% output matrices
outpred = zeros(time, npredsp + 1); % open predator matrix
outprey = zeros(time, npreysp + 1); % open prey matrix
outpred(:, 1) = time; % input time
outprey(:, 1) = time; % input time
outpred(1, 2:npredsp+1) = npred; % input starting predator abundances
outprey(1, 2:npreysp+1) = nprey; % input starting prey abundances

for i = 2:time;
    
    %% predators eat some prey
    
    % expand for pred sp.
    meet = preyoccrep .* predoccrep; % check for predator-prey encounters
    
    % random numbers to check against attack rate
    thresh = reshape(rand(maxint, 1), [npreysites, npreysp, npredsp]);
    predate = thresh .* meet > (1 - attka); % find successful predation
    preyocc(find(sum(predate, 3) == 1)) = 0; % kill eaten prey
    
    % number each prey sp. (cols) eat by each pred sp. (rows)
    neat = permute(sum(predate, 1), [3, 2, 1]);
    
    % background predator mortality
    % need to make this a fitness based mortality
    
    predoccp(find(rand(npredsites, npredsp) .* predoccp > ...
        (1 - mort(:, 1))),1) = 0;
    predocc = permute(predoccp, [1 3 2]); % z-dimension predator species
    predocc3D = repelem(predocc, 2, 1, 1);
    predocc = permute(predocc3D, [1, 3, 2]);
    
    %% prey and predators have offspring and randomly fill sites
    
    % prey
    openprey = find(sum(preyocc, 2) == 0); % indices for open sites
    preyoff = ceil(birt .* sum(preyocc, 1)); % number of prey offspring
    maxoff = min([length(openprey), sum(preyoff)]); % max no. offspring
    sampvec = repelem(1:npreysp, preyoff);
    
    if length(openprey) > sum(preyoff) % more sites than offspring
        
        % randomly select sites to fill
        openprey = datasample(openprey, maxoff, 'Replace', false);
        
    elseif length(openprey) < sum(preyoff) % more offspring than sites
        
        % randomly select offspring to populate site
        sampvec = datasample(sampvec, maxoff, 'Replace', false);
        
    end
    
    lind = openprey' + (sampvec - 1) * size(preyocc, 1); % get linear indices
    preyocc(lind) = 1; % fill in prey
    
    % predators
    openpred = find(sum(predoccp, 2) == 0); % indices for open sites
    predoff = sum(ceil(ceff .* neat), 2); % number of predator offspring
    maxoff = min([length(openpred) sum(sum(predoff))]);
    sampvec = repelem(1:npredsp, predoff);
    
    if length(openpred) > sum(predoff) % more sites than offspring
        
        % randomly select sites to fill
        openpred = datasample(openpred, maxoff, 'Replace', false);
        
    elseif length(openpred) < sum(predoff) % more offspring than sites
        
        % randomly select offspring to populate site
        sampvec = datasample(sampvec, maxoff, 'Replace', false);
        
    end
    
    lind = openpred' + (sampvec - 1) * size(predoccp, 1); % get linear indices
    predoccp(lind) = 1; % fill in prey
    
    outpred(i, 2:npredsp + 1) = sum(predoccp);
    outprey(i, 2:npredsp + 1) = sum(preyocc);
end

plot(1:time, outpred(:,2), 1:time, outprey(:,2))