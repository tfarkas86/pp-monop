%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metacommunity predator-prey model with evolving partners
% model aims to look at effects of predator evolution on coexistence of
% prey, and or predators
%
% writen by Tim Farkas, 2015/16 after Mark Urban model 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize script

clear all

%% input model parameters

npreysp = 1; % number of prey species
npredsp = 1; % number of predator species

npredsites = 100; % maximum number of predators 
maxprey = 9; % size of predator territory (max number of prey to attack)

mort = [0.1]; % baseline predator mortality rates
%mort = [0.1] .* ones(1, npredsp);

ceff = [0.1]; % predator-prey conversion efficiencies, cols prey, row pred
%ceff = [0.1] .* ones(npredsp, npreysp);

attkm = [0.5]; % attack rate matrix, cols prey, rows predators
%attkm = [0.5] .* ones(npredsp, npreysp); 

birt = [0.1]; % prey birth rates
%birt = [0.1] .* ones(1, npreysp);

% starting number of prey and predators
% nprey = [200]
% npred = [20]

%% derived quantities

npreysites = npredsites * maxprey; % number of microsites per habitat
nprey = ceil(npreysites / 2 / npreysp * ones(1, npreysp)); 
npred = ceil(npredsites / 10 / npredsp * ones(1, npredsp)); 
maxint = npredsp * npreysp * npreysites; % max p-p interactions

% make 3D attack rate array, where rows now are sites, cols still prey, 
% and z-dimension is predator species (a for array)
attka = repelem(permute(attkm, [3 2 1]), npreysites, 1, 1); 

%% initialize matrices to follow occupancy/abundance in sites

preyocc = zeros(npreysites, npreysp); % occupancy of prey at their sites
predoccp = zeros(npredsites, npredsp); % occupancy of predators

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

%% start time loop
time = 1000;

% output matrices
outpred = zeros(time, npredsp + 1); % open predator matrix
outprey = zeros(time, npreysp + 1); % open prey matrix
outpred(:, 1) = time; % input time
outprey(:, 1) = time; % input time
outpred(1, 2:npredsp+1) = npred; % input starting predator abundances
outprey(1, 2:npreysp+1) = nprey; % input starting prey abundances

for i = 2:time;
    
    %% reformat predator and prey matrices
    
    predocc = permute(predoccp, [1 3 2]); % z-dimension predator species
    predocc3D = repelem(predocc, maxprey, 1, 1);
    predocc = permute(predocc3D, [1, 3, 2]);
   
    predoccrep = repmat(predocc3D, 1, npreysp, 1); % expand for all prey sp.
   
    % expand for pred sp.
    preyoccrep = repmat(preyocc, 1, 1, npredsp);
    
    %% predators eat some prey
      
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

hold on
plot(1:time, outprey(:,2))
plot(1:time, outpred(:,2))
hold off