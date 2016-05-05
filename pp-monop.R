# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# metacommunity predator-prey model with evolving partners
# model aims to look at effects of predator evolution on coexistence of
# prey, and or predators
#
# writen by Tim Farkas, 2015/16 after Mark Urban model 2009
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Model Parameters ----

npreysp <- 2 # number of prey species
npredsp <- 2 # number of predator species

birt <- rep(0.1, npreysp) # prey birth rate

mort <- rep(0.1, npredsp) # predator mortality rate

# ceff <- c(0.1, 0.1) # conversion efficiencies, cols prey, row preds
ceff <- matrix(0.1, nrow=npredsp, ncol=npreysp)

# attkm <-  c(.9; .1) # attack rate matrix, cols prey, rows predators
attk <- matrix(0.5, nrow=npredsp, ncol=npreysp)
attk <- .5

nalls <- 20 # number of alleles

nsites <- 50 # number of prey sites
maxpred <- 50 # max number of predators per species

afqprey <- rep(.5, npreysp) # allele frequencies for prey
afqpred <- rep(.5, npredsp) # allele frequencies for predators

# make predator and prey matrices ----

# prey matrix has rows a sites to be filled with prey
# column 1 specific species, other 20 specify alleles at 20 loci
preymat <- matrix(0, nrow=nsites, ncol=nalls + 1, 
                  dimnames=list(sites=NULL, attributes=c("spp", 1:20)))

# predator array has unlimiting spaces for predators (rows), alleles (cols)
# and 3rd dimension is for multiple predator species
predmat <- array(0, dim=c(maxpred, nalls + 1, npredsp), 
                 dimnames=list(predator=NULL, allele=NULL, species=NULL))

# specify number of prey and predators
# number of prey is 1/10 number of sites
nprey <- rep(ceiling(nsites/npreysp/10), npreysp)
npred <- rep(ceiling(sum(nprey)/10), npredsp)

# fill matrices ----

ipreysp <- rep(1:npreysp, nprey) # initial prey species identities
preymat[1:length(ipreysp), 1] <- ipreysp # add species to prey matrix

ipredsp <- matrix(1, nrow=npred, npredsp) # initial predator species 
predmat[1:nrow(ipredsp), 1, 1:ncol(ipredsp)] <- ipredsp # add predators


# makeMats function creates initial matrices for predators and prey
# and take as arguments npredsp and npreysp (number of predator and prey 
# species), nsites (number of microsites occupiable by prey), imaxpred 
# (arbitrarily large matrix length to hold predator individuals), nprey 
# (number of initial prey individuals), npred (number of initial predators), 
# and nalls (number of alleles)

makeMats <- function(npredsp=1, npreysp=1, nsites=50, imaxpred=50, 
                     nprey=rep(ceiling(nsites/npreysp/10), npreysp),
                     npred=rep(ceiling(sum(nprey)/10), npredsp), 
                     nalls=20) {
  
  # preymat is 2D matrix with locations either empty or filled by species
  # each individual with 20 loci
  preymat <- matrix(0, nrow=nsites, ncol=nalls + 1, 
                    dimnames=list(sites=NULL, attributes=c("spp", 1:nalls)))
  # predator matrix is 3D array with rows as locations, columns as occupancy 
  # (0 or 1) and 20 loci, and depth as predator species
  predmat <- array(0, dim=c(imaxpred, nalls + 1, npredsp), 
                   dimnames=list(predator=NULL, attributes=c("occ", 1:nalls), 
                                 species=NULL))
  
  ipreysp <- rep(1:npreysp, nprey) # initial prey species identities
  preymat[1:length(ipreysp), 1] <- ipreysp # add species to prey matrix
  
  ipredsp <- matrix(1, nrow=npred, npredsp) # initial predator species 
  predmat[1:nrow(ipredsp), 1, 1:ncol(ipredsp)] <- ipredsp # add predators
  
  mats <- list(preymat=preymat, predmat=predmat)
  return(mats)
  
}

npreds <- sum(mats$predmat[,1,]) # how many predators of each species
npreys <- sum(mats$preymat[mats$preymat[,1] > 0,1]) # how many total prey?

# encrt is function to calculate proportion of prey that will be encountered 
# by at least one predator. uses recursion to subtract out probability of 
# predators overlapping when penal = TRUE

encrt <- function(mat, arr=1, penal=FALSE) {
  
  pred <- sum(mat$predmat[,1,])
  prey <- sum(mat$preymat[mat$preymat[,1] > 0, 1])
  locs <- nrow(mat$preymat)
  
  if (penal==TRUE) {
    
    if (arr == pred) {
      return(choose(pred, arr) * (prey / locs) ^ arr) }
    
    else {
      return(choose(pred, arr) * (prey / locs) ^ arr - 
               encrt(mat=mat, arr=arr + 1))} }
  
  else { 
    return(pred * prey / locs)
  }
  
}

mat <- makeMats(nprey=20, npred=5)

# killprey function takes a matrix of predators and prey, compares phenotypes 
# using allele information, integrates a baseline attack rate (attk) and kills
# prey in the preymat by setting all values to 0, returning a full object
# with both predator and prey matrices

killprey <- function(mat, attk=.5) {
  
  nenc <- ifelse(encrt(mat=mat) < 1, rbinom(1, 1, prob=encrt(mat=mat)), 
                 round(encrt(mat=mat))) # number of predator-prey encounters
  
  preyind <- which(mat$preymat[,1] > 0, arr.ind=TRUE) # which locations are inhabited by prey
  predind <- which(mat$predmat[,1,] > 0, arr.ind=TRUE) # by predators
  
  emat <- cbind(prey=sample(preyind, nenc), pred=sample(predind, nenc)) # encounter matrix
  amat <- attk +  rowMeans((mat$preymat[emat[,1], 2:21])) -
    rowMeans((mat$predmat[emat[,2], 2:21,])) # attack rate matrix
  emat <- cbind(emat, kill=amat > runif(nenc, 0, 1)) # find successful attacks
  
  mat$preymat[emat[,1],] <- mat$preymat[emat[,1],] * (1 - emat[,3]) 
  
  return(mat)
  
}












# attkmin = 0.1; % minimum attack rate
# attkmax = 0.9; % maximum attack rate
# 
# birt = [.1, .3]; % prey birth rates
# %birt = [0.1] .* ones(1, npreysp);
# 
# % starting number of prey and predators
# % nprey = [200]
# npred = [3, 3];
# 
# nalls = 10; % number of alleles
# 
# mutr = 0.01; % mutation rate
# 
# %% derived quantities
# 
# npreysites = npredsites * maxprey; % number of microsites per habitat
# % nprey = ceil(npreysites / 2 / npreysp * ones(1, npreysp)); 
# % npred = ceil(npredsites / 10 / npredsp * ones(1, npredsp)); 
# maxint = npredsp * npreysp * npreysites; % max p-p interactions
# 
# % make 3D attack rate array, where rows now are sites, cols are prey, 
# % and z-dimension is predator species (a for array)
# attka = repelem(permute(attkm, [3 2 1]), npreysites, 1, 1);
# 
# attkran = attkmax - attkmin; % range of attack rates
# 
# %% initialize matrices to follow occupancy/abundance in sites
# 
# preyocc = zeros(npreysites, npreysp); % occupancy of prey at their sites
# predoccp = zeros(npredsites, npredsp); % occupancy of predators
# 
# % set initial locations of prey and predators
# preyloc = datasample(1:npreysites, sum(nprey), 'Replace', false); % prey
# sampvec = repelem(1:npreysp, nprey);
# lind = preyloc + (sampvec - 1) * size(preyocc, 1); % get linear indices
# preyocc(lind) = 1; % fill in prey
# 
# predloc = datasample(1:npredsites, sum(npred), 'Replace', false); % predators
# sampvec = repelem(1:npredsp, npred);
# lind = predloc + (sampvec - 1) * size(predoccp, 1); % get linear indices
# predoccp(lind) = 1; % fill in predators
# 
# % make prey-site sized matrix showing which sites are threatened by the
# % presence of a predator
# 
# %% make matrices to follow alleles and phenotypes
# 
# % this has to be the most inefficient way of doing this ...
# zone = [iphen * nalls, (1 - iphen) * nalls]; % get ones and zeros vector
# ind = repmat([1 3], 1, npredsp) + repelem(0:(npredsp - 1), 2); % change index
# zone = repelem(repmat([1 0], 1, npredsp), zone(ind)); % fill ones and zeros
# allpred = vec2mat(zone, nalls); % alleles cols, predators z-dim
# allpred = permute(allpred, [3, 2, 1]);
# allpred = repelem(allpred, npredsites, 1, 1)
# 
# 
# %% start time loop
# time = 2000;
# 
# % output matrices
# outpred = zeros(time, npredsp + 1); % open predator matrix
# outprey = zeros(time, npreysp + 1); % open prey matrix
# outpred(:, 1) = time; % input time
# outprey(:, 1) = time; % input time
# outpred(1, 2:npredsp+1) = npred; % input starting predator abundances
# outprey(1, 2:npreysp+1) = nprey; % input starting prey abundances
# 
# for i = 2:time;
# 
# %% reformat predator and prey matrices
# 
# predocc = permute(predoccp, [1 3 2]); % z-dimension predator species
# predocc3D = repelem(predocc, maxprey, 1, 1);
# predocc = permute(predocc3D, [1, 3, 2]);
# 
# predoccrep = repmat(predocc3D, 1, npreysp, 1); % expand for all prey sp.
# 
# % expand for pred sp.
# preyoccrep = repmat(preyocc, 1, 1, npredsp);
# 
# %% predators eat some prey
# 
# meet = preyoccrep .* predoccrep; % check for predator-prey encounters
# 
# % random numbers to check against attack rate
# thresh = reshape(rand(maxint, 1), [npreysites, npreysp, npredsp]);
# predate = thresh .* meet > (1 - attka); % find successful predation
# preyocc(find(sum(predate, 3) == 1)) = 0; % kill eaten prey
# 
# % number each prey sp. (cols) eat by each pred sp. (rows)
# neat = permute(sum(predate, 1), [3, 2, 1]);
# 
# % background predator mortality
# % need to make this a fitness based mortality
# 
# predoccp(find(rand(npredsites, npredsp) .* predoccp > ...
#               (1 - mort(:, 1))),1) = 0;
# predocc = permute(predoccp, [1 3 2]); % z-dimension predator species
# predocc3D = repelem(predocc, 2, 1, 1);
# predocc = permute(predocc3D, [1, 3, 2]);
# 
# %% prey and predators have offspring and randomly fill sites
# 
# % prey
# openprey = find(sum(preyocc, 2) == 0); % indices for open sites
# preyoff = ceil(birt .* sum(preyocc, 1)); % number of prey offspring
# maxoff = min([length(openprey), sum(preyoff)]); % max no. offspring
# sampvec = repelem(1:npreysp, preyoff);
# 
# if length(openprey) > sum(preyoff) % more sites than offspring
# 
# % randomly select sites to fill
# openprey = datasample(openprey, maxoff, 'Replace', false);
# 
# elseif length(openprey) < sum(preyoff) % more offspring than sites
# 
# % randomly select offspring to populate site
# sampvec = datasample(sampvec, maxoff, 'Replace', false);
# 
# end
# 
# lind = openprey' + (sampvec - 1) * size(preyocc, 1); % get linear indices
# preyocc(lind) = 1; % fill in prey
# 
# % predators
# openpred = find(sum(predoccp, 2) == 0); % indices for open sites
# predoff = sum(ceil(ceff .* neat), 2); % number of predator offspring
# maxoff = min([length(openpred) sum(sum(predoff))]);
# sampvec = repelem(1:npredsp, predoff);
# 
# if length(openpred) > sum(predoff) % more sites than offspring
# 
# % randomly select sites to fill
# openpred = datasample(openpred, maxoff, 'Replace', false);
# 
# elseif length(openpred) < sum(predoff) % more offspring than sites
# 
# % randomly select offspring to populate site
# sampvec = datasample(sampvec, maxoff, 'Replace', false);
# 
# end
# 
# lind = openpred' + (sampvec - 1) * size(predoccp, 1); % get linear indices
# predoccp(lind) = 1; % fill in prey
# 
# outpred(i, 2:npredsp + 1) = sum(predoccp);
# outprey(i, 2:npreysp + 1) = sum(preyocc);
# end
# 
# % plot abundance over time
# close
# hold on
# plot(1:time, outprey(:, 2), '-r')
# plot(1:time, outprey(:, 3), '-b')
# plot(1:time, outpred(:, 2), '--k')
