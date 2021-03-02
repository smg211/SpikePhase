function [pval, phz_shuff, coh_shuff] = coh_bootstrap(phz, nshuffs)
% Get the bootstrapped p-value for the coherence of a dataset given phase

coherence_real = get_coherence(phz);

if size(phz, 2) == 1
  phz = phz';
end

coh_shuff = [];
phz_shuff = {};
for shuffle = 1:nshuffs
  
  % Jitter each spk by a random amomunt +/- pi
  jitter = rand(1, length(phz)); % uniform distribution b/w [0, 1]
  jitter_centered = (jitter - 0.5); % jitter from [-.5, .5]
  jitter_amt = jitter_centered*2*pi; % jitter from [-.5*cycle, .5*cycle]
  
  phz_jitter = phz + jitter_amt;
  
  % Make sure phase is between 0 and 2*pi
  phz_jitter(sign(phz_jitter) == -1) = phz_jitter(sign(phz_jitter) == -1) + 2*pi;
  phz_jitter(phz_jitter > 2*pi) = phz_jitter(phz_jitter > 2*pi) - 2*pi;
  
  % Get the phase of this shuffle
  phz_shuff{shuffle} = phz_jitter;
  
  % Get the shuffled coherence:
  coh_shuff = [coh_shuff get_coherence(phz_jitter);];
end

pval = 1 - (sum(coherence_real > coh_shuff)/nshuffs);