function coh_norm = get_coherence_norm(phz, n_subsmp, size_subsmp)
% Get the average coherence of repeated subsamples of the input phase

if any(phz < 0 | phz > 2*pi)
  warning('get_coherence_norm is supposed to receive phase in radians as input, but this input has values outside of [0, 2*pi] \n');
end

n_phz_in = length(phz);
if n_phz_in >= size_subsmp
  
  % initialize the vector for saving the coherence of each repetition
  coh_rep = [];
  for rep = 1:n_subsmp
    % get the sub sample of inputs to use for calculating coherence
    % for this repetition
    i_phz_rep = randperm(n_phz_in, size_subsmp);
    phz_rep = phz(i_phz_rep);
    coh_rep(rep) = get_coherence(phz_rep);
  end % end coherence calculation repetition loop
  
  % take the median of the distribution of subsampled coherence
  coh_norm = median(coh_rep, 'omitnan');
else
  warning('The size of the desired sub sample is larger than the number of inputs. Cannot get count-normalized coherence.');
  coh_norm = nan;
end % end if n_spikes2ana > size_subsmp