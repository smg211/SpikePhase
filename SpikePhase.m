function spike_phase = SpikePhase(cfg, dat, spike_metrics)

% SPIKEPHASE takes in an ACS dataset and returns the following metrics for
% each unit, u, and each condition, cond (fieldnames corresponding to each
% metric in the output in parentheses):
%   1. Raw coherence of all spikes (coh)
%   2. Coherence normalized by a user defined number of spikes (coh_Nnormspk)
%   3. Coherence normalized by the minimum number of spikes across
%      all conditions (coh_condnormspk)
%   4. P-value from a bootstrapped permutation test of all spikes (p_boot)
%   5. P-value from a rayleigh test for unimodality of all spikes (p_ray)
%   6. The number of spikes in each phase bin (n_spk_phz)
%   7. The percent of spikes in each phase bin (p_spk_phz)
%   8. The preferred phase, aka the phase at the center of the bin with the most spikes (phzpref)
%
% Use as
%   spike_phase = SpikePhase(cfg, dat, spike_metrics)
% where the input values are
%   cfg                        = structure containing the configuration options
%   dat                        = structure containing the data to calculate metrics on
%   spike_metrics (optional)   = structure containing the waveform classification (broad-spiking vs.
%                                narrow spiking) for each unit, for plotting purposes only
% and the output value is
%   spike_phase                = structure containing the above listed metrics, organized into 
%                                unit x condition (x phase bin) matrices where u is index of the unit in the input dat structure
%
% The input dat structure must contain:
%   dat.condition([ACS_cond_ix]).stimwave_eq.phi_offset  = phase offset of the ACS waveform (phase at the start of stim)
%   dat.condition([ACS_cond_ix]).stimwave_eq.freq        = frequency of the ACS waveform
%   dat.condition([ACS_cond_ix]).stimwave_eq.t_offset    = time offset of the ACS waveform (-1* time in seconds at the start of stim)
%   dat.unit.condition.ts                                = the spike times for each unit, during each condition
%
% The input cfg structure should contain:
%   cfg.doplot                 = T/F (default = false), whether to plot phase histograms for
%                                each unit, each condition and a 1-figure summary of the results for the
%                                entire dataset
%   cfg.dosave                 = T/F (default = false), whether to save the spike_phase output and figures
%   cfg.doclose                = T/F (default = false), whether to close the figures automatically once they are plotted
%   cfg.n_phz_bins             = int (default = 10), the number of phase bins to use for calculating 
%                                number and percent of spikes by phase, phase preference, and Rayleigh P-value
%   cfg.nshuffs                = int (default = 1000), the number of shuffles to do to get the bootstrapped P-value
%   cfg.n_spk_persmp           = [1 x N] (default = 100), the number of spikes to include in each spike-count normalized subsample
%                                if multiple values are given, spike-count normalized values are calculated for each
%                                number of spikes
%   cfg.n_reps_subsmp          = int (default = 250), the number of subsamples of size n_spk_persmp to take for each unit when calculating
%                                the spike-count normalized coherence
%
% If cfg.doplot = true, then cfg may contain:
%   cfg.phzhist_norm           = 'count' (default) or 'probability', whether to plot the raw number of spikes or probability of spikes in
%                                each phase bin
%   cfg.sumhist_method         = 'hist' or 'line' (default) or , whether to plot the distribution of phase preferences as overlapping bars
%                                for all and significantly entrained neurons or to plot them as lines
%
% If cfg.dosave = true, then cfg should contain:
%   cfg.dat_out_fname          = string, specifying the exact location and name where to save the spike_phase output
%   cfg.plot_out_dir           = string, specifying the directory where the output plots should be saved
%
%
% By Sandon Griffin (sandon.griffin@ucsf.edu)

% Get the configuration inputs

% General Saving & Plotting Settings
dosave          = ft_getopt(cfg, 'dosave', false);
doclose         = ft_getopt(cfg, 'doclose', false);
doplot          = ft_getopt(cfg, 'doplot', true);
dat_out_fname   = ft_getopt(cfg, 'dat_out_fname', []);
plot_out_dir    = ft_getopt(cfg, 'plot_out_dir', []);

% Get function-specific settings
n_spk_persmp    = ft_getopt(cfg, 'n_spk_persmp', 100);
n_reps_subsmp 	= ft_getopt(cfg, 'n_reps_subsmp', 250);
n_phz_bins      = ft_getopt(cfg, 'n_phz_bins', 10);
nshuffs         = ft_getopt(cfg, 'nshuffs', 1000);
phzhist_norm    = ft_getopt(cfg, 'phzhist_norm', 'count');
sumplot_method  = ft_getopt(cfg, 'sumhist_method', 'line');

% Create Directory for Output Figures
if dosave && exist(plot_out_dir) ~= 7
  mkdir(plot_out_dir)
end

% Create the colormap
cmap = cbrewer('qual', 'Set1', length(dat.condition));
if length(dat.condition) == 3
  cmap = cmap([1 3 2], :);
end



% Determine the phase bin edges and centers
phz_bin_edges = linspace(0, 2*pi, n_phz_bins+1); % edges of phase bins
phz_bin_centers = phz_bin_edges(2:end) - mean(abs(diff(phz_bin_edges)))/2;



% get the wave equation values from the dat structure
c_stim = find(strcmpi({dat.condition.name}, 'stim')); % index of the condition associated with stim
phase_offset = dat.condition(c_stim).stimwave_eq.phi_offset;
freq = dat.condition(c_stim).stimwave_eq.freq;
t_offset = dat.condition(c_stim).stimwave_eq.t_offset;

%% Calculate the metrics 

% initialize arrays for storing output
coh             = nan(length(dat.unit), length(dat.condition), length(n_spk_persmp)); % Raw Coherence of all spikes
coh_Nnormspk    = nan(length(dat.unit), length(dat.condition), length(n_spk_persmp)); % Coherence normalized by a user defined number of spikes
coh_condnormspk = nan(length(dat.unit), length(dat.condition)); % Coherence normalized by the minimum number of spikes across all conditions
p_boot          = nan(length(dat.unit), length(dat.condition)); % P-value from a bootstrapped permutation test of all spikes
p_ray           = nan(length(dat.unit), length(dat.condition)); % P-value from a rayleigh test for unimodality of all spikes
n_spk_phz       = nan(length(dat.unit), length(dat.condition), length(phz_bin_centers)); % number of spikes in each phase bin 
p_spk_phz       = nan(length(dat.unit), length(dat.condition), length(phz_bin_centers)); % percentage of spikes in each phase bin
phzpref         = nan(length(dat.unit), length(dat.condition)); % preferred phase

for u = 1:length(dat.unit) % for each unit
  % initialize some stuff for this unit
  spike_phz = {};
  nspk = [];
  
  for cond = 1:length(dat.condition)
    % determine the number of spikes for this unit during this condition
    nspk(cond) = length(dat.unit(u).condition(cond).ts);
    
    % get the phase of each spike (between 0 and 2*pi)
    spike_phz{cond} = rem(phase_offset+2*pi*freq*(dat.unit(u).condition(cond).ts+t_offset), 2*pi);
    spike_phz{cond}(sign(spike_phz{cond}) == -1) = spike_phz{cond}(sign(spike_phz{cond}) == -1) + 2*pi;
    
    % Determine the nummber and percent of spikes in each phase bin
    n_spk_phz(u,cond,:) = histcounts(spike_phz{cond}, 'BinEdges', phz_bin_edges);
    if nspk(cond) > 0 % if there are any spikes in this condition during the sampled time
      p_spk_phz(u,cond,:) = (n_spk_phz(u,cond,:)./nspk(cond)).*100;
    else % if there are no spikes in this condition during the sampled time
      % set every bin to zero
      p_spk_phz(u,cond,:) = n_spk_phz(u,cond,:);
    end
    
    % Determine the preferred phase
    [~, i_pref_phz] = max(p_spk_phz(u,cond,:));
    phzpref(u,cond) = phz_bin_centers(i_pref_phz);
    
    
    if nspk(cond) > 1
      % Get the coherence and P-value for all spikes in this condition
      coh(u, cond) = get_coherence(spike_phz{cond});

      % Calculate the bootstrapped P-value
      p_boot(u,cond) = coh_bootstrap(spike_phz{cond}, nshuffs);
      
      % Calculate the Rayleigh P-value
      [p_ray(u,cond), ~]  = circ_rtest(phz_bin_edges(2:end), squeeze(n_spk_phz(u,cond,:)), mean(diff(phz_bin_edges)));
    else
      coh(u,cond) = nan;
      p_boot(u,cond) = nan;
      p_ray(u,cond) = nan;
    end
  end
  
  % determine the minimum number of spikes across all conditions for this unit
  nspk_min = min(nspk);
  for cond = 1:length(dat.condition)
    % Get coherence in pre, stim and post from repeated samples of
    % n number of spikes where n = the minimum number of spikes from
    % any given condition for this unit
    coh_condnormspk(u,cond) = get_coherence_norm(spike_phz{cond}, n_reps_subsmp, nspk_min);
    
    % Get coherence in pre, stim and post from repeated samples of
    % n number of spikes where n is defined by the user
    for n = 1:length(n_spk_persmp)
      coh_Nnormspk(u, cond, n) = get_coherence_norm(spike_phz{cond}, n_reps_subsmp, n_spk_persmp(n));
    end
  end % end condition loop
end % end unit loop

%% Plot phase histograms for each unit, each condition
if doplot
  % set the number of units per plot
  u_perplot = 7;
  nfig = 0;
  for u = 1:length(dat.unit)
    if rem(u,u_perplot) == 1 % if its time to create a new figure
      % start a new figure
      fig = figure('units','normalized','outerposition',[0 0 1 1]);
      
      % Hide the figure from being displayed
      if doclose; set(fig, 'Visible', 'off'); end
      nfig = nfig+1;
      
      % Deterine the index of the units on this figure
      if u+u_perplot < length(dat.unit)
        unit_range_fig = [u u+u_perplot];
      else
        unit_range_fig = [u length(dat.unit)];
      end
      
      % reset the unit index for this figure
      u_fig = 1;
    else
      % advance the unit index for this figure up by one
      u_fig = u_fig+1;
    end
    for cond = 1:length(dat.condition)
      subplot(3, u_perplot, (cond-1)*u_perplot+u_fig);
      hold on;

      % Plot the phase histogram for this unit, this condition
      %     polarhistogram('BinEdges', phase_bin_edges, 'BinCounts', squeeze(n_spk_phz(u,cond,:)))
      if strcmp(phzhist_norm, 'probability')
        h = histogram('BinEdges', [phz_bin_edges 2*pi+phz_bin_edges(2:end)], 'BinCounts', [squeeze(p_spk_phz(u,cond,:)); squeeze(p_spk_phz(u,cond,:))]);
      elseif strcmp(phzhist_norm, 'count')
        h = histogram('BinEdges', [phz_bin_edges 2*pi+phz_bin_edges(2:end)], 'BinCounts', [squeeze(n_spk_phz(u,cond,:)); squeeze(n_spk_phz(u,cond,:))]);
      end
      h.FaceColor = cmap(cond, :);
      h.FaceAlpha = 1;
      
      % Plot a sine wave overtop of the phase histogram
      a = gca;
      %         a.XTickLabel = [0:2:6].*57.3;
      i = 0:0.01:4*pi;
      x = 2*pi*(1/(2*pi))*i;
      plot(i, 0.2*a.YLim(2)*sin(x)+0.2*a.YLim(2), '-k', 'LineWidth', 3);
      
      a.XTick = linspace(0, 4*pi, 5);
      a.XTickLabels = round(a.XTick*57.3);
      a.XLim = [0 4*pi];
      
      % Set the title
      title_str = [dat.all.U_idstr{u} ' (' dat.unit(u).tag ')'];
      
      % Include the waveform classification in the plot if it was provided
      if exist('spike_metrics') == 1 && isfield(spike_metrics.unit, 'wf_class')
        title_str = [title_str ' : ' spike_metrics.unit(u).wf_class];
      end
      a.Title.String = title_str;
      a.Title.FontSize = 14;
      
      % Color the title in red if either of the P-values are significant 
      if p_boot(u,cond) <= 0.05 || p_ray(u,cond) <= 0.05
        a.Title.Color = [1 0 0];
      end
      
      % Add coherence and p value to the plot as a legend above the plot
      dummy(1) = plot(0, 0, '.k');
      dummy(2) = plot(0, 0, '.k');
      dummy(3) = plot(0, 0, '.k');
      leg = legend(dummy);
      leg.String = {['Coh = ' num2str(round(coh(u,cond), 3))], ['p_ray = ' num2str(round(p_ray(u,cond), 2))], ['p_boot = ' num2str(round(p_boot(u,cond), 2))]};
      leg.Location = 'NorthOutside';
      leg.FontSize = 12;
    end
    
    if rem(u,u_perplot) == 0 || u == length(dat.unit)
      if dosave
        % save the previous image
        print(fig, [plot_out_dir num2str(n_phz_bins) 'phaseBins_' num2str(unit_range_fig(1)) '_' num2str(unit_range_fig(2)) '.png'], '-dpng')
      end
      if doclose; close(fig); end
    end
  end
  
  %% Plot number of units entrained for each condition and preferred phase histograms for each condition
  fig = figure('units','normalized','outerposition',[0 0 1 1]);
  
  % Hide the figure from being displayed
  if doclose; set(fig, 'Visible', 'off'); end
  
  % Determine the index of singlle and multi-units
  i_su = find(strcmpi({dat.unit.tag}, 'SU'));
  i_mu = find(strcmpi({dat.unit.tag}, 'MU'));
  i_unit = {i_su, i_mu};
  unit_str = {'SINGLE', 'MULTI'};
  
  % plot the average and SEM coherence in each condition separately for
  % single and multi-units
  for i = 1:length(i_unit)
    subplot(2,2,(i-1)*2+1); hold on;
    dat_avg = nanmean(coh(i_unit{i}, :), 1);
    dat_sem = nanstd(coh(i_unit{i}, :), 1)./sqrt(sum(~isnan(coh(i_unit{i}, :))));
    if length(dat_avg) == length(dat_sem) && any(~isnan(dat_avg)) && any(~isnan(dat_sem))
      b1 = barwitherr(dat_sem, dat_avg, 'FaceColor', 'flat');
      for cond = 1:length(dat.condition)
        b1.CData(cond,:) = cmap(cond, :);
      end
    end
    a = gca;
    a.YLabel.String = 'Coherence';
    a.XLabel.String = 'Condition';
    a.XTick = 1:length(dat.condition);
    a.XTickLabel = {dat.condition.name};
    a.Title.String = [unit_str{i} ' Unit Coherence'];
    a.Title.FontSize = 20;
    a.FontSize = 16;
  end
  
  % plot histogram of preferred phase for all neurons (transparent or dotted),
  % and significantly coherent neurons (solid)
  for cond = 1:length(dat.condition)
    subplot(length(dat.condition), 2, (cond-1)*2+2); hold on;
    
    phzpref_keep = phzpref([i_su i_mu], :);
    
    if strcmp(sumplot_method, 'hist')
      % Make the histogram for all units transparent
      h_all = histogram(phzpref_keep(:,cond), 'BinEdges', phz_bin_edges);
      h_all.FaceAlpha = 0.1;
      h_all.FaceColor = cmap(cond,:);

      % make the histogram for significant units solid
      h_sig = histogram(phzpref_keep(p_boot([i_su i_mu],cond) <= 0.05,cond), 'BinEdges', phz_bin_edges);
      h_sig.FaceColor = cmap(cond,:);
      h_sig = histogram(phzpref_keep(p_ray([i_su i_mu],cond) <= 0.05,cond), 'BinEdges', phz_bin_edges);
      h_sig.FaceColor = cmap(cond,:);
    elseif strcmp(sumplot_method, 'line')
      h_all = histcounts(phzpref_keep(:,cond), 'BinEdges', phz_bin_edges);
      plot([phz_bin_centers 2*pi+phz_bin_centers], [h_all h_all], ':', 'Color', cmap(cond,:), 'LineWidth', 2);
      
      h_sig = histcounts(phzpref_keep(p_boot([i_su i_mu],cond) <= 0.05 | p_ray([i_su i_mu],cond) <= 0.05,cond), 'BinEdges', phz_bin_edges);
      plot([phz_bin_centers 2*pi+phz_bin_centers], [h_sig h_sig], '-', 'Color', cmap(cond,:), 'LineWidth', 3)
    end
    
    % Plot the sine wave over top
    a = gca;
    i = [0:0.01:2*pi];
    x = 2*pi*(1/(2*pi))*i;
    y = 0.125*a.YLim(2)*sin(x)+0.125*a.YLim(2);
    plot([i 2*pi+i], [y y], '-k', 'LineWidth', 2);
    
    a = gca;
    a.YLabel.String = 'No. of SU + MU';
    a.XLabel.String = 'Phase';
    title(['Preferred Phase: ' upper(dat.condition(cond).name)]);
    a.XTick = linspace(0, 4*pi, 5);
    a.XTickLabels = round(a.XTick*57.3);
    a.XLim = [0 4*pi];
    a.FontSize = 16;
  end
  
  if dosave
    print(fig, [plot_out_dir 'Summary_' num2str(n_phz_bins) 'phaseBins.png'], '-dpng')
    if doclose; close; end
  end
end

%% Save Data
if dosave
  spike_phase = [];
  spike_phase.cfg = cfg;
  spike_phase.info.U_id = dat.all.U_id;
  spike_phase.info.U_idstr = dat.all.U_idstr;
  spike_phase.coh = coh;
  spike_phase.coh_condnormspk = coh_condnormspk;
  for n = 1:length(n_spk_persmp)
    eval(['spike_phase.coh_Nnormspk' num2str(n_spk_persmp(n)) ' = squeeze(coh_Nnormspk(:, :, n));']);
  end
  spike_phase.pval_bs = p_boot;
  spike_phase.pval_rl = p_ray;
  spike_phase.n_spk_phz = n_spk_phz;
  spike_phase.p_spk_phz = p_spk_phz;
  spike_phase.phzpref = phzpref;
  spike_phase.phzbin_edges = phz_bin_edges;
  spike_phase.phzbin_centers = phz_bin_centers;
  
  save(dat_out_fname, 'spike_phase', '-v7.3');
end
