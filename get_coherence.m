function coherence = get_coherence(phz)
% Get coherence given phase

% Delete nans
phz(isnan(phz)) = [];

xs = cos(phz);
ys = sin(phz);

x_mean = mean(xs);
y_mean = mean(ys);

coherence = (x_mean^2 + y_mean^2)^.5;

% 
% % phz = normrnd(pi/2, pi/2, 1, 100); % high coherence example
% phz = 2*pi*rand(1, 100); % low coherence example
% 
% % get the phase of each spike
% phz(sign(phz) == -1) = phz(sign(phz) == -1) + 2*pi;
% 
% figure; 
% subplot(1, 2, 1);
% hold on;
% cfg = [];
% cfg.plot_sine = 'continuous';
% cfg.repeat_x = true;
% cfg.color = [0.3020, 0.6863, 0.2902];
% cfg.plot_legend = false;
% plot_hist(cfg, phz);
% 
% ax = gca;
% ax.FontSize = 30;
% ax.YLim = [0 20];
% 
% subplot(1, 2, 2);
% hold on;
% plot(xs, ys, 'o', 'MarkerEdgeColor', [0.3020, 0.6863, 0.2902], 'MarkerSize', 8);
% axis square;
% ax = gca;
% ax.YLim = [-1 1];
% plot(0, 0, 'or', 'MarkerSize', 8, 'LineWidth', 2);
% plot(x_mean, y_mean, 'x', 'Color', [0.3020, 0.6863, 0.2902], 'MarkerSize', 10, 'LineWidth', 3);
% 
% ax = gca;
% ax.FontSize = 30;
% 
