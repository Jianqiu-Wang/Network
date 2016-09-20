% compare_cheb(lambda, c, Nbin)
%
% Compare a histogram of the eigenvalue list lambda with a scaled
% Chebyshev density plot.

function compare_cheb(lambda, c, D, Nbin, M, K)

  if nargin < 3, Nbin = 100; end
  h = 2/Nbin;
  xx = linspace(-1,1,1000);
  xx = xx(2:end-1);
  length(D)
  xx_m = linspace(-1,1,length(D));
  xx_m = xx_m(M*K+1:end-M*K); % Cut both ends
  yy = plot_cheb(c*h, xx);
  hist(lambda, Nbin);
  hold on; 
  plot(xx, yy, 'r-', 'linewidth', 2);
  plot(xx_m, D(end-M*K:-1:M*K+1), 'g-', 'linewidth', 1); % Flip D
  legend('Histogram','KPM','MEM');
  hold off
  axis tight
