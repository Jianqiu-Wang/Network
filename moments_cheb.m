% c = moments_cheb(A, N, num_samples, kind)
%
% Estimate Chebyshev moments c_0, ... c_{N-1} where
%   c_n = sum_i T_n(lambda_i) = tr(T_n(A))
% using an estimator based on num_samples Gaussian vectors.

function [c,sigma,num_samples] = moments_cheb(A, N, num_samples, kind)

  if nargin < 4, kind = 1;          end
  if nargin < 3, num_samples = 100; end
  if nargin < 2, N = 10;            end

  m = length(A);
  Z = randn(m,num_samples);
  N
  c = zeros(N,1);
  sigma = zeros(N,1);
  dev = zeros(num_samples,1);

  c(1) = m;
  c(2) = trace(A);
  sigma(1) = c(1) * 1e-10;
  % sigma(2) = c(2) * 1e-10 + 1e-16;

  P0 = Z;
  P1 = kind*(A*Z);
  for j = 1:num_samples
      dev(j) = Z(:,j)'*P1(:,j);
  end
  sigma(2) = std(dev);
    
  for np = 3:N
    Pn = 2*(A*P1) - P0;
    for j = 1:num_samples
      dev(j) = Z(:,j)'*Pn(:,j);
    end
    c(np) = mean(dev);
    sigma(np) = std(dev);
    P0 = P1;
    P1 = Pn;
  end
