%% [D, mu] = maxent(mu_hat, sigma, lambda, Nbin)%% Compute Chebyshev moments density of states by MEMfunction [D,mu] = maxent(mu_hat,sigma,lambda,Nbin)  % calcualte Jackson kernal
  kpm = filter_jackson(mu_hat);  %% mu_hat is a normolized matrix  % Normalize input data  %  scale_factor = mu_hat(1);  sigma  = sigma  / scale_factor;  mu_hat = mu_hat / scale_factor;  % Initialization
  %
  M = length(mu_hat);
  K = 4;               % Extrapolation 4 <= K <= 10
  I = 4;               % For accuracy improvement 4 <= I
  L = M*K*I;           % Total pixels

  % Apply filter for extrapolated sum -- then work with modified moments
  %
  mu = zeros(M*K,1);  mu(1:M) = mu_hat;  mu = filter_jackson(mu);  mu_hat(1:M) = mu(1:M);  mu(1:M) = filter_jackson(mu(1:M));  %% what the hell are you doing ......

  chi2    = sum( (((mu_hat-mu(1:M))./sigma)).^2 )
  alpha   = chi2;
  lambda0 = (mu(1:M)-mu_hat)./(alpha*(sigma.^2));
  xi      = mu_hat - mu(1:M) + alpha*(sigma.^2).*lambda0;
  d0      = dnought(mu,L,M);

  step = 1/2;

  % termination criterion, according to paper
  %  Outer loop: cond1 (mentioned above eq(21))
  %  Middle loop: cond2 (init. false)
  cond1 = (chi2 <= M);


  fprintf('Enter main loop\n');
  while not(cond1)

    fprintf('Run with alpha = %e\n', alpha);
    % Save old values in case we need to reduce alpha step and go back

    % Newton iteration
    %
    while 1

      % Take a Newton step
      xi = mu_hat - mu(1:M) + alpha*(sigma.^2).*lambda0;
      fprintf('Newton step (rnorm: %e)\n', norm(xi));
      H = zeros(M);
      for r=1:M
        for c=1:M
          H(r,c) = ( mu(r+c-1) + mu(abs(r-c)+1) )/2;
        end
        H(r,r) = H(r,r) + alpha*(sigma(r))^2;
      end
      lambda0 = lambda0-H\xi;

      d = dnought2d(d0,L,M,lambda0);  % eq(23) function at the bottom. Get D with new lambda.
      mu = d2mu(d,L,M,K);  % New mu derived from eq(18). Function at the bottom.

      %This line calculates relative entropy eq(19) with L discrete samples by the trapezoidal rule.
      S = sum( (d - d0 - d.*log(d./d0))*(pi/(length(d)-1)) );

      chi2 = sum( (((mu_hat-mu(1:M))./sigma)).^2 );  % Calculate chi_square
      Qp = S - chi2/(2*alpha);
      Qd = Qp + sum( (xi.^2)./(2*alpha*(sigma.^2)) );

      fprintf('  -> Qd = %e, Qp = %e, Chi2 = %e\n', Qd, Qp,chi2);
      % Termination criterion (below (28))
      if (abs((Qd-Qp)/Qd)<=0.02)	% Within 2% diffrence
        fprintf('Newton iteration terminated\n');
        break;
      end

    end
    chi2_old    = chi2;
    xi_old      = xi;
    lambda0_old = lambda0;
    mu_old      = mu;
    alpha_old = alpha;


    cond1 = (chi2 <= M);  % termination criterion

      % Retrieve old values if new chi_square is larger than previous one, reduce alpha step
    if (chi2 > chi2_old*10)
      fprintf('Chi Big\n');
     	step = step/2;  % Reduce alpha step
     	alpha = alpha_old*(1-step);    % step size alpha*(1- 1/2 1/4 ...1/2n)
     	chi2 = chi2_old;
    	xi = xi_old;
    	lambda0 = lambda0_old;
	   	mu = mu_old;
    else

      % Saved values %
      alpha = alpha*(1-step);
    end


  end

  D = d;  %  Fetch final D(phi) output
  compare_cheb(1-lambda,kpm,(D*scale_factor*200),Nbin,M,K);
  figure(2)
  hold on
  semilogy(D*scale_factor*200)
  hold off



  %  SUBFUNCTIONS BELOW   %

  %  Get density of states function from mu
  function d0 = dnought(mu,L,M)
    d0 = zeros(L+1,1);
    for i = 1:(L+1)
      phi = pi*((i-1)+1/2)/(L+1);  % match L+1 sampling
      d0(i)= 1 + 2*sum( mu(2:M).*cos(phi*(1:M-1)') );  % eq(8)
    end
  end

  %  Get new D(phi) from D0(phi) and lambda
  function d = dnought2d(d0,L,M,lambda0)
    d = zeros(L+1,1);
    for i = 1:(L+1)
      phi = pi*((i-1)+1/2)/(L+1);  % match L+1 sampling
      d(i)= d0(i)*exp(-sum( lambda0.*cos(phi*(0:M-1)') )); % eq(23)
    end
    fprintf('# of Negative D =  %d\n', sum(d<0));
  end

  %  Transform D(phi) to mu
  function mu = d2mu(d,L,M,K)
    mu = zeros(M*K,1); % M*K point extrapolation for mu
    phil = pi*((0:L)'+1/2)/(L+1); % L+1 point sampling in eq(18)
    for i = 1:(M*K)
      mu(i) = sum( cos((i-1)*phil).*d ); % eq(18)
    end
  end

end
