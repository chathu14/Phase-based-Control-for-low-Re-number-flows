function growth_rate = growth_rate_func(t,y)

% Assuming you have time (t) and signal (y)
t = t(:); y = y(:);  % ensure column vectors

% Take the absolute value to get the envelope
env = abs(y);

% Fit exponential: env ≈ A * exp(σ * t)
fitFunc = @(b, t) b(1) * exp(b(2) * t);  % b(2) is growth rate σ
initialGuess = [max(env), 0.1];  % [A, sigma]

% Use nonlinear least squares fit
opts = optimset('Display','off');
params = lsqcurvefit(fitFunc, initialGuess, t, env, [], [], opts);

growth_rate = params(2);