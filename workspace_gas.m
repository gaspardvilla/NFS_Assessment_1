%% Clean
close all
clearvars
clc


%% Not moving part

% Initialization of some parameters
L = 2.0;
rho = 1.0;
A_in = 0.5;
A_out = 0.1;

% Stopping conditions 
tol = 10^(-6);
nb_iter_max = 10^4;

% Linear dependence with x of the area
A = @(x) A_in + ((A_out - A_in) * x / L);

% Pressure
p_0 = 10;
p_out = 0;
dp = p_0 - p_out;


%% What we can change

% Initialization
N = 21; % Nb of nodes
alpha_u = 0.1;
alpha_p = 0.1;

r_u = tol + 1;
r_p = tol + 1;

x_p = linspace(0, L, N);
x_u = x_p(1:N-1) + (L / (2 * (N-1)));

% Area definition
A_p = A(x_p);
A_u = A(x_u);

% Take a guess
m_dot = 1;


%% Hard part

% Initialization of the velocity and the pressure
u_old = m_dot ./ (rho * A_u);
p_star = p_0 - ((dp * x_p) / L);

nb_iter = 0;



%% For Loop is coming

while (r_u > tol) && (r_p > tol) && (nb_iter < nb_iter_max)

    % Solvers
    [u_star, M_u, b_u] = solver_u(N-1, A_u, A_p, rho, u_old, p_star, p_0);
    [p_prime, b_p] = solver_p(N, A_u, A_p, rho, u_old, u_star);
    
    % Correctors
    d = get_d(N, rho, u_old, A_u, A_p);
    u_calc = u_star;
    p_calc = p_prime;
    
    u_calc(1) = u_star(1) + (d(1) * (p_prime(1) - p_prime(2)));
    p_calc(1) = p_0 - ((1/2) * rho *((u_calc(1) * (A_u(1) / A_p(1)))^2));
    
    for i = 2 : N - 1
        u_calc(i) = u_star(i) + (d(i) * (p_prime(i) - p_prime(i+1)));
        p_calc(i) = p_star(i) + p_prime(i);
    end
    
    p_calc(N) = p_star(N) + p_prime(N);
    
    % Residual and RHS
    r_u = norm((M_u * u_calc) - b_u) / norm(diag(M_u) .* u_calc);
    r_p = norm(b_p);
    
    % Update wothg under-relaxation
    u_old = (alpha_u * u_calc) + ((1 - alpha_u) * u_old);
    p_star = (alpha_p * p_calc) + ((1 - alpha_p) * p_star);
    
    % Update iterators
    nb_iter = nb_iter + 1;
    if mod(nb_iter, 100) == 0
        fprintf('Ieration: %4.2f\n', nb_iter);
    end

end


%% True solutions

true_u = @(x) (A_out ./ A(x)) * sqrt(2 * p_0 / rho);
true_p = @(x) p_0 * (1 - ((A_out ./ A(x)).^2));




%% Plots
figure(1);
plot(x_u, true_u(x_u));
hold on
plot(x_u, u_old);
hold off


figure(2);
plot(x_p, true_p(x_p));
hold on
plot(x_p, p_star);
hold off











