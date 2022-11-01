%% Clean
close all
clearvars
clear all
clc


%% Not moving part

% Display the plots
plot_figure = False;

% Initialization of some parameters
L = 2.0;
rho = 1.0;
A_in = 0.5;
A_out = 0.1;

% Stopping conditions 
tol = 10^(-6);
nb_iter_max = 10^6;

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

r_u = [tol + 1];
r_p = [tol + 1];

x_p = linspace(0, L, N)';
x_u = x_p(1:N-1) + (L / (2 * (N-1)));

% Area definition
A_p = A(x_p);
A_u = A(x_u);

% Take a guess
m_dot = 1;


%% Hard part

% Initialization of the velocity and the pressure
u_old = (m_dot ./ (rho * A_u));
p_star = (p_0 - ((dp * x_p) / L));

nb_iter = 0;


%% True solutions
true_u = @(x) (A_out ./ A(x)) * sqrt(2 * p_0 / rho);
true_p = @(x) p_0 * (1 - ((A_out ./ A(x)).^2));



%% For Loop is coming

while ((r_u(end) > tol) || (r_p(end) > tol)) && (nb_iter < nb_iter_max)

    % Solvers
    [u_star, M_u, b_u] = solver_u(N, u_old, A_u, A_p, p_0, p_star, rho);
    [p_prime, M_p, b_p] = solver_p(N, rho, M_u, A_u, u_star);


    % Correctors
    d = get_d(N, M_u, A_u);
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
    r_u(end + 1) = norm((M_u * u_calc) - b_u) / norm(diag(M_u) .* u_calc);
    r_p(end + 1) = norm(b_p);
    
    % Update wothg under-relaxation
    u_old = (alpha_u * u_calc) + ((1 - alpha_u) * u_old);
    p_star = (alpha_p * p_calc) + ((1 - alpha_p) * p_star);
    
    % Update iterators
    nb_iter = nb_iter + 1;
    if mod(nb_iter, 100) == 0
        fprintf('Iteration: %4.f\n', nb_iter);
    end

end


%% Final computations

% Flow rate
m_flow = @(x) rho .* x .* A_u;
flow_rate_pred = m_flow(u_old);
flow_rate_true = m_flow(true_u(x_u));

% Global mass balance
mass = abs(flow_rate_pred(1) - flow_rate_pred(end)) / mean(flow_rate_pred);

% Relative errors
rel_err = @(u, v) norm(u - v) / norm(u);
vel_error = rel_err(true_u(x_u), u_old);
pressure_error = rel_err(true_p(x_p), p_star);
flow_error = rel_err(flow_rate_true, flow_rate_pred);


%% Final prints
fprintf('Final r_u = %4.10f\n', r_u(end));
fprintf('Final r_p = %4.10f\n', r_p(end));
fprintf('Relative error for velocity = %4.10f\n', vel_error);
fprintf('Relative error for pressure = %4.10f\n', pressure_error);
fprintf('Relative error for flow rate = %4.10f\n', flow_error);
fprintf('Global mass balance = %4.25f\n', mass);


%% Plots

if plot_figure

    % Residuals plot
    figure(1);
    ru = semilogy((1:1:nb_iter), r_u(2:end));
    hold on;
    grid on;
    rp = semilogy((1:1:nb_iter), r_p(2:end));
    yline(tol, 'r--', 'Stopping condition', 'LabelHorizontalAlignment', 'left');
    last_iter = strcat('Last iteration:', {' '}, num2str(nb_iter));
    xline(nb_iter, 'k', last_iter);
    xlabel('Iteration');
    ylabel('Residual value');
    legend([ru, rp], {'r_u', 'r_p'})
    hold off;
    
    % Numerical solution against analytical solution
    figure(2);
    analytical_u = plot(x_u, true_u(x_u));
    hold on;
    grid on;
    numerical_u = plot(x_u, u_old, 'o-');
    xlabel('Position [m]');
    ylabel('Velocity value [m/s]');
    legend([analytical_u, numerical_u], {'Analytical solution u', 'Numerical solution u_{final}'})
    hold off;
    
    
    % Numerical solution against analytical solution
    figure(3);
    analytical_p = plot(x_p, true_p(x_p));
    hold on;
    grid on;
    numerical_p = plot(x_p, p_star, 'o-');
    xlabel('Position [m]');
    ylabel('Pressure value [Pa]');
    legend([analytical_p, numerical_p], {'Analytical solution p', 'Numerical solution p_{final}'})
    hold off;
    
    % Trial
    figure(4);
    yyaxis left;
    analytical_u = plot(x_u, true_u(x_u), '--');
    hold on;
    grid on;
    numerical_u = plot(x_u, u_old, 'o-');
    ylabel('Velocity value [m/s]');
    yyaxis right;
    analytical_p = plot(x_p, true_p(x_p), '--');
    numerical_p = plot(x_p, p_star, 'o-');
    ylabel('Pressure value [Pa]');
    xlabel('Position [m]');
    legend([analytical_u, numerical_u, analytical_p, numerical_p], {'Analytical solution u', 'Numerical solution u_{final}', 'Analytical solution p', 'Numerical solution p_{final}'})
    hold off;

end











