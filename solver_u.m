function [u_star, M_u, b_u] = solver_u(N, u_old, A_u, A_p, p_0, p_star, rho)
    % Test function of the solver for u

    % Initialization
    n = N - 1;
    M_u = zeros(n, n);
    b_u = zeros(n, 1);

    % 1 case
    F_A = rho * u_old(1) * A_u(1);
    F_B = rho * A_p(2) * ((u_old(1) + u_old(2)) / 2);
    M_u(1, 1) = F_B + ((1/2) * F_A * ((A_u(1) / A_p(1))^2));
    b_u(1) = ((p_0 - p_star(2)) * A_u(1)) + (F_A * (A_u(1) / A_p(1)) * u_old(1));

    % n case
    M_u(n, n) = rho * u_old(n) * A_u(n);
    M_u(n, n-1) = (-rho) * ((u_old(n-1) + u_old(n)) / 2) * A_p(n);
    b_u(n) = (p_star(n) - p_star(n+1)) * A_u(n);

    % 2 to n-1 cases
    for i = 2 : n-1
        M_u(i, i) = rho * ((u_old(i) + u_old(i+1)) / 2) * A_p(i + 1);
        M_u(i, i-1) = (-rho) * ((u_old(i-1) + u_old(i)) / 2) * A_p(i);
        b_u(i) = (p_star(i) - p_star(i+1)) * A_u(i);
    end

    % Compute the final result
    u_star = M_u \ b_u;
end