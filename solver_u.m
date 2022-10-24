function [u_star, M_u, b_u] = solver_u(N, A_u, A_p, rho, u_old, p_star, p_0)
    % This function solve the linear system M_u * u = b_u associated with
    % the discretized momentum equation.

    % Remark: it's not the same N than the nb of integration points, it's
    % the size of the sample.

    % Initialization of M_u and b_u
    M_u = zeros(N, N);
    b_u = zeros(N, 1);

    % 1 case
    F_A = rho * u_old(1) * A_u(1);
    F_B = rho * A_p(2) * ((u_old(1) + u_old(2)) / 2);
    M_u(1, 1) = F_B + ((1/2) * F_A * ((A_u(1) / A_p(1))^2));
    b_u(1) = ((p_0 - p_star(2)) * A_u(1)) + (F_A * u_old(1) * (A_u(1) / A_p(1)));

    % N case
    M_u(N, N) = rho * u_old(N) * A_u(N);
    M_u(N, N-1) = -rho * A_p(N-1) * ((u_old(N-1) + u_old(N)) / 2);
    b_u(N) = (p_star(N-1) - p_star(N)) * A_u(N);

    % 2 to N-1 cases
    for i = 2 : N-1
        M_u(i,i) = rho * A_p(i+1) * ((u_old(i) + u_old(i+1)) / 2);
        M_u(i, i-1) = -rho * A_p(i) * ((u_old(i-1) + u_old(i)) / 2);
        b_u(i) = (p_star(i) - p_star(i+1)) * A_u(i);
    end

    % Solve the linear system
    u_star = linsolve(M_u, b_u);

end