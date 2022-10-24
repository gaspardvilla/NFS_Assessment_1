function [p_prime, b_p] = solver_p(N, A_u, A_p, rho, u_old, u_star)
    % This function solve the linear system M_p * p = b_p associated with
    % the discretized momentum equation.

    % Initialization of M_u and b_u
    M_p = zeros(N, N);
    b_p = zeros(N, 1);

    % 1 and N cases
    M_p(1,1) = 1;
    M_p(N,N) = 1;
    b_p(1) = 0;
    b_p(N) = 0;

    % Compute the vector d
    d = get_d(N, rho, u_old, A_u, A_p);

    % 2 to N-1 cases
    for i = 2 : N-1
        M_p(i,i-1) = -rho * A_u(i-1) * d(i-1);
        M_p(i,i+1) = -rho * A_u(i) * d(i);
        M_p(i,i) = -(M_p(i,i-1) + M_p(i,i+1));
        b_p(i) = (rho * A_u(i-1) * u_star(i-1)) - (rho * A_u(i) * u_star(i));
    end

    % Solve the linear system
    p_prime = linsolve(M_p, b_p);

end