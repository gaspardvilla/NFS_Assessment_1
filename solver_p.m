function [p_prime, M_p, b_p] = solver_p(N, rho, M_u, A_u, u_star)
    % Test function of the solver for p

    % Initialization
    M_p = zeros(N, N);
    b_p = zeros(N, 1);
    d = test_d(N, M_u, A_u);

    % 1 and N cases
    M_p(1, 1) = 1;
    M_p(N, N) = 1;

    b_p(1) = 0;
    b_p(N) = 0;
    
    % 2 to N-1 cases
    for i = 2 : N-1
        M_p(i, i-1) = (-rho) * d(i-1) * A_u(i-1);
        M_p(i, i+1) = (-rho) * d(i) * A_u(i);
        M_p(i, i) = -(M_p(i, i-1) + M_p(i, i+1));
        b_p(i) = (rho * A_u(i-1) * u_star(i-1)) - (rho * A_u(i) * u_star(i));
    end

    % Compute the final result
    p_prime = M_p \ b_p;
end