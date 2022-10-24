function d = get_d(N, rho, u_old, A_u, A_p)
    % Compute the vector d
    d = zeros(N-1, 1);

    F_A = rho * u_old(1) * A_u(1);
    F_B = rho * A_p(2) * ((u_old(1) + u_old(2)) / 2);
    a_1 = F_B + ((1/2) * F_A * ((A_u(1) / A_p(1))^2));
    d(1) = A_u(1) / a_1;

    a_n = rho * u_old(N-1) * A_u(N-1);
    d(N-1) = A_u(N-1) / a_n;

    for i = 2 : N-2
        a_e = rho * A_p(i+1) * ((u_old(i) + u_old(i+1)) / 2);
        d(i) = A_u(i) / a_e;
    end
end