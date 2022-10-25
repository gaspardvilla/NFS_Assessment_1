function d = get_d(N, M_u, A_u)
    % Test function for the computation of d
    
    % Initialization
    n = N - 1;
    d = zeros(n, 1);
    
    % For loop
    for i = 1 : n
        d(i) = A_u(i) / M_u(i,i);
    end

end