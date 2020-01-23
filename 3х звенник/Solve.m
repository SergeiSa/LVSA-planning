    function [g, theta, d_theta] = Solve(tau, q, n, w_max, dt, g_desired)
        cvx_begin
        variables g(n) theta(n) d_theta(n-1)
        minimize norm(g_desired - g)
        subject to
        0 <= g <= 1;
        (theta - q) - tau.*g == 0;
        for index = 1:(n-1)
            d_theta(index) == (theta(index+1) - theta(index)) / dt;
        end
        -w_max <= tau(1:(n-1)).*d_theta <= w_max;
        cvx_end
    end