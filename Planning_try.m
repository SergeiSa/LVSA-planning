
k_desired = [5; 5; 20; 100; 100];
g_desired = 1./k_desired;

tau = [2; 4; 6; 8; 10];
q = [0.1; 0.2; 0.3; 0.4; 0.5];

n = length(k_desired);

w_max = 100;
dt = 0.001;

cvx_begin

variables g(n) theta(n) 

% for i = 1:(n-1)
%     d_theta(i) = (theta(i+1) - theta(i)) / dt;
% end

minimize norm(g_desired - g)
subject to
    0 <= g <= 1;
    (theta - q) - tau.*g == 0;
%     tau.*d_theta <= w_max;
cvx_end

k = 1./g;
disp(['k: ', mat2str(k)])