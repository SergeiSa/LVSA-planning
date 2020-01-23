
dt = 0.01;
n = 50;

q_max = 1;
q_min = 0;
tf = dt * n;

M = [1 0 0 0;
     1 tf tf^2 tf^3;
     0 1 0 0
     0 1 2*tf 3*tf^2];
a = M \ [q_min; q_max; 0; 0];

q = zeros(n, 1);
for i = 1:n
    t = (i-1)*dt;
    q(i) = dot(a, [1 t t^2 t^3]);
end


k_desired = [5; 5; 20; 100; 100];
g_desired = 1./k_desired;

tau = [2; 4; 6; 8; 10];
q = [0.1; 0.2; 0.3; 0.4; 0.5];

n = length(k_desired);

w_max = 500;

cvx_begin

variables g(n) theta(n) d_theta(n-1) 



minimize norm(g_desired - g)
subject to
    0 <= g <= 1;
    (theta - q) - tau.*g == 0;
    
    for i = 1:(n-1)
        d_theta(i) == (theta(i+1) - theta(i)) / dt;
    end
    
    tau(1:(n-1)).*d_theta <= w_max;
cvx_end

k = 1./g;
disp(['k: ', mat2str(k)])