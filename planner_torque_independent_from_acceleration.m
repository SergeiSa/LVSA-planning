close all; clear; clc;

dt = 0.01;
n = 50;

tf = dt * n;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

q_max = 1;
q_min = 0;

M_q = [1 0 0 0;
       1 tf tf^2 tf^3;
       0 1 0 0
       0 1 2*tf 3*tf^2];
a_q = M_q \ [q_min; q_max; 0; 0];


k_max = 100;
k_min = 10;

M_k = [1 0 0 0;
       1 tf tf^2 tf^3;
       0 1 0 0
       0 1 2*tf 3*tf^2];
a_k = M_k \ [k_min; k_max; 0; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%

mass = 2;
w_max = 300;
%%%%%%%%%%%%%%%%%%%%%%%%%%

q = zeros(n, 1);
k_desired = zeros(n, 1);
g_desired = zeros(n, 1);
tau = zeros(n, 1);
tau1 = zeros(n, 1);
for i = 1:n
    t            = (i-1)*dt;
    
    time(i)      = t;
    q(i)         = dot(a_q, [1 t t^2 t^3]);
    k_desired(i) = dot(a_k, [1 t t^2 t^3]);
    g_desired(i) = 1 / k_desired(i);
    
    %q = a1 + a2*t + a3*t*t + a4*t*t*t
    %dq/dt = a2 + 2*a3*t + 3*a4*t*t
    %ddq/ddt = 2*a3 + 6*a4*t
    ddq = 2*a_q(3) + 6*a_q(4) * t;
    tau(i) = ddq*mass + 50*exp(-100*(t - tf/2)^2);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%

weight = 0;


cvx_begin

variables g(n) theta(n) d_theta(n-1) 

minimize norm(g_desired - g) + weight*norm(theta(1) - q(1))
subject to
    0 <= g <= 1;
    (theta - q) - tau.*g == 0;
    
    for i = 1:(n-1)
        d_theta(i) == (theta(i+1) - theta(i)) / dt;
    end
    
    -w_max <= tau(1:(n-1)).*d_theta <= w_max;
    
cvx_end

k = 1./g;
%disp(['k: ', mat2str(k)])
w = tau(1:(n-1)).*d_theta;

inertia = 0.05;
dd_theta = diff(d_theta) / dt;
u_theta = dd_theta*inertia + tau(1:length(dd_theta));

inertia_k = 0.1;
damping_k = 0.1;
dk = diff(k) / dt;
u_k = inertia_k*dk + damping_k*k(1:length(dk));


FontSize = 20;

figure_handle = figure('Color', 'w');

subplot(2, 2, 1)
plot(time(1:length(u_theta)), u_theta, '-', 'LineWidth', 1, 'Color', 'k'); hold on;
plot(time(1:length(u_k)), u_k, '--', 'LineWidth', 2, 'Color', 'r'); hold on;
plot(time, tau, '-.', 'LineWidth', 1.5, 'Color', 'b'); hold on;
legend({'$u_{\theta}$', '$u_{k}$', '$\tau$'}, 'Interpreter', 'latex', 'FontSize', FontSize)
grid on; grid minor;
ax = gca;
ax.GridAlpha = 0.6;
ax.LineWidth = 0.5;
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.2;
ax.FontName = 'Times New Roman';
ax.FontSize = FontSize;
xlabel_handle = xlabel('$$t$$, s');
xlabel_handle.Interpreter = 'latex';
ylabel_handle = ylabel('$$u$$, $$\tau$$');
ylabel_handle.Interpreter = 'latex';


figure_handle = figure('Color', 'w');

subplot(2, 2, 1)
plot(time, q, '--', 'LineWidth', 2, 'Color', 'r'); hold on;
plot(time, theta, 'LineWidth', 1.2, 'Color', 'k')
legend({'$q$', '$\theta$'}, 'Interpreter', 'latex', 'FontSize', FontSize)
grid on; grid minor;
ax = gca;
ax.GridAlpha = 0.6;
ax.LineWidth = 0.5;
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.2;
ax.FontName = 'Times New Roman';
ax.FontSize = FontSize;
xlabel_handle = xlabel('$$t$$, s');
xlabel_handle.Interpreter = 'latex';
ylabel_handle = ylabel('$$q$$, $$\theta$$, (rad)');
ylabel_handle.Interpreter = 'latex';



subplot(2, 2, 2)
plot(time, k); hold on;
plot(time, k_desired, '--', 'LineWidth', 2); 
legend({'$k$', '$k^*$'}, 'Interpreter', 'latex', 'FontSize', FontSize)
grid on; grid minor;
ax = gca;
ax.GridAlpha = 0.6;
ax.LineWidth = 0.5;
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.2;
ax.FontName = 'Times New Roman';
ax.FontSize = FontSize;

subplot(2, 2, 3)
plot(time(1:(n-1)), w, '-', 'LineWidth', 1, 'Color', 'k'); hold on;
plot(time, w_max*ones(size(time)), '--', 'LineWidth', 2, 'Color', 'r'); 
plot(time, -w_max*ones(size(time)), '-.', 'LineWidth', 1.5, 'Color', 'b');  
legend({"$w$", "$w_{max}$", "$-w_{max}$"}, 'Interpreter', 'latex', 'FontSize', FontSize)
grid on; grid minor;
ax = gca;
ax.GridAlpha = 0.6;
ax.LineWidth = 0.5;
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.2;
ax.FontName = 'Times New Roman';
ax.FontSize = FontSize;

subplot(2, 2, 4)
plot(time, tau); hold on;
legend({'$tau$'}, 'Interpreter', 'latex', 'FontSize', FontSize)
grid on; grid minor;
ax = gca;
ax.GridAlpha = 0.6;
ax.LineWidth = 0.5;
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.2;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

