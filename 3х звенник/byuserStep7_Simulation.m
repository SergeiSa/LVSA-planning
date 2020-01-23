close all;

%Create user interface object for SRD
SRD = SRDuserinterface;

%Load SimulationEngine and set up the simulation parameters
SimulationEngine = SRD.GetSimulationEngine();
SimulationEngine.CustomSolverType = 'Taylor';
%Can use 'Euler', 'Taylor', 'Runge', 'Implicit Euler', 'DAE Taylor', 'DAE Runge';
SimulationEngine.IC.q = SimulationEngine.IC.q + rand(3, 1)*0.5;
SimulationEngine.IC.v = SimulationEngine.IC.v + rand(3, 1)*0.5;


%Load InverseKinematicsEngine
InverseKinematicsEngine = SRD.GetInverseKinematicsEngine();

SimulationEngine.Time = InverseKinematicsEngine.TimeEnd - 0.0;

ControlInput = @InverseKinematicsEngine.EvaluatePolynomialApproximation;
% ControlInput = SimulationEngine.GetPlugs('Constant_IC_ControlInput');

%Set up a controller

%PD controller example
% Controller = SimulationEngine.GetPDcontroller('Varying gains PD', 'Kp', eye(SimulationEngine.dof)*500, ...
%                                                                   'Kd', eye(SimulationEngine.dof)*100);
%Can use .GetPDcontroller with 'PD', 'Varying gains PD', 
%'Computed torque PD' 'Computed torque PD wpinv'
%
%or .GetLQRcontroller with 'LQR', 'ProjectedLQR', 'LQR wpinv'
%or .GetMPcontroller with 'MP', 'One step MP'

%LQR example
Controller = SimulationEngine.GetLQRcontroller('LQR', 'unified_Q', 10000, 'unified_R', 1, ...
    'ILQR_TimeStep', 0.01);

%LQR example
% Controller = SimulationEngine.GetMPcontroller('MP', 'unified_Q', 10000, 'unified_R', 1);
% Controller = SimulationEngine.GetMPcontroller('MP', 'unified_Q', 10000, 'unified_R', 1, ...
%     'NumberOfPredictionSteps', 10, 'MP_PredictionTimeStep', 0.005);

% LimitController = SimulationEngine.LimitControllerOutput(LQRController, 150, -180);

%Simulate
tic
Res = SimulationEngine.Simulation(ControlInput, Controller);
%Can use .Simulation() and .SimulationStateSpace()
toc

FontSize = 20;


figure_handle = figure('Color', 'w');

subplot(2, 2, 1)
plot(Res.SimulationOutput.Time, Res.SimulationOutput.Position(:, 1), '-', 'LineWidth', 1.2, 'Color', [0.75 0.1 0]); hold on;
plot(Res.SimulationOutput.Time, Res.SimulationOutput.Position(:, 2), '--', 'LineWidth', 1.5, 'Color', [0.2 0.05 0.2]); hold on;
plot(Res.SimulationOutput.Time, Res.SimulationOutput.Position(:, 3), ':', 'LineWidth', 2, 'Color', [0.0 0.5 0]); hold on;
legend({'$q_1$', '$q_2$', '$q_3$'}, 'Interpreter', 'latex', 'FontSize', FontSize)
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
ylabel_handle = ylabel('$$q_i$$, (rad)');
ylabel_handle.Interpreter = 'latex';


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%



w_max = 100;

tau = Res.SimulationOutput.ControlActions(1:40:end, 2);
q = Res.SimulationOutput.Position(1:40:end, 2);
time = Res.SimulationOutput.Time(1:40:end);

n = size(tau, 1);
dt = time(2) - time(1);

%%%%%%%%%%%%%%%%%%

tf = time(end);

k_max = 2500;
k_min = 500;

M_k = [1 0 0 0;
       1 tf tf^2 tf^3;
       0 1 0 0
       0 1 2*tf 3*tf^2];
a_k = M_k \ [k_min; k_max; 0; 0];

k_desired = zeros(n, 1);
g_desired = zeros(n, 1);
for i = 1:n
    t            = time(i);
    k_desired(i) = dot(a_k, [1 t t^2 t^3]);
    g_desired(i) = 1 / k_desired(i);
end

%%%%%%%%%%%%%%%%%


cvx_begin

variables g(n) theta(n) d_theta(n-1) 

minimize norm(g_desired - g)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%



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
plot(time, k, ':', 'LineWidth', 2.5, 'Color', [0 0.4 0.4]); hold on;
plot(time, k_desired, '--', 'LineWidth', 1.5, 'Color', [1 0 0]); 
legend({'$k$', '$k^*$'}, 'Interpreter', 'latex', 'FontSize', FontSize)
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
ylabel_handle = ylabel('$$k$$, (Nm/rad)');
ylabel_handle.Interpreter = 'latex';

subplot(2, 2, 3)
plot(time(1:(n-1)), w, '-', 'LineWidth', 1, 'Color', [0.8 0.2 0]); hold on;
plot(time, w_max*ones(size(time)), '--', 'LineWidth', 2, 'Color', [0.0 0.5 0]); 
plot(time, -w_max*ones(size(time)), '-.', 'LineWidth', 1.5, 'Color', [0.0 0 0.5]);  
legend({"$w$", "$w_{max}$", "$-w_{max}$"}, 'Interpreter', 'latex', 'FontSize', FontSize)
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
ylabel_handle = ylabel('$$w$$, W');
ylabel_handle.Interpreter = 'latex';

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
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%Plot the output
figure_handle = SimulationEngine.PlotSimulationResults(Res.SimulationOutput, 'P, dP; V, U');

%If need - animate the resulting motion
% ToAnimate = true;
% if ToAnimate
%     Animation = SRDAnimation();
%     Animation.Animation_Accelerator = 100;
%     Animation.Animate(Res.SimulationOutput.Position);
% end