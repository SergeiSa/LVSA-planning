function byuserStep6_InverseKinematics()
close all;

%Create user interfase object for SRD
SRD = SRDuserinterface;

%Get initial position
InitialPosition = SRD.GetInitialPosition();
IC_Task = g_InverseKinematics_Task(InitialPosition);

ZeroOrderDerivativeNodes = {IC_Task(1), 0.1;
                            IC_Task(2), IC_Task(2)+0.30;  
                            IC_Task(3), 0}; 
FirstOrderDerivativeNodes = {0, 0; 
                             0, 0; 
                             0, 0}; 
SecondOrderDerivativeNodes = {0, 0; 
                              0, 0; 
                              0, 0}; 

TimeOfOneStage = 2;
TimeEnd = (size(ZeroOrderDerivativeNodes, 2) - 1)*TimeOfOneStage;
NodeTimes = (0:TimeOfOneStage:TimeEnd)';
Spline = TPSplineConstructorUI();
Spline.OutOfBoundariesBehaviour = 'LastValue';
Spline.GenerateSplines(NodeTimes, ZeroOrderDerivativeNodes, FirstOrderDerivativeNodes, SecondOrderDerivativeNodes);
Spline.Plot;
 
%%%%%%%%

%Solve inverse kinematics numerically
SRD.SetupNumericInverseKinematics('DesiredTask', @Spline.EvaluateAll, ...
    'TimeRange', [0; TimeEnd], ...
    'TimeStep', 0.01,...
    'PolynomialDegree', 6, ... 
    'NumberOfSegments', 20, ...
    'SolverType', 'lsqnonlin', ...
    'ToPlot', true, ...
    'Verbose', true);
%you can use 'lsqnonlin', 'fmincon' and 'ProjectionNonlinSolver'
end