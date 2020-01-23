close all; clear; %clear classes; %clc;

%Create user interfase object for SRD
SRD = SRDuserinterface;
%Create ground link
Ground = SRD.CreateGroundLink();

%Create other links
index = 1;
Shin =  SRDLinkWithJoint('abs_pivotY', index, 'datafile_Shin',  Ground, 1); index = index + 1;
Hip =   SRDLinkWithJoint('abs_pivotY', index, 'datafile_Hip',   Shin,   1); index = index + 1;
Torso = SRDLinkWithJoint('abs_pivotY', index, 'datafile_Torso', Hip,    1); index = index + 1;

%Set generalised coordinates; Use negative numbers to make clock-wise
%rotations be positive
Shin.SetUsedGenCoordinates(1);
Hip.SetUsedGenCoordinates(2);
Torso.SetUsedGenCoordinates(3);

LinkArray = [Ground; Shin; Hip; Torso]; %Create array of links
InitialPosition = [pi/7; -1.3*pi/3; 1*pi/5]; %Define initial position of the robot
% InitialPosition = [pi/6; -pi/6; pi/6]; %Define initial position of the robot
AxisLimits = [-1; 1; -1; 1; 0; 2]; %Set axis limits for the snapshot
ViewAngle = [0, 0];

SRD.CreateRobotStructure(LinkArray, InitialPosition, AxisLimits, ViewAngle); %Create the robot