RelativeBase = [0; 0; 0];        
RelativeFollower = [0; 0; 0.5];   
RelativeCoM = [0; 0; 0.25]; 
Mass = 0.05*2 + 2*2;
Inertia = eye(3);
Name = 'Shin';
save('datafile_Shin', 'RelativeBase', 'RelativeFollower', 'RelativeCoM', 'Mass', 'Inertia', 'Name');

RelativeBase = [0; 0; 0];        
RelativeFollower = [0; 0; 0.5];   
RelativeCoM = [0; 0; 0.25];  
Mass = 0.23 * 50 + 2*3;
Inertia = eye(3);
Name = 'Hip';
save('datafile_Hip', 'RelativeBase', 'RelativeFollower', 'RelativeCoM', 'Mass', 'Inertia', 'Name');

RelativeBase = [0; 0; 0];        
RelativeFollower = [0; 0; 0.7];   
RelativeCoM = [0; 0; 0.35]; 
Mass = 0.55 * 50 + 5;
Inertia = eye(3);
Name = 'Torso';
save('datafile_Torso', 'RelativeBase', 'RelativeFollower', 'RelativeCoM', 'Mass', 'Inertia', 'Name');


