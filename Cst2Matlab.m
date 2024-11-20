
clc
clear 
close all

%% =========== cst invokes =======================================================================
cst = actxserver('CSTStudio.application');
mws=invoke(cst,'Openfile','E:\n_omrani\C_shape\Data2\FullStructure.cst');
invoke(mws, 'Rebuild' );
units=invoke(mws,'Units');
invoke(units,'Geometry','um');
invoke(units,'Frequency','THz');
invoke(units,'Voltage','V');
invoke(units,'Resistance','Ohm');
invoke(units,'Inductance','H');
invoke(units,'TemperatureUnit','Kelvin');
invoke(units,'Time','ns');
invoke(units,'Current','A');
invoke(units,'Conductance','Siemens');
invoke(units,'Capacitance','F');
solver=invoke(mws,'Solver');
invoke(solver,'FrequencyRange',2,6);
mesh=invoke(mws, 'Mesh' );
invoke(mesh,'MeshType','PBA');
invoke(mesh, 'SetCreator', 'High Frequency');
    
%% ========== parameters =================================================================

nn=2^2;% 2bit phase
mm=2^3;% 3bit amplitude
dx=25;% unit cell dimension in um
dy=dx;
A=load('E:\n_omrani\C_shape\Data2\AC.mat');
Ph=load('E:\n_omrani\C_shape\Data2\Phimn.mat');
[Nx,Ny]=size(A.AC);

%% ========== calculations ===============================================================

for n=1:Nx
    for m=1:Ny
        if A.AC(m,n)==0.089625
            if Ph.phimn(m,n)==0
            i=1;
            elseif Ph.phimn(m,n)==pi/2
            i=2;
            elseif Ph.phimn(m,n)==pi
            i=3;
            else 
            i=4;
            end
        elseif A.AC(m,n)==0.17925
            if Ph.phimn(m,n)==0
            i=5;
            elseif Ph.phimn(m,n)==pi/2
            i=6;
            elseif Ph.phimn(m,n)==pi
            i=7;
            else 
            i=8;
            end
        elseif A.AC(m,n)==0.268875
            if Ph.phimn(m,n)==0
            i=9;
            elseif Ph.phimn(m,n)==pi/2
            i=10;
            elseif Ph.phimn(m,n)==pi
            i=11;
            else 
            i=12;
            end
        elseif A.AC(m,n)==0.3585
            if Ph.phimn(m,n)==0
            i=13;
            elseif Ph.phimn(m,n)==pi/2
            i=14;
            elseif Ph.phimn(m,n)==pi
            i=15;
            else 
            i=16;
            end
        elseif A.AC(m,n)==0.448125
            if Ph.phimn(m,n)==0
            i=17;
            elseif Ph.phimn(m,n)==pi/2
            i=18;
            elseif Ph.phimn(m,n)==pi
            i=19;
            else 
            i=20;
            end
        elseif A.AC(m,n)==0.53775
            if Ph.phimn(m,n)==0
            i=21;
            elseif Ph.phimn(m,n)==pi/2
            i=22;
            elseif Ph.phimn(m,n)==pi
            i=23;
            else 
            i=24;
            end

        elseif A.AC(m,n)==0.627375
            if Ph.phimn(m,n)==0
            i=25;
            elseif Ph.phimn(m,n)==pi/2
            i=26;
            elseif Ph.phimn(m,n)==pi
            i=27;
            else 
            i=28;
            end

        elseif A.AC(m,n)==0.717

            if Ph.phimn(m,n)==0
            i=29;
            elseif Ph.phimn(m,n)==pi/2
            i=30;
            elseif Ph.phimn(m,n)==pi/2
            i=31;
            else 
            i=32;
            end

        end
        if m*n~=1
            transform = invoke(mws, 'Transform');
            component=invoke(mws,'Component');
            invoke(transform, 'Reset');
            invoke(transform, 'Name', ['component',num2str(i)]);
            invoke(transform, 'Vector', num2str((m-1)*dx),num2str((n-1)*(-dy)),'0');
            invoke(transform, 'UsePickedPoints','False');
            invoke(transform, 'InvertPickedPoints','False');
            invoke(transform, 'MultipleObjects','True');
            invoke(transform, 'GroupObjects', 'False');
            invoke(transform, 'Repetitions', '1');
            invoke(transform, 'MultipleSelection','False');
            invoke(component,'New',['component',num2str(i),num2str(m),'-',num2str(n)]);
            invoke(transform,'Destination',['component',num2str(i),num2str(m),'-',num2str(n)]);
            invoke(transform,'Material','');
            invoke(transform, 'Transform','Shape','Translate');
            release(transform)
        end
    end
end
for j=2:mm*nn
   component=invoke(mws,'Component');
   invoke(component,'Delete', ['component',num2str(j)]);
end


