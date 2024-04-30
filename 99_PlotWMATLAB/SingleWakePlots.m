close all
clc
clear
load("WindFarmCS.mat")
load("WindFarmWF.mat")
%% Velocity deficit
%X-Y Direction
figure
plot(YCoordinates(2,:,1), c_0_vec(1)-Delta_U(2,:,1), LineStyle="none",Marker=".")
xlabel("Y/D")
ylabel("DELTA u/ u_a")
title("Velocity deficit for X=7D")

%X-Z Direction
figure
plot(c_0_vec(1)-Delta_U(2,:,1), ZCoordinates(2,:,1).*D./H, LineStyle="none",Marker=".")
xlabel("Z/D")
ylabel("DELTA u/ u_a")
title("Velocity deficit for X=7D-- XZ")

%Y-Z Direction
figure
plot3(YCoordinates(2,:,1), ZCoordinates(2,:,1), c_0_vec(1)-Delta_U(2,:,1), LineStyle="none",Marker=".")
xlabel("Y/D")
ylabel("Z/D")
zlabel("DELTA u/ u_a")
title("Velocity deficit for X=7D in Y and Z direction")

%Contour Y-Z Direction
%Reshape values
[Y,Z]=meshgrid(linspace(-2.98,3,1000), linspace(-2.1,3.85,1000));
% Interpolate DATA onto the grid using linear interpolation
DATA_grid = griddata(YCoordinates(2,:,1), ZCoordinates(2,:,1), c_0_vec(1)-Delta_U(2,:,1), Y, Z, 'linear');
% DATA_grid(isnan(DATA_grid))=0;

figure
surf(Y, Z, DATA_grid, EdgeColor="none")
colormap inferno;
xlabel("Y/D")
ylabel("Z/D")
zlabel("DELTA u/ u_a")
title("Velocity deficit for X=7D in Y and Z direction")
colorbar

%% Turbulence
%X-Y Direction
figure
plot(YCoordinates(2,:,1), Delta_TI(2,:,1), LineStyle="none",Marker=".")
xlabel("Y/D")
ylabel("TI+")
title("Rotor-added Turbulence for X=7D-- XY")

%Y-Z Direction
figure
plot3(YCoordinates(2,:,1), ZCoordinates(2,:,1), Delta_TI(2,:,1), LineStyle="none",Marker=".")
xlabel("Y/D")
ylabel("Z/D")
zlabel("TI+")
title("Rotor-added Turbulence for X=7D-- YZ")


%X-Z Direction
figure
plot(Delta_TI(2,:,1), ZCoordinates(2,:,1), LineStyle="none",Marker=".")
ylabel("Z/H")
xlabel("TI+")
title("Rotor-added Turbulence for X=7D-- XZ")

%Contour Y-Z Direction
%Reshape values
[Y,Z]=meshgrid(linspace(-2.98,3,1000), linspace(-2.1,3.85,1000));
% Interpolate DATA onto the grid using linear interpolation
DATA_grid = griddata(YCoordinates(2,:,1), ZCoordinates(2,:,1), Delta_TI(2,:,1), Y, Z, 'linear');
% DATA_grid(isnan(DATA_grid))=0;

figure
surf(Y, Z, DATA_grid, EdgeColor="none")
ylim([0,4])
colormap inferno;
xlabel("Y/D")
ylabel("Z/D")
zlabel("DELTA u/ u_a")
title("Rotor-added turbulence for X=7D in Y and Z direction")
colorbar
zlim([-0.1, 0.3])
x=1;