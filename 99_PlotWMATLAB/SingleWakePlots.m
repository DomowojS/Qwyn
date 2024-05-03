close all
clc
clear
load("WindFarmCS.mat")
load("WindFarmWF.mat")
%% Create Grid & Interpolate
Y=linspace(-1.49,1.5,1000);
Z=linspace(1.2, 4.1, 1000);
utmp=Delta_U(2,:,1);

ytmp=YCoordinates(2,:,1);
[~,id,~]=unique(YCoordinates(2,:,1));
u_interpY=interp1(ytmp(id),utmp(id),Y);

ztmp=ZCoordinates(2,:,1);
[~,id,~]=unique(ZCoordinates(2,:,1));
u_interpZ=interp1(ztmp(id),utmp(id),Z);


%% Velocity deficit
%Y Direction
figure
plot(Y, u_0_vec(1)-u_interpY, LineStyle="none",Marker=".")
xlabel("Y/D")
ylabel("DELTA u/ u_a")
title("Velocity deficit for X=7D")

%X-Z Direction
id=ZCoordinates(2,:,1)*D==H;
u=u_0_vec(1)-Delta_U(2,:,1);
y=YCoordinates(2,:,1);
figure
plot(u(id), y(id).*D)
ylabel("Z/D")
xlabel("DELTA u/ u_a")
xlim([0,8])
ylim([-120,120])
title("Velocity deficit for X=7D-- XZ")

%Y-Z Direction
figure
plot3(YCoordinates(2,:,1), ZCoordinates(2,:,1), u_ambient_zprofile(2,:,1)-Delta_U(2,:,1), LineStyle="none",Marker=".")
xlabel("Y/D")
ylabel("Z/D")
zlabel("DELTA u/ u_a")
title("Velocity deficit for X=7D in Y and Z direction")

%Contour Y-Z Direction
%Reshape values
[Y,Z]=meshgrid(linspace(-2.98,3,1000), linspace(-2.1,3.85,1000));
% Interpolate DATA onto the grid using linear interpolation
DATA_grid = griddata(YCoordinates(2,:,1), ZCoordinates(2,:,1), u_ambient_zprofile(2,:,1)-Delta_U(2,:,1), Y, Z, 'linear');
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