close all
clc
clear
load("WindFarmCS_CheckNewCoordinates.mat")
load("WindFarmWF_CheckNewCoordinates.mat")
Turbine=1;
heightpoint=50;
CrossPoint=50;

%Wind Speed Plot
z=zeros(1,Z_Res);
u=zeros(1,Z_Res);

for i=1:Z_Res
u(i)=u_ambient_zprofile(i);
z(i)=Z_Levels(i);
end

plot(u,z)

surf(XCoordinates(:,:,1,1),YCoordinates(:,:,1,1),Delta_U(:,:,heightpoint,1),"LineStyle","none")
colormap(flipud(inferno))
