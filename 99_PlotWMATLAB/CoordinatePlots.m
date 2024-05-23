close all
clc
clear
load("WindFarmCS_NewCoordinateSys.mat")
load("WindFarmWF_NewCoordinateSys.mat")
Turbine=1;

%2D Plot (X/Y)
hold on
plot(XCoordinates(1,:,1,Turbine),YCoordinates(1,:,1,Turbine),'LineStyle','none', 'Marker','.')
plot(XCoordinates(2,:,1,Turbine),YCoordinates(2,:,1,Turbine),'LineStyle','none', 'Marker','.')
plot(XCoordinates(3,:,1,Turbine),YCoordinates(3,:,1,Turbine),'LineStyle','none', 'Marker','.')
xlim([-3,17])
ylim([-2,5])
title("2D-Plot of Coordinates")

%3D Plot (X/Y/Z)
figure(2)
hold on
for i=1:numel(Z_Levels)
z=zeros(1,Y_Res)+Z_Levels(i);

plot3(XCoordinates(1,:,1,Turbine),YCoordinates(1,:,1,Turbine),z./80,'LineStyle','none', 'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
plot3(XCoordinates(2,:,1,Turbine),YCoordinates(2,:,1,Turbine),z./80,'LineStyle','none', 'Marker','.','MarkerFaceColor','g','MarkerEdgeColor','g')
plot3(XCoordinates(3,:,1,Turbine),YCoordinates(3,:,1,Turbine),z./80,'LineStyle','none', 'Marker','.','MarkerFaceColor','r','MarkerEdgeColor','r')
xlim([-3,17])
ylim([-2,5])
zlim([0,7])
end