close all
clc
clear
load("WindFarmCS_CheckNewCoordinates.mat")
load("WindFarmWF_CheckNewCoordinates.mat")
Turbine=1;
heightpoint=100;
CrossPoint=100;
AxialPoint=50;

%Wind Speed Plot
z=zeros(1,Z_Res);
u=zeros(1,Z_Res);

for i=1:Z_Res
u(i)=u_ambient_zprofile(i);
z(i)=Z_Levels(i)*D/H;
end

figure(1)
hold on
plot(u,z)
plot(linspace(min(u),u_ambient,2),[1,1],'r')
plot([u_ambient,u_ambient],linspace(0,1,2),'r')
title("Height/ Shear profile of the wind.")
xlabel("Wind speed in m/s")
ylabel("Height in 1/H")
xlim([0,max(u)])

%Crosssection Graph (In Y/X plane)
figure(2)
surf(XCoordinates(:,:,1,1),YCoordinates(:,:,1,1),Delta_U(:,:,heightpoint,1),"LineStyle","none")
colormap(flipud(inferno))
title(["Crossection in Y/X plane at Z=",num2str(z(heightpoint)*D/H),"H"])
view([0,0,2])
xlabel("x/D")
ylabel("y/D")
col=colorbar;
ylabel(col,'Delta u')

figure(3)
surf(XCoordinates(:,:,1,1),YCoordinates(:,:,1,1),u_ambient_zprofile(heightpoint)-Delta_U(:,:,heightpoint,1),"LineStyle","none")
colormap(inferno)
title(["Crossection in Y/X plane at Z=",num2str(z(heightpoint)*D/H),"H. Resulting speed."])
view([0,0,2])
xlabel("x/D")
ylabel("y/D")
col=colorbar;
ylabel(col,'u_{wake}')
clim([0, u_ambient]);

%Height Graphs (In Crosssection)
z=repmat(z,N,1);
x=XCoordinates(:,:,1);
for i=1:Z_Res
    u_u(:,i)=Delta_U(:,CrossPoint,i,1);
    if u_ambient_zprofile(i)>0
    u_tot(:,i)=u_ambient_zprofile(i)-u_u(:,i);
    else
    u_tot(:,i)=zeros(N,1);   
    end
end
figure(4)
surf(x,z,u_u,"LineStyle","none")
colormap(flipud(inferno))
title(["Crossection in Z/X plane at Y=",num2str(YCoordinates(1,CrossPoint,1,1)),"H. Speed deficit."])
view([0,0,2])
xlabel("x/D")
ylabel("z/D")
col=colorbar;
ylabel(col,'Delta u')
clim([0, u_ambient]);

figure(5)
surf(x,z,u_tot,"LineStyle","none")
colormap(inferno)
title(["Crossection in Z/X plane at Y=",num2str(YCoordinates(1,CrossPoint,1,1)),"H. Resulting speed."])
view([0,0,2])
xlabel("x/D")
ylabel("z/D")
col=colorbar;
ylabel(col,'u_{wake}')
clim([0, u_ambient]);

figure(6)
plot(u_tot(AxialPoint,:),z(AxialPoint,:))
title(["Wind velocity in Crossection in Z/X plane at Y=",num2str(YCoordinates(1,CrossPoint,1,1)),"H. At axial distance from turbine",num2str(x(AxialPoint,1)),"."])
xlabel("u_{wake}")
ylabel("z/D")
xlim([0,max(u_tot(AxialPoint,:))])
ylim([0,max(z(AxialPoint,:))])
