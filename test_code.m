%%%%%%%%%% Test 1: calculated surface displacements at a specified time after magma intrusion
clc;clear
[xx,yy] = meshgrid(-100:100,-100:100); %horizontal coordinates 
xloc=[xx(:),yy(:)];
volgeom=[0,0,0.1,1e-3]; % a heat source at depth of 100 m with volume of 1e-3 km^3
nu=0.25;                %Poisson's ratio
dt=5;                   %five years after the heat source emplacement  

du=TEC_source(volgeom,xloc,nu,dt);

%%%plot the displacment 
ux=reshape(du(1,:),size(xx));
uy=reshape(du(2,:),size(xx));
uz=reshape(du(3,:),size(xx));

figure;imagesc(ux);title('Disp in E-W direction (m)');colorbar;
figure;imagesc(uy);title('Disp in N-S direction (m)');colorbar;
figure;imagesc(uz);title('Disp in U-D direction (m)');colorbar;


%%%%%%%%%% Test 2: Calculated surface displacements at different time squences
clc;clear

[xx,yy] = meshgrid(-100:100,-100:100); %horizontal coordinates 
xloc=[xx(:),yy(:)];
volgeom=[0,0,0.1,1e-3]; % a heat source at depth of 100 m with volume of 1e-3 km^3
nu=0.25;                %Poisson's ratio
dt=5;                   %five years after the heat source emplacement  

dt=1:1:10;              %ten years after the heat source emplacement  
du=TEC_source(volgeom,xloc,nu,dt);


%%%plot temporal evolution of vertical displacments at the center point  
load U_disp;
line=101;col=101; %center pixel 

tz=[];
dt_len=length(dt);
for id=1:dt_len
    name_t=['t_',num2str(id)];

%   dx=reshape(U_disp.(name_t).ux,size(xx));
%   dy=reshape(U_disp.(name_t).uy,size(xx));
    dz=reshape(U_disp.(name_t).uz,size(xx));
    tz(id,:)=dz(line,col);  
    time=dt(id);
    fprintf('Displacement at the time: %f(year)\n', time)
end 

figure;plot(dt,tz,'r-diamond','LineWidth',4);
xlabel('Time (year)','FontName','Times New Roman','FontSize',18);
ylabel('Vertical displacements','FontName','Times New Roman','FontSize',18);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Temporal evolution of vertical displacement')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END
