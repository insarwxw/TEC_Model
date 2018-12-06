function du=TEC_source(volgeom,xloc,nu,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computes surface displacements due to thermal elastic contraction (TEC)
%%% of a heat source, such as an intruded magma body by a volcano eurption.
%%%  
%%% Oct. 12,2017, Xiaowen Wang, The University of Tokyo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%     volgeom = TMC source geometry: East, North, Depth, Volume
%               <length, length, length, volume>
%        units:    km      km      km      km^3
%
%        xloc = Matrix of grid coordinates <length>, stored columnwise,
%               i.e., [east, north], units: m
%          nu = Poisson's ratio
%          dt = elapsed time after the magma intrusion (year) 
%Outputs:
%          du = surface displacements in East, North, and Vertical directions
%               at the time "dt" after the magma intrusion 
%
%          U_disp: store the displacement at different time if dt is a vector;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other parameters related to the source 
% a:        thermal expansity 
% k:        thermal diffusivity (m^2/s)
% T:        magma temprature (K) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read input variables
x0=volgeom(1)*1e3;y0=volgeom(2)*1e3;z0=volgeom(3)*1e3; %source location (m)
V=volgeom(4)*1e9; % Convert km^3 to meter (m^3)
x=xloc(:,1); y=xloc(:,2); %coordinates of surface grid points (m)

%%%%%%%%%
a=2e-5;  
k=1e-5;
T=1200;
r=nthroot(3*V/(4*pi),3); %radius of the source volume 
Q=V*T; 
q=Q/V;           %heat source volume density
%%%%%%%%% 
dz=0.01;         %represent surface depth (meter)
z=ones(size(x)).*dz;

dt_len=length(dt);
Le=[];Ln=[];Lu=[];

for i=1:dt_len

  t=dt(i);
  name_t=['t_',num2str(i)];

  %%%function Thermal_disp_Volume
  [u0,e0_13,e0_23,e0_33]=TEC_fun_Volume(x0,y0,z0,x,y,z,t,a,nu,k,r,q);
  [u1,e1_13,e1_23,e1_33]=TEC_fun_Volume(x0,y0,z0,x,y,-z,t,a,nu,k,r,q);

  %%%Equation （1） of Furuya (2005)
  f1=2*z(:)*ones(1,3);
  f2=[e1_13,e1_23,-e1_33];
  u=u0+(3-4*nu).*u1-f1.*f2;
  %%%
  ux=u(:,1);
  uy=u(:,2);
  uz=u(:,3);
   
  U_disp.(name_t).ux=ux';
  U_disp.(name_t).uy=uy';
  U_disp.(name_t).uz=uz';
%%%
  Le=[Le;ux'];
  Ln=[Ln;uy'];
  Lu=[Lu;uz'];

end
  
save U_disp U_disp
du=[Le;Ln;Lu];

%%%Calculate the displacement velocity 
%B=dt';vx=nan(1,size(x,1));vy=nan(1,size(x,1));vz=nan(1,size(x,1));

%for j=1:length(ux)
%  le=Le(:,j);pe=polyfit(B,le,1);vx(j)=pe(1);
%  ln=Ln(:,j);pn=polyfit(B,ln,1);vy(j)=pn(1);
%  lu=Lu(:,j);pu=polyfit(B,lu,1);vz(j)=pu(1);
%end

%vu=[vx;vy;vz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

