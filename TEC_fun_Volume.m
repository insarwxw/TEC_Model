function [u,e_13,e_23,e_33]=TEC_fun_Volume(x0,y0,z0,x,y,z,t,a,v,k,r,q)
%%%%%%%%%
% x0,y0,z0 
% x,y,z:    coordinates of the surface grid points (meters)
% d:        depth of the heat source (m)
% t:        elapsed time after the magma intrusion  
% a:        thermal expansity 
% v:        poisson's ratio
% k:        thermal diffusivity (m^2/s)
% T0:       magma temprature 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%initialize the coordinate system 
C0=[x0,y0,-z0]; %% coordinats of the source point 
C=[x,y,z];      %% coordinates of the grid points 
len=length(C);

%%%distance between the source center and the grid points 
dR=C-ones(len,1)*C0;
R=sqrt(dR(:,1).^2+dR(:,2).^2+dR(:,3).^2);

%%%
m=(a*(1+v))/(1-v);
%%%convert the time elapse from year to seconds 
ts=t*365*86400;
s=4*k*ts;

X_l=(R-r)./sqrt(s); %lowrer limit of the integration
X_h=(R+r)./sqrt(s); %upper limit of the integration 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pha_R=nan(len,1);pha_RR=nan(len,1);
u=nan(len,3);e_33=nan(len,1);e_13=nan(len,1);e_23=nan(len,1);
for i=1:len
    i_l=X_l(i); %lowrer limit of the integration
    i_h=X_h(i); %upper limit of the integration 
    i_R=R(i); 

%%%%%fun_1: the integration part of the equation （16） 
    fun_1 = @(X) (X.^2./2-i_R.*X./sqrt(s)+i_R.^2./s).*(-X.*erfc(X)+exp(-X.^2)./sqrt(pi))+(X.^3.*erfc(X))./6-(exp(-X.^2).*(X.^2+1))./(6*sqrt(pi));
%%%%%fun_2: the integration part of the equation （17）
    fun_2 = @(X) ((X-i_R./sqrt(s)).^2.*(-X.*erfc(X)+exp(-X.^2)./sqrt(pi)))./2+(X.^3./6-i_R.^2.*X./(2*s)+i_R.^3./(2*sqrt(s.^3))).*erfc(X)-(exp(-X.^2).*(X.^2+1))./(6*sqrt(pi));
    
    v_fun1=fun_1(i_h)-fun_1(i_l);
    v_fun2=fun_2(i_h)-fun_1(i_l);

    pha_R(i)=(q*m*r^3)./(3*i_R.^2)+((q*m*sqrt(s.^3))./(2*i_R.^2)).*v_fun1; %equation (16)
    pha_RR(i)=-(2*q*m*r^3)./(3*i_R.^3)-((q*m*(sqrt(s.^3)))./i_R.^3).*v_fun2; %equation (17)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Calculate the displacements based on the equation (12)

    u(i,:)=(dR(i,:).*pha_R(i))./i_R;

%%%%Calculate the strain e_ii;e_ij based on the equation (13-14)

    dR_33=dR(i,3);
    e_33(i)=((1-(dR_33/i_R).^2)/i_R).*pha_R(i)+((dR_33/i_R).^2).*pha_RR(i);

    dR_13=dR(i,1).*dR(i,3);
    e_13(i)=(dR_13/i_R.^2).*(pha_RR(i)+pha_R(i)/i_R);

    dR_23=dR(i,2).*dR(i,3);
    e_23(i)=(dR_23./i_R.^2).*(pha_RR(i)+pha_R(i)/i_R);

end
%%%%%
end


