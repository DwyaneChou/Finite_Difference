clc
clear

g  = 9.80616;
dx = 0.01;
dt = 0.0001;
history_interval = 0.2; %seconds
run_time         = 40.0000001;%seconds
% run_time         = 0.800001;%seconds

min_z = 9980;
max_z = 10020;

x  = -180:dx:180;
x  = x';
nx = size(x,1);

% Initialize the fields
% z = ones(size(x))*12;
% z(abs(x)<=1)=18;
% u = zeros(size(x));

[u,z]=Gaussian(x,g);

e = 0.5*(u.^2.*z + z.^2);

output_count = 0;

z_sum(1)  = sum(z);
e_sum(1)  = sum(e);

figure('visible','off')
plot(x,z)
xlim([min(x) max(x)])
ylim([min_z max_z])
disp(num2str(output_count))
disp(['z = ',num2str(z_sum(1))])
disp(['e = ',num2str(e_sum(1))])
title(['Finite Difference CRK4 step ',num2str(output_count,'%05d')])
print(gcf,'-r300','-dpng',['.\output\MCV',num2str(output_count,'%05d'),'.png']);

output(output_count,u,z,e,dx,dt,history_interval)

% Integration
nt     = run_time/dt;
for i = 1:nt    
    [uTend1,zTend1] = tend(dx,u,z);
    
    u1 = u + 0.5*dt*uTend1;
    z1 = z + 0.5*dt*zTend1;
    
    [uTend2,zTend2] = tend(dx,u1,z1);
    
    u2 = u + 0.5*dt*uTend2;
    z2 = z + 0.5*dt*zTend2;
    
    [uTend3,zTend3] = tend(dx,u2,z2);
    
    u3 = u + dt*uTend3;
    z3 = z + dt*zTend3;
    
    [uTend4,zTend4] = tend(dx,u3,z3);
    
    du  = (uTend1 +2*uTend2 +2*uTend3 +uTend4)/6;
    dz  = (zTend1 +2*zTend2 +2*zTend3 +zTend4)/6;
        
    E      = 0.5*u.^2 + z;
    
    E3     = 0.5*u3.^2 + z3;
    
    an  = -innerProduct(du.^2,dz)*dt^2;
    bn  = (innerProduct(du.^2,z)...
         + 2*innerProduct(u.*du,dz)...
         + innerProduct(dz,dz))*dt;
    cn  = -2*(innerProduct(du,z.*u)...
            + innerProduct(dz,E) ...
            - innerProduct(E3,zTend4)...
            - innerProduct(z3.*u3,uTend4));
       
    if bn>=0
        beta_n = 2*cn/(bn+sqrt(bn^2-4*an*cn));
    else
        beta_n = 2*cn/(bn-sqrt(bn^2-4*an*cn));
    end
    
    tau_n = beta_n*dt;
%     tau_n = dt;
    
    u  = u + tau_n*du;
%     u  = u - Diff6(u,nx);
    z  = z + tau_n*dz;
%     z  = z - Diff6(z,nx);

    e  = 0.5*(u.^2.*z + z.^2);
    
    z_sum(i+1)    = sum(z);
    e_sum(i+1)    = sum(e);
    
    disp(['tau_n = ',num2str(tau_n,'%5.5e'),...
          ' de = '  ,num2str((e_sum(i+1) - e_sum(1))/sum(e_sum(1)),'%5.5e'),...
          ' dz = '  ,num2str((z_sum(i+1) - z_sum(1))/sum(z_sum(1)),'%5.5e')])
    
    if rem(i,history_interval/dt)==0
        figure('visible','off')
        plot(x,z)
        xlim([min(x) max(x)])
        ylim([min_z max_z])
        output_count = output_count + 1;
        disp(num2str(output_count))
        title(['Finite Difference CRK4 step ',num2str(output_count,'%05d')])
        print(gcf,'-r300','-dpng',['.\output\MCV',num2str(output_count,'%05d'),'.png']);
        output(output_count,u,z,e,dx,dt,history_interval)
    end
end

function [uTend,zTend] = tend(dx,u,z)
% u wind speed
% z geopotential height
% e energy

nx = size(u,1);

E  = 0.5*u.^2+z;
zu = z.*u;

% Compute cell tend
fu = E;
fz = zu;

% 6th order
Lu = CD6(fu,dx,nx);
Az = CD6(fz,dx,nx);

% Update u tend
uTend = -Lu;

% Update z tend
zTend = -Az;

% Energy Conservation check
% EnergyConservation = innerProduct(E,zTend)...
%                    + innerProduct(zu,uTend);
% disp(['EnergyConservation = ',num2str(EnergyConservation)])

end

function ab=innerProduct(a,b)
ab = sum(a.*b);
end

function output(output_count,u,z,e,dx,dt,history_interval)
f_out = 'output.nc';
output_precision = 'NC_DOUBLE';
west_east = size(u,1);

if output_count==0
    mode           = netcdf.getConstant('NETCDF4');
    mode           = bitor(mode,netcdf.getConstant('CLOBBER'));
    ncid           = netcdf.create(f_out,mode);
    disp(['ncid = ',num2str(ncid)])
        
    % Define Dimensions
    time_dimID             = netcdf.defDim(ncid,'time'            ,netcdf.getConstant('NC_UNLIMITED'));
    west_east_dimID        = netcdf.defDim(ncid,'west_east'       ,west_east);
    
    % Define Attribute
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Model_Source'     ,'MCV_1d written by Zhou Lilong')
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'dx'               ,dx)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'dt'               ,dt)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history_interval' ,history_interval)
    
    % Define Variables
    u_id = netcdf.defVar(ncid,'u',output_precision,[west_east_dimID,time_dimID]);
    netcdf.putAtt(ncid,u_id,'units','m/s');
    netcdf.putAtt(ncid,u_id,'description','wind speed on the left point');
    
    z_id = netcdf.defVar(ncid,'z',output_precision,[west_east_dimID,time_dimID]);
    netcdf.putAtt(ncid,z_id,'units','m^2/s^2');
    netcdf.putAtt(ncid,z_id,'description','geopotential height on the left point');
    
    e_id = netcdf.defVar(ncid,'e',output_precision,[west_east_dimID,time_dimID]);
    netcdf.putAtt(ncid,e_id,'units','m^4/s^4');
    netcdf.putAtt(ncid,e_id,'description','energy on the left point');
    
    % Put Variables
    netcdf.putVar(ncid, u_id    ,[0,0],[west_east,1],u   );
    netcdf.putVar(ncid, z_id    ,[0,0],[west_east,1],z   );
    netcdf.putVar(ncid, e_id    ,[0,0],[west_east,1],e   );
    
    netcdf.close(ncid)
else
    ncid = netcdf.open(f_out,'WRITE');
        
    u_id    = netcdf.inqVarID(ncid,'u');
    z_id    = netcdf.inqVarID(ncid,'z');
    e_id    = netcdf.inqVarID(ncid,'e');
    netcdf.reDef(ncid)
    
    netcdf.putVar(ncid, u_id,[0,output_count],[west_east,1],u);
    netcdf.putVar(ncid, z_id,[0,output_count],[west_east,1],z);
    netcdf.putVar(ncid, e_id,[0,output_count],[west_east,1],e);
    
    netcdf.close(ncid)
end

end

function [u,z]=Gaussian(lambda,g)
mu    = 0;
sigma = 20/sqrt(2*pi);

u = zeros(size(lambda));
z = 10000 - g + exp(-(lambda-mu).^2/(2*sigma^2))*g;
end

% 6th order center difference operator
function DF = CD6(F,dx,nx)
DF(4:nx-3,1) = (-F(1:nx-6) + 9*F(2:nx-5) - 45*F(3:nx-4) + 45*F(5:nx-2) - 9*F(6:nx-1) + F(7:nx))/(60*dx);
DF(nx-2    ) = (-F(nx-5  ) + 9*F(nx-4  ) - 45*F(nx-3  ) + 45*F(nx-1  ) - 9*F(nx    ) + F(1   ))/(60*dx);
DF(nx-1    ) = (-F(nx-4  ) + 9*F(nx-3  ) - 45*F(nx-2  ) + 45*F(nx    ) - 9*F(1     ) + F(2   ))/(60*dx);
DF(nx      ) = (-F(nx-3  ) + 9*F(nx-2  ) - 45*F(nx-1  ) + 45*F(1     ) - 9*F(2     ) + F(3   ))/(60*dx);
DF(1       ) = (-F(nx-2  ) + 9*F(nx-1  ) - 45*F(nx    ) + 45*F(2     ) - 9*F(3     ) + F(4   ))/(60*dx);
DF(2       ) = (-F(nx-1  ) + 9*F(nx    ) - 45*F(1     ) + 45*F(3     ) - 9*F(4     ) + F(5   ))/(60*dx);
DF(3       ) = (-F(nx    ) + 9*F(1     ) - 45*F(2     ) + 45*F(4     ) - 9*F(5     ) + F(6   ))/(60*dx);
end

% 4th order center difference operator
function DF = CD4(F,dx,nx)
DF(3:nx-2,1)   = (F(1:nx-4) - 8*F(2:nx-3) + 8*F(4:nx-1) - F(5:nx))/(12*dx);
DF(nx-1)       = (F(nx-3)   - 8*F(nx-2)   + 8*F(nx)     - F(1)   )/(12*dx);
DF(nx)         = (F(nx-2)   - 8*F(nx-1)   + 8*F(1)      - F(2)   )/(12*dx);
DF(1)          = (F(nx-1)   - 8*F(nx  )   + 8*F(2)      - F(3)   )/(12*dx);
DF(2)          = (F(nx  )   - 8*F(1   )   + 8*F(3)      - F(4)   )/(12*dx);
end

% 3th order difference with 4th order diffusion
function DF = D3D4(F,dx,nx)
a = -157;
b = 1432;
c = -7195;
d = 0;
e = 7205;
f = -1448;
g = 163;
dh= 9600*dx;

DF(4:nx-3,1) = (a*F(1:nx-6) + b*F(2:nx-5) + c*F(3:nx-4) + d*F(4:nx-3) + e*F(5:nx-2) + f*F(6:nx-1) + g*F(7:nx))/dh;
DF(nx-2    ) = (a*F(nx-5  ) + b*F(nx-4  ) + c*F(nx-3  ) + d*F(nx-2  ) + e*F(nx-1  ) + f*F(nx    ) + g*F(1   ))/dh;
DF(nx-1    ) = (a*F(nx-4  ) + b*F(nx-3  ) + c*F(nx-2  ) + d*F(nx-1  ) + e*F(nx    ) + f*F(1     ) + g*F(2   ))/dh;
DF(nx      ) = (a*F(nx-3  ) + b*F(nx-2  ) + c*F(nx-1  ) + d*F(nx    ) + e*F(1     ) + f*F(2     ) + g*F(3   ))/dh;
DF(1       ) = (a*F(nx-2  ) + b*F(nx-1  ) + c*F(nx    ) + d*F(1     ) + e*F(2     ) + f*F(3     ) + g*F(4   ))/dh;
DF(2       ) = (a*F(nx-1  ) + b*F(nx    ) + c*F(1     ) + d*F(2     ) + e*F(3     ) + f*F(4     ) + g*F(5   ))/dh;
DF(3       ) = (a*F(nx    ) + b*F(1     ) + c*F(2     ) + d*F(3     ) + e*F(4     ) + f*F(5     ) + g*F(6   ))/dh;

end

% 2th order center difference operator
function DF = CD2(F,dx,nx)
DF(2:nx-1,1)   = (F(3:nx,1) - F(1:nx-2))/(2*dx);
DF(1)          = (F(2)      - F(nx    ))/(2*dx);
DF(nx)         = (F(1)      - F(nx-1  ))/(2*dx);
end

% 6th order diffusion
function DF = Diff6(F,nx)
DF(4:nx-3,1) = -(F(1:nx-6) - 6*F(2:nx-5) + 15*F(3:nx-4) - 20*F(4:nx-3,1) + 15*F(5:nx-2) - 6*F(6:nx-1) + F(7:nx))/(64);
DF(nx-2    ) = -(F(nx-5  ) - 6*F(nx-4  ) + 15*F(nx-3  ) - 20*F(nx-2    ) + 15*F(nx-1  ) - 6*F(nx    ) + F(1   ))/(64);
DF(nx-1    ) = -(F(nx-4  ) - 6*F(nx-3  ) + 15*F(nx-2  ) - 20*F(nx-1    ) + 15*F(nx    ) - 6*F(1     ) + F(2   ))/(64);
DF(nx      ) = -(F(nx-3  ) - 6*F(nx-2  ) + 15*F(nx-1  ) - 20*F(nx      ) + 15*F(1     ) - 6*F(2     ) + F(3   ))/(64);
DF(1       ) = -(F(nx-2  ) - 6*F(nx-1  ) + 15*F(nx    ) - 20*F(1       ) + 15*F(2     ) - 6*F(3     ) + F(4   ))/(64);
DF(2       ) = -(F(nx-1  ) - 6*F(nx    ) + 15*F(1     ) - 20*F(2       ) + 15*F(3     ) - 6*F(4     ) + F(5   ))/(64);
DF(3       ) = -(F(nx    ) - 6*F(1     ) + 15*F(2     ) - 20*F(3       ) + 15*F(4     ) - 6*F(5     ) + F(6   ))/(64);
end

% 4th order diffusion
function DF = Diff4(F,nx)
DF(4:nx-3,1) = (-F(1:nx-6) + 12*F(2:nx-5) - 39*F(3:nx-4) + 56*F(4:nx-3,1) - 39*F(5:nx-2) + 12*F(6:nx-1) - F(7:nx))/(16*6);
DF(nx-2    ) = (-F(nx-5  ) + 12*F(nx-4  ) - 39*F(nx-3  ) + 56*F(nx-2    ) - 39*F(nx-1  ) + 12*F(nx    ) - F(1   ))/(16*6);
DF(nx-1    ) = (-F(nx-4  ) + 12*F(nx-3  ) - 39*F(nx-2  ) + 56*F(nx-1    ) - 39*F(nx    ) + 12*F(1     ) - F(2   ))/(16*6);
DF(nx      ) = (-F(nx-3  ) + 12*F(nx-2  ) - 39*F(nx-1  ) + 56*F(nx      ) - 39*F(1     ) + 12*F(2     ) - F(3   ))/(16*6);
DF(1       ) = (-F(nx-2  ) + 12*F(nx-1  ) - 39*F(nx    ) + 56*F(1       ) - 39*F(2     ) + 12*F(3     ) - F(4   ))/(16*6);
DF(2       ) = (-F(nx-1  ) + 12*F(nx    ) - 39*F(1     ) + 56*F(2       ) - 39*F(3     ) + 12*F(4     ) - F(5   ))/(16*6);
DF(3       ) = (-F(nx    ) + 12*F(1     ) - 39*F(2     ) + 56*F(3       ) - 39*F(4     ) + 12*F(5     ) - F(6   ))/(16*6);
end

% 2nd order diffusion
function DF = Diff2(F,nx)
a = 2;
b = -27;
c = 270;
d = -490;
e = 270;
f = -27;
g = 2;
dh= 180*4;

DF(4:nx-3,1) = -(a*F(1:nx-6) + b*F(2:nx-5) + c*F(3:nx-4) + d*F(4:nx-3) + e*F(5:nx-2) + f*F(6:nx-1) + g*F(7:nx))/dh;
DF(nx-2    ) = -(a*F(nx-5  ) + b*F(nx-4  ) + c*F(nx-3  ) + d*F(nx-2  ) + e*F(nx-1  ) + f*F(nx    ) + g*F(1   ))/dh;
DF(nx-1    ) = -(a*F(nx-4  ) + b*F(nx-3  ) + c*F(nx-2  ) + d*F(nx-1  ) + e*F(nx    ) + f*F(1     ) + g*F(2   ))/dh;
DF(nx      ) = -(a*F(nx-3  ) + b*F(nx-2  ) + c*F(nx-1  ) + d*F(nx    ) + e*F(1     ) + f*F(2     ) + g*F(3   ))/dh;
DF(1       ) = -(a*F(nx-2  ) + b*F(nx-1  ) + c*F(nx    ) + d*F(1     ) + e*F(2     ) + f*F(3     ) + g*F(4   ))/dh;
DF(2       ) = -(a*F(nx-1  ) + b*F(nx    ) + c*F(1     ) + d*F(2     ) + e*F(3     ) + f*F(4     ) + g*F(5   ))/dh;
DF(3       ) = -(a*F(nx    ) + b*F(1     ) + c*F(2     ) + d*F(3     ) + e*F(4     ) + f*F(5     ) + g*F(6   ))/dh;
end