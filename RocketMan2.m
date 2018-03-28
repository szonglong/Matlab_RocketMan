%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%static variables
ve1=2500; ve2=3550;
dm1=251; dm2=4.1; md1=20000; md2=4500;
G=6.67*10^-11; RE=6371000; ME=5.972*10^24;
A= 113;


%initialise values at t=0
z0=0;
v0=0;
m0=407000; m1=15500;
tau = 0.1; %interval
maxstep = 6000; % no. of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assigning variables to init values
z=z0; v=v0; mf1=m0; a=0; mf2=m1;


for istep=1:maxstep
    iter=istep;
    t = (istep-1)*tau; 
    m=mf1+md1;
    rho=(1.225*exp(-(z)/8780));
    xplot(istep) = t;
    yplot(istep) = z;
%     yplot(istep) = a/9.81;
%     yplot(istep) = 0.5*rho*v^2;
   

    % while mf > 0
    if( mf1 > 0 )
        a_prime = ve1*(dm1/tau)/m - G*ME/((z+RE)^2) - (0.5*0.82*rho*A*v^2)/m;
        if (a_prime>0)
            a = a_prime;
        else
            a = 0;
        end
        mf1 = mf1-dm1
    else
        break
    end
    v = v+a*tau;
    z = z+v*tau;
end

for istep=iter:maxstep
    t = (istep-1)*tau; 
    m=mf2+md2;
    rho=(1.225*exp(-(z)/8780));
    xplot(istep) = t;
    yplot(istep) = z;
%     yplot(istep) = a/9.81;
%     yplot(istep) = 0.5*rho*v^2;
   

    if( mf2 > 0 )
        a = ve2*(dm2/tau)/m - G*ME/((z+RE)^2) - (0.5*0.82*rho*A*v^2)/m;
        mf2 = mf2-dm2;
        f=istep
    else
        a= -G*ME/((z+RE)^2) - (0.5*0.82*rho*A*v^2)/m;
    end
    v = v+a*tau;
    z = z+v*tau;
    
end


plot1=plot(xplot,yplot,'-');
xlabel('t'); ylabel('z');












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changelog:
% 1.0 Initialise program. Basic acceleration
% 1.1 Added in dry mass and gravitational potential
% 1.2 Added in real world values (Falcon Heavy)
% 1.3 Added in drag force
% 1.4 Verified with Falcon Heavy
% 1.5 Added in stage 2

% 2.0 Make 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%