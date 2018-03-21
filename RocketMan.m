%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%static variables
ve=2500;
dm=251; md=20000;
G=6.67*10^-11; RE=6371000; ME=5.972*10^24;
A= 113;


%initialise values at t=0
z0=0;
v0=0;
m0=407000;
tau = 0.1; %interval
maxstep = 2000; % no. of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assigning variables to init values
z=z0; v=v0; mf=m0; a=0;


for istep=1:maxstep
    
    t = (istep-1)*tau; 
    m=mf+md;
    rho=(1.225*exp(-(z)/8780));
    xplot(istep) = t;
%     yplot(istep) = z;
    yplot(istep) = a/9.81;
    yplot2(istep) = 0.5*rho*v^2;
   

    % while mf > 0
    if( mf > 0 )
        a_prime = ve*(dm/tau)/m - G*ME/((z+RE)^2) - (0.5*0.82*rho*A*v^2)/m;
        if (a_prime>0)
            a = a_prime;
        else
            a = 0;
        end
        mf = mf-dm;
    else
        a= -G*ME/((z+RE)^2) - (0.5*0.82*rho*A*v^2)/m;
    end
    v = v+a*tau;
    z = z+v*tau;
end

plot1=plot(xplot,yplot,'-');
%plot2=plot(xplot,yplot2,'x');
xlabel('t'); ylabel('z');












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changelog:
% 1.0 Initialise program. Basic acceleration
% 1.1 Added in dry mass and gravitational potential
% 1.2 Added in real world values (Falcon Heavy)
% 1.3 Added in drag force
% 1.4 Verified with Falcon Heavy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
