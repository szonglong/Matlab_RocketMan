%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%static variables
ve=100;
rho=1;
dm=1; md=60;
G=10; RE=50; ME=10;


%initialise values at t=0
z0=0;
v0=0;
m0=960;
tau = 0.01; %interval
maxstep = 2000; % no. of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assigning variables to init values
z=z0; v=v0; mf=m0; a=0


for istep=1:maxstep
    
    t = (istep-1)*tau; 
    xplot(istep) = t;
    yplot(istep) = z;
    m=mf+md

    % If ball reaches ground (y<0), break out of the loop
    if( mf > 0 )
        a = ve*dm/m - G*ME/((z+RE)^2);
        mf = mf-dm;
    else
        a= -G*ME/(RE^2)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
