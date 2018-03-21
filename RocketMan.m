%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%static variables
ve=10;
rho=1;
dm=1;

%initialise values at t=0
z0=0;
v0=0;
m0=540;
tau = 0.01; %interval
maxstep = 1000; % no. of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assigning variables to init values
z=z0; v=v0; m=m0;


for istep=1:maxstep
    
    t = (istep-1)*tau; 
    xplot(istep) = t;
    yplot(istep) = z;
%     yplot2(istep) = y;
    

    a = ve*dm/m;
    m = m-dm;
    v = v+a*tau;
    z = z+v*tau
    
    
    
    % If ball reaches ground (y<0), break out of the loop
    if( m <= 0 )
        xplot(istep+1) = t; % Record last values computed
        yplot(istep+1) = z;
        break; % Break out of the for loop
    end

end

plot1=plot(xplot,yplot,'-');
xlabel('t'); ylabel('z');












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changelog:
% 1.0 Initialise program. Basic acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
