%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%static variables
ve=4400;
rho=1;
dm=251; md=20000;
G=6.67*10^-11; RE=6371000; ME=5.972*10^24;


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
    m=mf+md
    xplot(istep) = t;
    yplot(istep) = z;
   

    % If ball reaches ground (y<0), break out of the loop
    if( mf > 0 )
        a_prime = ve*(dm/tau)/m - G*ME/((z+RE)^2);
        if (a_prime>0)
            a = a_prime;
        else
            a = 0;
        end
        mf = mf-dm;
    else
        a= -G*ME/((z+RE)^2);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
