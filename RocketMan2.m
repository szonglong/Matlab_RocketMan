%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%static variables
ve1=2500; ve2=3550;%exhaust velocity
dm1=251; dm2=4.1; md2=4500; %md1=20000; %dry mass and mass flow rate
G=6.67*10^-11; RE=6371000; ME=5.972*10^24; %constants
A= 113; %area of rocket 


%initialise values at t=0
z0=0; %height above earth
v0=0; %intitial radial velocity
m0=407000; m1=15500; %initial fuel masses
vt0=460; st0=0; %tangential components
tau = 0.1; %interval
maxstep = 6000; % no. of iterations
alpha=0; %angle of thrust vectoring
at=0; %initial tangential velocity
gamma=2.7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assigning variables to init values
z=z0; v=v0; mf1=m0; a=0; mf2=m1; vt=vt0; st=st0;


for istep=1:maxstep
    t = (istep-1)*tau; 
    m=mf1+mf2+md2; %total mass
    rho=(1.225*exp(-(z)/8780)); %density function
    
    aFd=(0.5*0.82*rho*A*v^2)/m; %magnitude of frictional acceleration
    aFi1=ve1*(dm1/tau)/m; %magnitude of stage 1 thrust acceleration
    aFi2=ve2*(dm2/tau)/m; %magnitude of stage 2 thrust acceleration
    aFg=G*ME/((z+RE)^2); %magnitude of gravitational acceleration
    aFc=vt^2/(RE+z); %effective acceleration due to required centripetal force
    
    theta=atan2(vt-vt0,v); %angle of rocket
    phi=theta+alpha; %angle of thrust vector
    
    
    
    xplot(istep) = t;
    yplot(istep) = v;
   

    % while mf > 0
    if( mf1 > 0 )
        iter=istep;
        a_prime = aFi1 + aFc - aFg - aFd;
        aw= -(vt)*(1/(RE+z))*v -(vt)*dm1/m;
        at= aw+ aFi1;
        if (a_prime>0)
            a = a_prime;
        else
            a = 0;
        end
        mf1 = mf1-dm1;
    else
        if( mf2 > 0 )
            alpha=t*0.03;
            a = aFi2*cos(gamma) + aFc - aFg - aFd;
            mf2 = mf2-dm2;
            aw=-(vt)*(1/(RE+z))*v -(vt)*dm2/m;
            at= aw + aFi2*sin(gamma);
        else
            a= aFc -aFg - aFd;
        end
    end
    v = v+a*tau;
    z = z+v*tau;
    vt = vt + at*tau;
    st = st+vt*tau;
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

% 2.0 Make 2D: a_r, v_t
% 2.1 a_t
% 2.2 make single loop
% 2.3 consolidated variables
% 2.4 split into radial and tangential
% 2.5 realised that theta is problematic, try simple first/  Last one with
% (problematic) theta
% 2.6 Simplified version - thrust indep of angle. Single gamma on stage 2.
% Need to solve gamma st v=0 at end of burn

%to do:
% alpha(t) st vr = 0, vt = vreq, z = zreq; -> single variable v3
% v4 -> Powell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%