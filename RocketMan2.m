%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array=zeros;
index = 1;
for i=1.6:0.1:2.1
    array(index,:)=fun(i);
    index = index +1;
    i
end
array

function f=fun(theta)
%static variables
ve1=2500; ve2=4600;%3550;%exhaust velocity
dm1=251; dm2=4.1; md2=4500; %md1=20000; %dry mass and mass flow rate
G=6.67*10^-11; RE=6371000; ME=5.972*10^24; %constants
A= 113; %area of rocket 


%initialise values at t=0
z0=0; %height above earth
vz0=0; %intitial radial velocity
vt0=0; %initial tangential velocity
vrt0=0; %initial rocket tangential velocity
m0=407000; m1=15500; %initial fuel masses
vw0=460; st0=0; %tangential components
tau = 0.1; %interval
maxstep = 24000; % no. of iterations
alpha=0; %angle of thrust vectoring
at=0; %initial tangential velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assigning variables to init values
z=z0; vz=vz0; mf1=m0; az=0; mf2=m1; st=st0; vt=vt0; vrt=vrt0;

%v2 "hand of God" theta
%theta=2.039239112185;

for istep=1:maxstep
    t = (istep-1)*tau; 
    m=mf1+mf2+md2; %total mass
    rho=(1.225*exp(-(z)/8780)); %density function
    vreq = sqrt(G*ME/(RE+z));
    
    aFd=(0.5*0.82*rho*A*(vz^2+(vt-vw0*(RE/(RE+z))+vw0)^2))/m; %magnitude of frictional acceleration
    aFi1=ve1*(dm1/tau)/m; %magnitude of stage 1 thrust acceleration
    aFi2=ve2*(dm2/tau)/m; %magnitude of stage 2 thrust acceleration
    aFg=G*ME/((z+RE)^2); %magnitude of gravitational acceleration
    aFc=(vt+vw0)^2/(RE+z); %effective acceleration due to required centripetal force
    
    xplot(istep) = t;
    yplot(istep) = z;
    yplot2(istep)= vz;%(vt+vw0)-vreq;
   

    % stage 1
    if( mf1 > 0 ) %while there is fuel in stage 1
        a_prime = aFi1 + aFc - aFg - aFd;
        art= 0;
        if (a_prime>0) %if upward force > downward force
            az = a_prime;
        else
            az = 0;
        end
        mf1 = mf1-dm1;
    else
        if( mf2 > 0 )
            az = aFi2*cos(theta) + aFc - aFg - aFd;
            mf2 = mf2-dm2;
            art = aFi2*sin(theta) - aFd*sin(theta);
        else
            break %when fuel runs out, cut the loop
%             az=  aFc -aFg - aFd;
%             art=- aFd*sin(theta);
        end
    end
    
    vz = vz+az*tau;
    z = z+vz*tau;
    
    vw = 460- ((vw0*RE)/((RE+z)));
    vrt= vrt + art*tau;
    vt = vrt + vw;
    st = st+vt*tau;
end


% plot1=plot(xplot,yplot,'-',xplot,yplot2,'-');
% xlabel('t'); ylabel('z');
f=vz;
end

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
% 2.7 vt is earth moving under rocket.. all equations and variables
% redefined
% 2.8 remove aw, explicit vw
% 2.9 revised some expressions, should be a basic correct one
% 2.10 returns array of vz for range of theta

%to do:
% alpha(t) st vr = 0, vt = vreq, z = zreq; -> single variable v3
% v4 -> Powell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%