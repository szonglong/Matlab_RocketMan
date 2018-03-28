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
A= 113; atc=9;


%initialise values at t=0
z0=0;
v0=0;
m0=407000; m1=15500;
vt0=460; st0=0;
tau = 0.1; %interval
maxstep = 6000; % no. of iterations
alpha=0; 
at=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assigning variables to init values
z=z0; v=v0; mf1=m0; a=0; mf2=m1; vt=vt0; st=st0;


for istep=1:maxstep
    t = (istep-1)*tau; 
    m=mf1+mf2+md2;
    rho=(1.225*exp(-(z)/8780));
    
    aFd=(0.5*0.82*rho*A*v^2)/m;
    aFi1=ve1*(dm1/tau)/m;
    aFi2=ve2*(dm2/tau)/m;
    aFg=G*ME/((z+RE)^2);
    aFc=vt^2/(RE+z);
    
    theta=atan2(vt-vt0,v);
    phi=theta+alpha;
    
    aw=(vt0)*(RE/(RE+z)^2)*v;
    
    xplot(istep) = t;
    yplot(istep) = at;
   

    % while mf > 0
    if( mf1 > 0 )
        iter=istep
        a_prime = aFi1*cos(phi) + aFc - aFg - aFd*cos(theta);
        at= -aw % +aFi1*sin(phi)% - aFd*sin(theta);
        if (a_prime>0)
            a = a_prime;
        else
            a = 0;
        end
        mf1 = mf1-dm1;
    else
        if( mf2 > 0 )
            alpha=t*0.0
            a = aFi2*cos(phi) + aFc - aFg - aFd*cos(theta);
            mf2 = mf2-dm2;
            at= -aw % +aFi2*sin(phi) %- aFd*sin(theta) ;
        else
            a= aFc -aFg - aFd*cos(theta);
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

%to do:
% alpha(t) st vr = 0, vt = vreq, z = zreq; -> single variable v3
% v4 -> Powell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%