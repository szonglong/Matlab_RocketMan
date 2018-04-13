%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RocketMan! 
% Describes the trajectory of a rocket shot vertically upward
% V2 is a one parameter solver for 
% 
% NOTE: clear variables before each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long


% Init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
init_theta_range = 1.5:0.1:2.6;
init_ve2_range = 4400:200:5400;

init_theta_fat=max(init_theta_range)-min(init_theta_range);
init_ve2_fat=max(init_ve2_range)-min(init_ve2_range);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tiny_count = 0;

theta_range=init_theta_range;%2.02574;
ve2_range=init_ve2_range;%4622;
theta_fat=init_theta_fat;
ve2_fat=init_ve2_fat;


while tiny_count < 10
    theta_fat=max(theta_range)-min(theta_range);
    ve2_fat=max(ve2_range)-min(ve2_range);


    %generating arrays
    xxzrange=zeros();
    xxtrange=zeros();
    yyzrange=zeros();
    yytrange=zeros();

    xxzindex = 1; yyzindex = 1; xxtindex = 1; yytindex = 1;

    for mode=["vz" "vt"]
    array=zeros();
    index2 = 1;    
        for i=ve2_range
            index = 1;
            for j = theta_range
                array(index,index2)=fun(j,i,mode);
                if index>1
                    if array(index,index2)*array(index-1,index2) < 0
                        if mode == "vz"
                            xxzrange(xxzindex,:) = i;
                            xxzindex = xxzindex + 1;
                            yyzrange(yyzindex,:) = j;
                            yyzindex = yyzindex + 1;
                            break
                        end
                        if mode == "vt"
                            xxtrange(xxtindex,:) = i;
                            xxtindex = xxtindex + 1;
                            yytrange(yytindex,:) = j;
                            yytindex = yytindex + 1;

                            break
                        end
                    end
                end
                index = index +1;

            end
            index2 = index2 + 1;
        end

    end


    global lt lz
    lt = lagranp(xxtrange,yytrange);
    lz = lagranp(xxzrange,yyzrange);

    new_ve2 = Fit_and_solve()
    if 1.5 < polyval(lt, new_ve2) && polyval(lt, new_ve2) < 3.14
        new_theta = polyval(lt, new_ve2)
        theta_range = (new_theta-0.1*theta_fat):(theta_fat/20):(new_theta+0.1*theta_fat);
    else
        break
    end
    ve2_range = (new_ve2-0.1*ve2_fat):(ve2_fat/20):(new_ve2+0.1*ve2_fat);
    
    tiny_count = tiny_count + 1;
end



function ans=Fit_and_solve()
    format long
    TOL=1e-8;
    Error=1e10;
    Left=4619; Right=4626; % Initialization of variables

    fLeft = f2(Left); fRight = f2(Right);
    while (Error>TOL)
        Middle = (Left+Right)/2; fMiddle = f2(Middle);
        if ( fLeft * fMiddle <= 0 )
            Right=Middle; fRight=fMiddle;
        else
            Left=Middle; fLeft=fMiddle;
        end
        Error=abs( (Right-Left)/Middle );
    end
    ans=Middle;
end

function f=f2(x)
    global lt lz
    f=polyval(dif_poly_coeff(lt,lz),x);
end

function coeff_diff = dif_poly_coeff(x1, x2)
    x1_order = length(x1);
    x2_order = length(x2);
    if x1_order > x2_order
         max_order = size(x1);
    else
         max_order = size(x2);
    end
    new_x1 = padarray(x1,max_order-size(x1),0,'pre');
    new_x2 = padarray(x2,max_order-size(x2),0,'pre');
    coeff_diff = new_x1 - new_x2;
    return
end

function f=fun(theta,ve2,mode)
    %static variables
    ve1=2500; %ve2=4700;%3550;%exhaust velocity
    dm1=251; dm2=4.1; md2=4000; %md1=20000; %dry mass and mass flow rate
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
    maxstep = 8000; % no. of iterations
    alpha=0; %angle of thrust vectoring
    at=0; %initial tangential velocity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %assigning variables to init values
    z=z0; vz=vz0; mf1=m0; az=0; mf2=m1; st=st0; vt=vt0; vrt=vrt0;

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
% 
%         xplot(istep) = t;
%         yplot(istep) = (vt+vw0)-vreq;
%         yplot2(istep)= vz;%(vt+vw0)-vreq;


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
%                 break %when fuel runs out, cut the loop
                az=  aFc -aFg - aFd;
                art=- aFd*sin(theta);
            end
        end

        vz = vz+az*tau;
        z = z+vz*tau;

        vw = 460- ((vw0*RE)/((RE+z)));
        vrt= vrt + art*tau;
        vt = vrt + vw;
        st = st+vt*tau;
    end


    % plot1=plot(xplot,yplot,'r-',xplot,yplot2,'b-');
    % xlabel('t'); ylabel('z');
    if mode == "vz"
        f=vz;
    end
    if mode == "vt"
        f=(vt+vw0)-vreq;
    end
end

function l = lagranp(x,y)
    %Input : x = [x1 ... xN], y = [y1 ... yN]
    %Output: l = Lagrange polynomial coefficients of degree N-1
    l = 0;
    for m = 1:length(x)
        P = 1;
        for k = 1:length(x)
            if k ~= m, P = conv(P,[1 -x(k)])/(x(m)-x(k)); end
        end
        L(m,:) = P; %cardinal function
        l = l + y(m)*P; %Lagrange polynomial
    end
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
% 2.11 returns matrix of vz for range of theta and ve2

% 3.0.0 band comparison method (brute force): forming grandarrays
% 3.0.1 point generator
% 3.0.2 interpolation and root solving for one example
% 3.0.3 arranging to neater
% 3.0.4 2 variable solver to get required vr and vt
% 3.0.5 Optimised solver fpr 2 variable

% v4 -> Powell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%