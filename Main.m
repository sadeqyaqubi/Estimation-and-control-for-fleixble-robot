clc; clear all; close all;
addpath(genpath('C:\Users\knsaya\Desktop\Toolboxes\yalmip\YALMIP-master'))
% PRELIMINARIES


rhob = 7800;
Eb = 2e11;
bb = 0.04;
hb = 0.06;
lb = 1;
Ab = bb*hb;
mb = rhob*Ab*lb;
Im = (1/3)*mb*lb^2;
Mp = 2;


modeno = 1;
intstp = 0.01;
g = 9.8;
Im = (1/3)*mb*lb^2;
Ib = (1/12)*bb*hb^3;


% CONTROLLER SETUP

T = 0.000075;
tf = 5;
stp = ceil(tf/T);

uthres = [200 , 200];

cfl = 1 ;
cfv = 0;

x(1,:) = zeros(1,2*(2+modeno));
x(2,:) = zeros(1,2*(2+modeno));
xs(1,:) = zeros(1,2*(2+modeno));
xs(2,:) = zeros(1,2*(2+modeno));
xs2(1,:) = zeros(1,2*(2+modeno));
xs2(2,:) = zeros(1,2*(2+modeno));
xs2(1:2,2+modeno) = 0.1;
xs3 = xs2;
wlref(1:4) = zeros(1,4) ;
u(1,:) = zeros(1,1);

ar = pi/5;
wr = 0.8*pi;

for k = 1:stp+1
    theta_r(k) = ar*cos(wr*k*T);
    thetadot_r(k) =  -ar*wr*sin(wr*k*T);
    thetaddot_r(k) =  -ar*wr^2*cos(wr*k*T);
end

for k = 1:3
    u(k,:) = 0;

end

for k = 1:3
    u(k,:) = 0;
    uob(k) = 0;
    fb(k,:) = 0;
end

% MODE SHAPE CALCULATION

xivec = 0:intstp: lb;

c1 = 1;
phi = zeros(modeno,lb/intstp + 1);
phiprime = zeros(modeno,lb/intstp + 1);
phisegond = zeros(modeno,lb/intstp + 1);
phikol = zeros(modeno,lb/intstp + 1);
phinel = zeros(modeno,lb/intstp + 1);
ximat = 0:intstp:lb;
lambda = .2;
alphac = 1;
omtildint = 0;
%% MAIN LOOP
for k=3:stp
    %%
    % SIMULATIONS
    t(k) = k*T;
    x_k = x(k-1,:);
    x_k1 = x(k-2,:);

    xs_k = xs(k-1,:);
    xs_k1 = xs(k-2,:);

    %
    %     xs2_k = xs2(k-1,:);
    %     xs2_k1 = xs2(k-2,:);

    u_k = [u(k-1,:) , zeros(1, modeno), fb(k-1)].';
    u_k1 = [u(k-2,:) , zeros(1, modeno), fb(k-2)].';


    phi0 = phi;
    phisegond0 = phisegond;
    phikol0 = phikol;
    phinel0 = phinel;

    [natf_k, avec_k, pf_k] = modeshape (lb,modeno,intstp,Eb,Ib,Mp,x_k(2+modeno),g);
    [phi_k, Imassv_k, Ixiq_k, Inel_k] = intcal (lb, avec_k, natf_k, pf_k);
    [natf_k1, avec_k1, pf_k1] = modeshape (lb,modeno,intstp,Eb,Ib,Mp,x_k1(2+modeno),g);
    [phi_k1, Imassv_k1, Ixiq_k1, Inel_k1] = intcal (lb, avec_k1, natf_k1, pf_k1);


    c1 = 1;
    for k1 = 1:modeno
        for k2 = 0:intstp: lb
            phi(k1,c1) = avec_k(k1)*(cos(natf_k(k1)*k2)-cosh(natf_k(k1)*k2)) + pf_k(k1)*avec_k(k1)*(sin(natf_k(k1)*k2)-sinh(natf_k(k1)*k2));

            phiprime(k1,c1) = avec_k(k1)*natf_k(k1)*(-sin(natf_k(k1)*k2)-sinh(natf_k(k1)*k2))...
                + pf_k(k1)*avec_k(k1)*natf_k(k1)*(cos(natf_k(k1)*k2)-cosh(natf_k(k1)*k2));

            phisegond(k1,c1) =  avec_k(k1)*(natf_k(k1)^2)*(-cos(natf_k(k1)*k2)-cosh(natf_k(k1)*k2))...
                + pf_k(k1)*avec_k(k1)*(natf_k(k1)^2)*(-sin(natf_k(k1)*k2)-sinh(natf_k(k1)*k2));

            phikol(k1,c1) =  avec_k(k1)*(natf_k(k1)^3)*(sin(natf_k(k1)*k2)-sinh(natf_k(k1)*k2))...
                + pf_k(k1)*avec_k(k1)*(natf_k(k1)^3)*(-cos(natf_k(k1)*k2)-cosh(natf_k(k1)*k2));

            phinel(k1,c1) =  avec_k(k1)*(natf_k(k1)^4)*(cos(natf_k(k1)*k2)-cosh(natf_k(k1)*k2))...
                + pf_k(k1)*avec_k(k1)*(natf_k(k1)^4)*(sin(natf_k(k1)*k2)-sinh(natf_k(k1)*k2));

            c1 = c1 + 1;
        end
        c1 = 1;
    end


    x(k,:) = syseq(x_k, x_k1, u_k, u_k1, phi_k, Imassv_k, Ixiq_k, Inel_k, phi_k1, Imassv_k1, Ixiq_k1, Inel_k1,....
        mb, lb, rhob, Ab, g, Mp, Eb, Ib, modeno, T, cfl);
    ks = 1;
    xs(k,:) = syseqs(xs_k, xs_k1, u_k, u_k1, phi_k, Imassv_k, Ixiq_k, Inel_k, phi_k1, Imassv_k1, Ixiq_k1, Inel_k1,....
        mb, lb, rhob, Ab, g, Mp, Eb, Ib, modeno, T, cfl, ks, x_k);
    kl = 10;
    kil = 5;
    %     xs2(k,:) = syseqs2(xs2_k, xs2_k1, u_k, u_k1, phi_k, Imassv_k, Ixiq_k, Inel_k, phi_k1, Imassv_k1, Ixiq_k1, Inel_k1,....
    %         mb, lb, rhob, Ab, g, Mp, Eb, Ib, modeno, T, cfl, kl, kil, x_k, uob(k-1), omtildint);

    xs3_k = x_k;
    xs3_k(2+modeno) = xs3(k-1,2+modeno);
    xs3_k(2*(2+modeno)) = xs3(k-1,2*(2+modeno));
    xs3(k,:) = x(k,:);


    %% obs

    %     ewl(k) = xs2(k,2+modeno) - x(k,2+modeno) ;

    x_k = xs(k,:);
    x_k1 = xs(k-1,:);

    theta_k = xs_k(1);
    eta_k = xs_k(2:1+modeno);
    wlr_k = xs_k(2+modeno);

    thetadot_k = xs_k(3+modeno);
    etadot_k = xs_k(4+modeno:2*(2+modeno)-1);
    wldotr_k = xs_k(2*(2+modeno));

    theta_k1 = xs_k1(1);
    eta_k1 = xs_k1(2:1+modeno);
    wlr_k1 = xs_k1(2+modeno);

    thetadot_k1 = xs_k1(3+modeno);
    etadot_k1 = xs_k1(4+modeno:2*(2+modeno)-1);
    wldotr_k1 = xs_k1(2*(2+modeno));

    x0_k = x(k,:);
    x0_k1 = x(k-1,:);

    theta0_k = x_k(1);
    eta0_k = x_k(2:1+modeno);
    wlr0_k = x_k(2+modeno);

    thetadot0_k = x_k(3+modeno);
    etadot0_k = x_k(4+modeno:2*(2+modeno)-1);
    wldotr0_k = x_k(2*(2+modeno));


    % CONTROL
    M_k = zeros(2+modeno,2+modeno);
    B_k = zeros(2+modeno,2+modeno);
    c1= 1;
    for k1 = 1:modeno
        for k2 = 1:modeno
            m_kp_k (c1) = Imassv_k(k1,k2)*eta_k(k1)*eta_k(k2);
            m_kdotp_k (c1) = Imassv_k(k1,k2)*eta_k(k1)*etadot_k(k2);
            c1 = c1 + 1;
        end
    end

    M_k = zeros(2+modeno,2+modeno);
    B_k = zeros(2+modeno,2+modeno);

    c1= 1;
    for k1 = 1:modeno
        for k2 = 1:modeno
            m_kp_k (c1) = Imassv_k(k1,k2)*eta_k(k1)*eta_k(k2);
            m_kdotp_k (c1) = Imassv_k(k1,k2)*eta_k(k1)*etadot_k(k2);
            c1 = c1 + 1;
        end
    end

    wl_k = x_k(2:1+modeno)*phi_k(:,end);

    wlgal(k) = wl_k;
    wdotl_k =  x_k(4+modeno:2*(2+modeno)-1)*phi_k(:,end);
    wxixixil_k = x_k(2:1+modeno)*phi_k(:,end);

    M_k(1,1) = Im+ Mp*lb^2 +  cfl*rhob*Ab*sum(m_kp_k);
    M_k(1,2:1+modeno) = cfl*rhob*Ab*Ixiq_k + cfl*Mp*lb*phi_k(:,end).';
    M_k(2:1+modeno,1) = cfl*Ixiq_k;
    M_k(2:1+modeno, 2:1+modeno) = Imassv_k;
    M_k(2+modeno,1) = lb;
    M_k(2+modeno,2+modeno) = 1;

    H_k(1) = 2*rhob*Ab*thetadot_k*sum(m_kdotp_k) + 2*Mp*(wl_k*wdotl_k) + .5*mb*g*lb*cos(theta_k) + Mp*g*lb*cos(theta_k) - Mp*g*sin(theta_k)*wl_k;
    H_k(2:1+modeno) = cfl*(Eb*Ib/rhob*Ab)*Inel_k*eta_k.' + cfl*thetadot_k^2*Imassv_k*eta_k.';
    H_k(2+modeno) = - wlr_k*thetadot_k^2 + g*cos(theta_k) - 0*(Eb*Ib/Mp)*wxixixil_k;

    B_k = zeros(2+modeno,2+modeno);
    B_k(1,1) = 1;
    B_k(end,end) = 1;



    [natfo_k, aveco_k, pfo_k] = modeshape (lb,modeno,intstp,Eb,Ib,Mp,x_k(2+modeno),g);
    [phio_k, Imassvo_k, Ixiqo_k, Inelo_k] = intcal (lb, aveco_k, natfo_k, pfo_k);


    %

    c1 = 1;
    for k1 = 1:modeno
        for k2 = 0:intstp: lb
            phio(k1,c1) = aveco_k(k1)*(cos(natfo_k(k1)*k2)-cosh(natfo_k(k1)*k2)) + pfo_k(k1)*aveco_k(k1)*(sin(natfo_k(k1)*k2)-sinh(natfo_k(k1)*k2));

            phiprimeo(k1,c1) = aveco_k(k1)*natfo_k(k1)*(-sin(natfo_k(k1)*k2)-sinh(natfo_k(k1)*k2))...
                + pf_k(k1)*aveco_k(k1)*natfo_k(k1)*(cos(natfo_k(k1)*k2)-cosh(natfo_k(k1)*k2));

            phisegondo(k1,c1) =  aveco_k(k1)*(natfo_k(k1)^2)*(-cos(natfo_k(k1)*k2)-cosh(natfo_k(k1)*k2))...
                + pfo_k(k1)*aveco_k(k1)*(natfo_k(k1)^2)*(-sin(natfo_k(k1)*k2)-sinh(natfo_k(k1)*k2));

            phikolo(k1,c1) =  aveco_k(k1)*(natfo_k(k1)^3)*(sin(natfo_k(k1)*k2)-sinh(natfo_k(k1)*k2))...
                + pfo_k(k1)*aveco_k(k1)*(natfo_k(k1)^3)*(-cos(natfo_k(k1)*k2)-cosh(natfo_k(k1)*k2));

            phinelo(k1,c1) =  aveco_k(k1)*(natfo_k(k1)^4)*(cos(natfo_k(k1)*k2)-cosh(natfo_k(k1)*k2))...
                + pfo_k(k1)*aveco_k(k1)*(natfo_k(k1)^4)*(sin(natfo_k(k1)*k2)-sinh(natfo_k(k1)*k2));

            c1 = c1 + 1;
        end
        c1 = 1;
    end

    w0o(k) = x(k,2:1+modeno)*phio(:,1);
    wdot0o(k) = x(k,4+modeno:2*(2+modeno)-1)*phio(:,1);
    wdotxi0o(k) = x(k,4+modeno:2*(2+modeno)-1)*phiprimeo(:,1);
    wxixi0o(k) = x(k,2:1+modeno)*phisegondo(:,1);
    wxixixi0o(k) = x(k,2:1+modeno)*phikolo(:,1);

    wlo(k) = x(k,2:1+modeno)*phio(:,end);
    wdotlo(k) = x(k,4+modeno:2*(2+modeno)-1)*phio(:,end);
    wdotxilo(k) = x(k,4+modeno:2*(2+modeno)-1)*phiprimeo(:,end);
    wxixilo(k) = x(k,2:1+modeno)*phisegondo(:,end);
    wxixixilo(k) = x(k,2:1+modeno)*phikolo(:,end);

    c1 = 1;
    c2 = 1;
    for k1 = 1:modeno
        for k2 = 1:modeno
            m_kdotp_k (c1) = Imassv_k(k1,k2)*eta_k(k1)*etadot_k(k2);
            c1 = c1 + 1;
        end
        m_ixiq_k(k1) = sum(phi(k1,:).*ximat)*intstp*eta_k(k1);
        c2 = c2 + 1;
    end

    e0_k = theta_k - theta_r(k) + lambda*(thetadot_k - thetadot_r(k));
    xi0 = e0_k*(thetadot_k - thetadot_r(k)-lambda*thetaddot_r(k));
    mu0 = lambda*e0_k;
    alphac = .9*e0_k;
    xi00 = -xi0 - alphac^2  ;

    thetaddotr = xi00/mu0;

    Mrr = M_k(1,1);
    Mrf = M_k(1,2:1+modeno);
    Mfr = M_k(2:1+modeno,1);
    Mff = M_k(2:1+modeno,2:1+modeno);

    Hr = H_k(1);
    Hf = H_k(2:1+modeno);

    Brr = 1;

    M0 = Mrr - Mrf*(Mff\Mfr);
    H0 = Hr.' - Mrf*(Mff\Hf.');

    delta0 = M0*thetaddotr+H0;

    ew_k = wlr_k + lambda*wldotr_k + lb*e0_k;
    xi1 = ew_k*( (wldotr_k) + lb*(thetadot_k - thetadot_r(k)-lambda*thetaddot_r(k)) );
    alphaw = .95*ew_k;
    xi11 = -xi1 - alphaw^2;
    mu1 = lambda*ew_k;
    zddotr = xi11/mu1;

    zddotr = (-((wldotr_k) + lb*(thetadot_k - thetadot_r(k)-lambda*thetaddot_r(k))) -.81*sign(ew_k))/lambda
    if (abs(ew_k) > 1e-3)
        fb(k) = zddotr + H_k(2+modeno) ;
    else
        fb(k) = 0;
    end
        fb(k) = zddotr + H_k(2+modeno) ;

%                     fb(k) = 0;

    u(k) = delta0 ;

    % observer 2
    xs3_k = xs3(k,:);
    xs3_k1 = xs3(k-1,:);

    theta_k = xs3_k(1);
    eta_k = xs3_k(2:1+modeno);
    wlr_k = xs3_k(2+modeno);

    thetadot_k = xs3_k(3+modeno);
    etadot_k = xs3_k(4+modeno:2*(2+modeno)-1);
    wldotr_k = xs3_k(2*(2+modeno));

    theta_k1 = xs3_k1(1);
    eta_k1 = xs3_k1(2:1+modeno);
    wlr_k1 = xs3_k1(2+modeno);

    thetadot_k1 = xs3_k1(3+modeno);
    etadot_k1 = xs3_k1(4+modeno:2*(2+modeno)-1);
    wldotr_k1 = xs3_k1(2*(2+modeno));

    wlref(k) = x_k(2+modeno);
    wlref(k+1) = x_k(2+modeno) + T*x_k(2*(2+modeno));
    wlref(k+2) =  wlref(k+1)+T*x_k(end);

    f0_k = - wlr_k*thetadot_k^2 + g*cos(theta_k) ;



    lambda0 = 0.1;
    s(k+1) = (wlr_k+T*wldotr_k) + lambda0*wlr_k - wlref(k+1) - lambda0*wlref(k);



    xi = wlr_k + 2*T*wldotr_k + (T^2)*f0_k(end) + lambda0*wlr_k + lambda0*T*wldotr_k - wlref(k+2) - lambda0*wlref(k+1);
    bo = T^2;
    alphas = 0.99;
    uob(k) = (alphas*abs(s(k+1)) - xi) / bo;


    %%
    Progress = (t(k)/tf)*100

end
%% PLOTS
close all;
figure; plot(t,x(:,1)); hold all; plot(t,theta_r(1:end-1)); legend('System Response', 'Reference Signal'); xlabel('Time [s]'); ylabel ('\theta [rad]'); axis([0, stp*T, 1.2*min(x(:,1)), 1.2*max(x(:,1))])
figure; plot(t,x(:,3+modeno)); hold all; plot(t,thetadot_r(1:end-1)); legend('System Response', 'Reference Signal'); xlabel('Time [s]'); ylabel ('\theta\dot [rad]');  axis([0, stp*T, 1.2*min(x(:,3+modeno)), 1.2*max(x(:,3+modeno))])


figure; plot(t,u(1:end,1)); legend('Input torque'); xlabel('Time [s]'); ylabel ('Input torque [N.m]'); axis([0, stp*T, 60, 170])
figure; plot(t,fb); xlabel('Time [s]'); ylabel ('Boundary input [N]');axis([0, stp*T, 4, 12])

figure; plot(t,x(:,2+modeno)); hold all; plot(t,xs3(:,2+modeno));  xlabel('Time [s]'); ylabel ('\omega_l [rad]'); legend('System Response', 'Estimate'); axis([0, stp*T, 1.2*min(xs3(:,2+modeno)), 1.2*max(xs3(:,2+modeno))])
figure; plot(t,x(:,2+modeno)); xlabel('Time [s]'); ylabel ('\omega_l [rad]'); legend('System Response', 'Estimate'); axis([0, stp*T, 0.1, -.6])
 figure; plot(t,wlgal);xlabel('Time [s]'); ylabel ('galerkin [rad]'); axis([0, stp*T, -.6, -.6])

%%
etavec = xs3(:,2+modeno) ./ phi(:,end);

cd = 1;
for k1 = 1:max(stp)
    for k2 = 1:max(size(phi))
    progress2 = (k1/stp)*100 
    disp2(k1,k2)=etavec(k1)*phi(k2);
end
end
    disp2(k1,k2)=etavec*phi;

[X,T0] = meshgrid(xivec,t);
figure; surf(X,T0,disp2);xlabel('\xi [m]'); ylabel ('Time [s]'); zlabel ('\omega [m]'); axis([0, 1, 0, 5, -.1, .12])

%%
etavec = xs3(:,2+modeno) ./ phi(:,end);

disp2=etavec*phi;

etavec2 = x(:,2+modeno) ./ phi(:,end);
disp3 = etavec2*phi;

[X,T0] = meshgrid(xivec,t);
figure; surf(X,T0,abs(disp3-disp2)); xlabel('\xi [m]'); ylabel ('Time [s]'); zlabel ('Estimation error magnitude [m]'); axis([0, 1, 0, 5, 0, .12])
figure; surf(X,T0,abs(disp3)); xlabel('\xi [m]'); ylabel ('Time [s]'); zlabel ('Link displacement [m]'); axis([0, 1, 0, 5, 0, .12])
