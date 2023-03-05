function y = syseqs2 (x_k, x_k1, u_k, u_k1, phi_k, Imassv_k, Ixiq_k, Inel_k, phi_k1, Imassv_k1, Ixiq_k1, Inel_k1,....
        mb, lb, rhob, Ab, g, Mp, Eb, Ib, modeno, T, cfl, kl, kil, x0_k, uob_k, omtildint)
% u_k(end) = 0
theta_k = x_k(1);
eta_k = x_k(2:1+modeno);
% eta0_k = x0_k(2:1+modeno);

wlr_k = x_k(2+modeno);

thetadot_k = x_k(3+modeno);
etadot_k = x_k(4+modeno:2*(2+modeno)-1);
wldotr_k = x_k(2*(2+modeno));

theta_k1 = x_k1(1);
eta_k1 = x_k1(2:1+modeno);
wlr_k1 = x_k1(2+modeno);

thetadot_k1 = x_k1(3+modeno);
etadot_k1 = x_k1(4+modeno:2*(2+modeno)-1);
wldotr_k1 = x_k1(2*(2+modeno));

wlr0_k = x0_k(2+modeno);

Im = (1/3)*mb*lb^2;

kl;
%%
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
wdotl_k =  x_k(4+modeno:2*(2+modeno)-1)*phi_k(:,end);
wxixixil_k = x_k(2:1+modeno)*phi_k(:,end);

M_k(1,1) = Im+ Mp*lb^2 +  cfl*rhob*Ab*sum(m_kp_k);
%  Im
%  Mp
%  lb
%   cfl*rhob*Ab*sum(m_kp_k)
M_k(1,2:1+modeno) = cfl*rhob*Ab*Ixiq_k + cfl*Mp*lb*phi_k(:,end).';
M_k(2:1+modeno,1) = cfl*Ixiq_k;
M_k(2:1+modeno, 2:1+modeno) = Imassv_k;
M_k(2+modeno,1) = lb;
M_k(2+modeno,2+modeno) = 1;


% uob_k
H_k(1) = 2*rhob*Ab*thetadot_k*sum(m_kdotp_k) + 2*Mp*(wl_k*wdotl_k) + .5*mb*g*lb*cos(theta_k) + Mp*g*lb*cos(theta_k) - Mp*g*sin(theta_k)*wl_k;
H_k(2:1+modeno) = cfl*(Eb*Ib/rhob*Ab)*Inel_k*eta_k.' + cfl*thetadot_k^2*Imassv_k*eta_k.';% + ks*(eta_k-eta0_k).';
H_k(2+modeno) = - wlr_k*thetadot_k^2 + g*cos(theta_k) - 0*(Eb*Ib/Mp)*wxixixil_k + kl*(wlr_k-wlr0_k) + kil*omtildint; 

B_k = zeros(2+modeno,2+modeno);
B_k(1,1) = 1;
B_k(1,end) = 0*lb-0*cos(theta_k);
B_k(end,end) = 1;


%%
% M_k1 = zeros(2+modeno,2+modeno);
% B_k1 = zeros(2+modeno,2+modeno);
% 
% c1= 1;
% for k1 = 1:modeno
%     for k2 = 1:modeno
%         m_kp_k1 (c1) = Imassv_k1(k1,k2)*eta_k1(k1)*eta_k1(k2);
%         m_kdotp_k1 (c1) = Imassv_k1(k1,k2)*eta_k1(k1)*etadot_k1(k2);
%         c1 = c1 + 1;
%     end
% end
% 
% wl_k1 = x_k1(2:1+modeno)*phi_k1(:,end);
% wdotl_k1 = x_k1(5+modeno:end)*phi_k1(:,end);
% wxixixil_k1 = x_k1(2:1+modeno)*phi_k1(:,end);
% 
% M_k1(1,1) = Im+ Mp*lb^2 +  cfl*rhob*Ab*sum(m_kp_k1);
% M_k1(1,2:1+modeno) = cfl*rhob*Ab*Ixiq_k1 + cfl*Mp*lb*phi_k1(:,end);
% M_k1(2:1+modeno,1) = cfl*Ixiq_k1;
% M_k1(2:1+modeno, 2:1+modeno) = Imassv_k1;
% M_k1(2+modeno,1) = lb;
% M_k1(2+modeno,2+modeno) = 1;
% 
% H_k1(1) = 2*rhob*Ab*thetadot_k1*sum(m_kdotp_k1) + 2*Mp*(wl_k1*wdotl_k1) + .5*mb*g*lb*cos(theta_k1) + Mp*g*lb*cos(theta_k1) - Mp*g*sin(theta_k1)*wl_k1;
% H_k1(2:1+modeno) = cfl*(Eb*Ib/rhob*Ab)*Inel_k1*eta_k1.' + cfl*thetadot_k1^2*Imassv_k1*eta_k1.';
% H_k1(2+modeno) = - wlr_k1*thetadot_k1^2 + g*cos(theta_k1) - (Eb*Ib/Mp)*wxixixil_k1;
% 
% B_k1(1,1) = 1;
% B_k1(1,end) = lb-wl_k1*cos(theta_k1);
% B_k1(end,end) = 1;

%%
f0_k = x_k(3+modeno:2*(2+modeno)).';
% M_k
% H_k
% B_k
% u_k



f1_k = M_k\(B_k*u_k-H_k');
% f0_k1 = x_k1(2+modeno:2*(1+modeno)).';
% f1_k1 = M_k1\(B_k*u_k1-H_k1');
f_k = [f0_k ; f1_k];
% f_k1 = [f0_k1 ; f1_k1];
% y = x_k.' + (3*T/2)*f_k - (T/2)*f_k1;
%  f0_k
%  f1_k

 y = x_k.' + T*f_k ;

 y;


