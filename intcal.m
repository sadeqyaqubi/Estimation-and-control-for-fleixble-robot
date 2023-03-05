function [phi, Imassv, Ixiq, Inel] = intcal (L, avec, natf, pf)

intstp = 0.01;

c1 = 1;

modeno = max(size(avec));

ximat = 0:intstp:L;


for k1 = 1:modeno
    for k2 = 0:intstp: L 
        phi(k1,c1) = avec(k1)*(cos(natf(k1)*k2)-cosh(natf(k1)*k2)) + pf(k1)*avec(k1)*(sin(natf(k1)*k2)-sinh(natf(k1)*k2));
        phinel(k1,c1) =  avec(k1)*(natf(k1)^4)*(cos(natf(k1)*k2)-cosh(natf(k1)*k2))...
                + pf(k1)*avec(k1)*(natf(k1)^4)*(sin(natf(k1)*k2)-sinh(natf(k1)*k2));

        c1 = c1 + 1;
    end
    c1 = 1;
end
     
for k1 = 1:modeno
%     Iq(k1) = sum(phi(k1,:))*intstp;
%     I2q(k1) = sum(phi(k1,:).^2)*intstp;

    for k2 = 1:modeno
   
        Imassv(k1,k2) = sum(phi(k1,:).*phi(k2,:))*intstp;
%         Imv = Imassv;
%         if k1 == k2
%             Imv(k1,k1) = 0;
%         end
        Inel(k1,k2) = sum(phinel(k1,:).*phi(k2,:))*intstp;
    end

%     Imv(1,1) = 0;

    Ixiq(k1) = sum(phi(k1,:).*ximat)*intstp;

end
