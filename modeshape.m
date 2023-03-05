function [natf, avec, pf] = modeshape (L,modeno,intstp,E,I,M,thetadot,g)

c1 = 1;
c2 = 1;
c3 = 1;
detval(c1) = 0;
maxf = 40;


for omega = 0.1 : .001 : maxf

    Modmat = zeros(2,2);
    Modmat(1,1) = -cos(omega*L)-cosh(omega*L);
    Modmat(1,2) = -sin(omega*L)-sinh(omega*L);
    Modmat(2,1) = (-0^2)*(cos(omega*L)-cosh(omega*L))-(1)*(sin(omega*L)-sinh(omega*L));
    Modmat(2,2) = (-0^2)*(sin(omega*L)-sinh(omega*L))-(1)*(-cos(omega*L)-cosh(omega*L));

    natfe = 0;
    detval (c1+1) = det(Modmat);

    if detval(c1+1)*detval(c1) < 0
        natf(c2) = omega;
        c2 = c2 + 1;
        natfe = 1;
    end
    natfe;
    if (c2 == modeno + 1)
        break
    end

    c1 = c1 + 1;

end

for k = 1:c2-1
    omega = natf(k);

    m11 = -cos(omega*L)-cosh(omega*L);
    m12 = -sin(omega*L)-sinh(omega*L);


    p = -m11/m12;
    pf(k) = p;

    for xi = 0 :intstp : L
        intdec (c3) = ( cos(omega*xi)-cosh(omega*xi) + p*(sin(omega*xi) - sinh(omega*xi)))^2;
        c3 = c3 + 1;
    end

    intval = sum(intdec)*intstp;

    avec(k) = (1/intval)^(1/2);
end
natfe;
Modmat;