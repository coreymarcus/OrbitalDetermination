%%% Part ii

mu=4902.801076;
r=1738;
mur=mu/r;
R=1737.1;
load('moon2.mat')

% Unnormalize Coefficients
Cn=C;Sn=S;Jn=J;
for l=1:20
    Nl0=sqrt(1/(2*l+1));
    Jn(l)=J(l)/Nl0;
    for m=1:l
        Nlm = sqrt(factorial(l+m)/((2)*factorial(l-m)*(2*l+1)));
        Sn(l,m)=S(l,m)/Nlm;
        Cn(l,m)=C(l,m)/Nlm;
    end
end

deg=10;
ph=0.1;
la=0.1;
dUdr=0; dUdph=0; dUdla=0;
for l=1:deg
    Plm=legendre(l,sin(ph));
    for m=0:l
        if mod(m+1,2)==0
            Plm(m+1)=-Plm(m+1); % Remove Condon-Shortley Phase
        end
        if m==0
            if l==1
                cPlmp=0;
            else
                cPlmp = -l*sin(ph)*sec(ph)*Plm(m+1)+(l+m)*sec(ph)*Plmprev(m+1);
            end
            dUdr = dUdr + (R/r)^l*(l+1)*Plm(m+1)*(-Jn(l)*cos(m*la));
            dUdla = dUdla + (R/r)^l*m*Plm(m+1)*(Jn(l)*sin(m*la));
            dUdph = dUdph + (R/r)^l*cPlmp*(-Jn(l)*cos(m*la));
            %if m+1<=l
            %    dUdph = dUdph + (R/r)^l*(l+1)*(Plm(m+2)-m*tan(ph)*Plm(m+1))*(-Jn(l)*cos(m*la));
            %else
            %    %dUdph = dUdph + (R/r)^l*(l+1)*(0-m*tan(ph)*Plm(m+1))*(-Jn(l)*cos(m*la)); %?
            %end
        else
            if l~=1 && m~=l
                cPlmp = -l*sin(ph)*sec(ph)*Plm(m+1)+(l+m)*sec(ph)*Plmprev(m+1);
            else
                cPlmp = -l*sin(ph)*sec(ph)*Plm(m+1);
                %cPlmp=0;
            end
            dUdr = dUdr + (R/r)^l*(l+1)*Plm(m+1)*(Cn(l,m)*cos(m*la)+Sn(l,m)*sin(m*la));
            dUdla = dUdla + (R/r)^l*m*Plm(m+1)*(Sn(l,m)*cos(m*la)-Cn(l,m)*sin(m*la));
            dUdph = dUdph + (R/r)^l*cPlmp*(Cn(l,m)*cos(m*la)+Sn(l,m)*sin(m*la));
            %if m+1<=l
            %    dUdph = dUdph + (R/r)^l*(l+1)*(Plm(m+2)-m*tan(ph)*Plm(m+1))*(Cn(l,m)*cos(m*la)+Sn(l,m)*sin(m*la));
            %else
            %    %dUdph = dUdph + (R/r)^l*(l+1)*(0-m*tan(ph)*Plm(m+1))*(Cn(l,m)*cos(m*la)+Sn(l,m)*sin(m*la)); %?
            %end
        end
    end
    Plmprev=Plm;
end
dUdr = -mu/r^2*dUdr;
dUdla = mu/r*dUdla;
dUdph = mu/r*dUdph;

dUdla=1/(r*cos(ph))*dUdla;
dUdph=1/r*dUdph;

gradUspher=[dUdr;dUdla;dUdph];


