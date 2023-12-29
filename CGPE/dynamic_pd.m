re = reycol;
im = imycol;
alpha = 4.4;
sig = 0.3;
mu = 3*alpha / (2*sig);
%real scale
Rb = 15.0;
R0 = 5.0;

%%%%%%%%set up parameters%%%%%%%
%%%%time step%%%%%%%%
dt=0.01;
t_max=100.0;
%%%%%computing domain%%%%%
a=-15;
b=15;
%%%%%mesh size&&grid pints%%%%
N=400;
h=(b-a)/(N-1);
index=1:(N);
x=a:h:b;
y=a:h:b;
r=zeros(N,N);
%%%%%equation parameters%%%%%%
mul=(1/N)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%initialize the starting wave function%%%
psi0=zeros(N,N);
v2d=zeros(N,N);
HDP=zeros(N,N);
xi=zeros(1,N);
km2=zeros(N,N);
%the inverse space; after FFT
for i=1:(N/2)
    %xi include xi,yj on grid
    xi(i)=(i-1)*2*pi/(b-a);%half of the momentum space
    xi(i+N/2)=((i-N/2-1)*2*pi/(b-a));%the other half folded into -km
end
%axis for meshgrid plot
[xx,yy]=meshgrid(x, x);
[xxi,yyi]=meshgrid(xi, xi);
for i=1:N
    for j=1:N
        r(i,j) = sqrt(x(i)^2+y(j)^2);
        rr = r(i,j);
        retemp = interp1(xcol,re,rr,'pchip');
        imtemp = interp1(xcol,im,rr,'pchip'); 
        %initial trial wavefunction
        psi0(i,j)= retemp + imtemp*1i;
        %define external potential
        v2d(i,j) = x(i)*x(i) + y(j)*y(j);%V2D
        %-km^2 on grid
        km2(i,j)=(xi(i).^2)+(xi(j).^2);
    end
end

kin=complex(cos(km2.*dt*0.5),-sin(km2.*dt*0.5));
R = alpha*(1+tanh(10*(R0 - r)))/2;
psiini2=abs(psi0).^2;
t=0; %time
k=0; %time step index
%psi0 now is the value before this loop
%define the psi1 value after this loop
psi1=ones(N,N);
%define output file
fid = fopen('1dgr.dat','w');
snr0 = 70;
psi0_noise = awgn(psi0,snr0,'measured','dB'); 
for i=1:N
    for j=1:N
        if (r(i,j)>R0)
            psi0(i,j) = psi0_noise(i,j);
        end
    end
end
%TIME LOOP BEGIN
while (t<t_max)
    %evaluate the difference of psi between adjacent loop
    err=max(max(abs(psi0-psi1)));
    %calculate equation term
    psi02=abs(psi0).^2; 
    fftpsi0=fft2(psi0);
    grad=sum(sum(km2.*abs(fftpsi0).^2))*mul/sum(sum(psi02));%fft.2
    pot=sum(sum(v2d.*psi02))/sum(sum(psi02));
    cubic=sum(sum(psi02.^2))/sum(sum(psi02));
    %output variables
    energy=grad+pot+0.5*cubic;
    chem=grad+pot+cubic;
    mass=sum(sum(psi02))/sum(sum(psiini2));
%     %%%%%%%%%%%%%renew the output file
%     fprintf(fid,' %12.8f  %12.8f  %12.8f  %12.8f %12.8f\n',t,energy,chem,mass,err);
     fprintf('%12.6e %12.6e %12.6e\n',t,chem,mass);
    %time step forward
    k=k+1;
    t=t+dt;
    %%%%%%%%%%%%%%%%%%time splitting
    psi1=psi0;
    %%%first half step V/2
    vpot=v2d+abs(psi1).^2;
    totpot=complex(cos(vpot.*dt*0.5),-sin(vpot.*dt*0.5));
    psi1=totpot.*psi1;
    %%%first half step T/2
    psi1=ifft2(kin.*fft2(psi1));
    %%%middle pump&decay terms
    for i=1:N
        for j=1:N
            if (R(i,j)~=0)
                HDP(i,j) = exp(R(i,j)*dt)*sqrt( R(i,j) / (R(i,j) + (exp(2*R(i,j)*dt) - 1)*sig*abs(psi1(i,j))^2) );
            else
                HDP(i,j) = sqrt( 1.0 / (1.0 + 2.0*sig*dt*abs(psi1(i,j))^2) );
            end
        end
    end
    psi1=HDP.*psi1;
    %%%last half step T/2
    psi1=ifft2(kin.*fft2(psi1));
    %%%last half step V/2
    vpot=v2d+abs(psi1).^2;
    totpot=complex(cos(vpot.*dt*0.5),-sin(vpot.*dt*0.5));
    psi1=totpot.*psi1;
    %renew psi0
    psi0=psi1;%./sqrt(sum(sum(abs(psi1).^2))*h^2);
    %subplot(2,1,1)
    psi02=abs(psi0).^2;
    mesh(xx, yy,psi02 )
    xlabel('x');
    ylabel('y');
   % axis([-8 8 -8 8 0 0.025])
    view(2)
%     subplot(2,1,2)
%     mesh(xxi, yyi, abs(fft2(psi0)).^2)
%     axis([-8 8 -8 8 0 0.02])
%     view(2)
    pause(0.01);
end %WHILE
fclose(fid);