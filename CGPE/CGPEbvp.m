
function y = CGPEbvp
%equation parameter
global alpha
global sig
global R0
global mu
global Rb

alpha = 4.4;
sig = 0.3;
mu = 3*alpha / (2*sig);
%real scale
Rb = 15.0;
R0 = 2.0;

% The differential equation set is seen in page 6 in https://arxiv.org/pdf/1310.2388.pdf
% four-point Lobatto method:4
xmesh = linspace(0,1,4);
solinit = bvpinit(xmesh,@guess);

options = bvpset('RelTol',1e-11,'NMax',40000,'Stats','on');


sol = bvp5c(@bvpfcn,@bcfun,solinit,options);

figure
sloy = abs(sol.y(1,:)).*abs(sol.y(1,:))+abs(sol.y(3,:)).*abs(sol.y(3,:));
slox = sol.x;
plot(Rb.*slox,sloy,'b');
y = sol;

%print the TF results
s=linspace(0,Rb,500);
psi0_temp=zeros(1,500);
for i = 1:500
    if (s(i)<sqrt(mu))
        psi0_temp(i)=sqrt(mu-s(i)^2);
    else
        psi0_temp(i)=0;
    end
end
hold on
plot(s(:),psi0_temp(:).^2)
axis([0 15 0 35])
%title('The densities of the steady states')
ylabel('|\phi(r)|^{2} ')
xlabel('r')
legend('数值解','初值')

end




function dydx = bvpfcn(x,y)
%ode form
%equation parameter
global alpha
global sig
global R0
global mu
global Rb

if (x==0)
    dydx = Rb*[ y(2)
                 0.5*(-y(6)*y(1) + (y(1)^2+y(3)^2)*y(1) - (alpha - sig*(y(1)^2+y(3)^2))*y(3))
                 y(4)
                 0.5*(-y(6)*y(3) + (y(1)^2+y(3)^2)*y(3) + (alpha - sig*(y(1)^2+y(3)^2))*y(1))
                 (alpha - sig*(y(1)^2+y(3)^2))*(y(1)^2+y(3)^2)*y(7)
                 0
                 1];
else
    dydx = Rb*[ y(2)
                 -y(2)/y(7) - y(6)*y(1) + y(7)^2*y(1) + (y(1)^2+y(3)^2)*y(1) - (alpha*(1+tanh(10*(R0 - y(7))))/2 - sig*(y(1)^2+y(3)^2))*y(3)
                 y(4)
                 -y(4)/y(7) - y(6)*y(3) + y(7)^2*y(3) + (y(1)^2+y(3)^2)*y(3) + (alpha*(1+tanh(10*(R0 - y(7))))/2 - sig*(y(1)^2+y(3)^2))*y(1)
                 (alpha*(1+tanh(10*(R0 - y(7))))/2 - sig*(y(1)^2+y(3)^2))*(y(1)^2+y(3)^2)*y(7)
                 0
                 1];
end
end
        


function res = bcfun(ya,yb)
res=[ya(2)
    yb(1)
    ya(4)
    ya(3)
    yb(3)
    yb(5)
    ya(7)];
end


function g = guess(x)
%equation parameter
global alpha
global sig
global R0
global mu
global Rb

f3 = @(r)(alpha.*(1+tanh(10.*(R0-Rb.*x)))./2-sig.*(mu-(r.*Rb).^2)).*(mu-(r.*Rb).^2).*r;

if (x*Rb<sqrt(mu))
g=[sqrt(2)/2*sqrt(mu-(x^2)*(Rb^2))                     %g=[sqrt(mu-(x^2)*(Rb^2))
   sqrt(2)/2*(-Rb*Rb*x/sqrt(mu-(x^2)*(Rb^2)))          %-Rb*Rb*x/sqrt(mu-(x^2)*(Rb^2))
   sqrt(2)/2*sqrt(mu-(x^2)*(Rb^2))                     % 0
   sqrt(2)/2*(-Rb*Rb*x/sqrt(mu-(x^2)*(Rb^2)))          % 0
   integral(f3,0,Rb*x,'AbsTol',1e-11)                      
   mu
   x];
else
  g=[ 0
    0
    0
    0
    integral(f3,0,sqrt(mu),'AbsTol',1e-11)
    mu
    x];
end
end



