clear;clc;close all;
load PRC_aeroelastic_Heave.mat;flag  = 2;
global PS_x;
global grad_PS_x;
global Thet;
if flag == 1
    PS_x      = PS_Stiff;%PS_Heave;
    grad_PS_x = PS_Grad_Stiff;%PS_Grad_Heave;
    Thet      = phase_refined;
else
    PS_x      = PS_Heave;%PS_Heave;
    grad_PS_x = PS_Grad_Heave;%PS_Grad_Heave;
    Thet      = phase_refined;
end
global omega;
T         = 5.242880000000000;
omega     = 2*pi/T;
frac      = 0.9;
T1        = frac*T;
options   = odeset('RelTol',1e-12,'AbsTol',[1e-12 1e-12]);
eps       = 0.000001;
lambda    = 0.000001;

for j = 1:100 % Solving BVP using newton iteration
    j
    lambdap     = lambda + eps;
    yp0         = [0;lambdap];
    lambdam     = lambda - eps;
    ym0         = [0;lambdam];
    [~, yp]     = ode45(@spike_timing_without_psi,[0 T1],yp0);
    [t, ym]     = ode45(@spike_timing_without_psi,[0 T1],ym0);
    J           = (yp(end,1)-ym(end,1))/(2*eps)';%clear yp ym

    y0          = [0;lambda];
    [t, y]      = ode45(@spike_timing_without_psi,[0 T1],y0);
    F           = (y(end,1)-2*pi);
    dF          = J;
    dlambda     = dF\(-F)
    lambda      = lambda + dlambda;
    clear J
    if abs(dlambda) < 5*10^-8
        break;
    end
end
global tspan u
tspan           = 0:0.001:T1;%linspace(0,T1,4000);
y0              = [0;lambda];
[t, y]          = ode45(@spike_timing_without_psi,tspan,y0); %simulates euler-lagrange equation with correct value of lambda!
y_temp          = wrapTo2Pi(y(:,1));
u               = y(:,2).*interp1(Thet,PS_x,y_temp)./(2);
theta           = y(:,1);
lambda          = y(:,2);
theta_dot = zeros(length(u),1);
lambda_dot = zeros(length(u),1);
for i= 1:length(t)
    yy_temp = spike_timing_without_psi(t,y(i,:));
    theta_dot(i) = yy_temp(1);lambda_dot(i) = yy_temp(2);
end
u_dot = lambda_dot.*interp1(Thet,PS_x,y(:,1))/(2) + y(:,2).*interp1(Thet,grad_PS_x,y(:,1)).*theta_dot/2;
save('uamp090.mat','t','u','u_dot');
Kappa = sum(u.^2)+sum(y(:,2).*(theta_dot-omega*ones(length(theta_dot),1)-interp1(Thet,PS_x,y(:,1)).*u));

fs = 18;
theta_dot = theta_dot - theta_dot(1);
lambda_dot = lambda_dot - lambda_dot(1);

figure;subplot(221);
plot(t/T,u,'k-','Linewidth',2);
grid on;grid minor;
xlabel('$t/T$', 'Interpreter', 'LaTeX','Fontsize',fs)
ylabel('$u(t)$', 'Interpreter', 'LaTeX','Fontsize',fs);
set(gca,'Fontsize',fs,'Fontname','Times');
box off;
print('-depsc','ctrl1.eps');close;
figure;subplot(221);
plot(t/T,theta_dot,'k-','Linewidth',2);
grid on;grid minor;
xlabel('$t/T$', 'Interpreter', 'LaTeX','Fontsize',fs)
ylabel('$\dot{\theta}(t)-\omega_n$', 'Interpreter', 'LaTeX','Fontsize',fs);
set(gca,'Fontsize',fs,'Fontname','Times');
box off;
print('-depsc','ctrl2.eps');close;
figure;subplot(221);
plot(t/T,lambda_dot,'k-','Linewidth',2);
grid on;grid minor;
xlabel('$t/T$', 'Interpreter', 'LaTeX','Fontsize',fs)
ylabel('$\dot{\lambda}(t)$', 'Interpreter', 'LaTeX','Fontsize',fs);
set(gca,'Fontsize',fs,'Fontname','Times');
box off;
print('-depsc','ctrl3.eps');close;


function df=spike_timing_without_psi(~,y)
df=zeros(2,1);
global omega
global PS_x;
global grad_PS_x;
global Thet;

y_temp = wrapTo2Pi(y(1));
u    = y(2)*interp1(Thet,PS_x,y_temp)/(2);%prc(y(1))
df(1)= omega+interp1(Thet,PS_x,y_temp)*u;%prc(y(1))
df(2)=-u*y(2)*interp1(Thet,grad_PS_x,y_temp);%dprc(y(1))
end
