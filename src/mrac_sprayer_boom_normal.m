%simulation for sbiagro2023 article
% based on Cui et al., 2019 Hydraulic-drive roll movement control of a spray boom using adaptive robust
% control strategy

clear;

pkg load control;
pkg load signal;

output_precision(3); % 2 decimals

clc;
% Parameters - config - nominal system g1
c2 = -0.0058;
c1 = 67.29;
c0 = 174.2;
b1 = 0.834;
b0 = 1.993;

% reference model

%km = 0.8; % m(s) = 0.8 / (s + 0.8)
%am1 = 0.8;

km = 1;
bm1 = 0.8; % 0.8 / (s^2 + 1,6s + 0,8)
am1 = 1.6;
am2 = 0.8;

%ms = tf([0 0 km], [1 am1]);
ms = tf([0 1 bm1], [1 am1 am2]);

[Am, Bm, Cm, Dm] = tf2ss(ms);

g1 = tf([c2 c1 c0],[1 b1 b0]);

dt = 0.01;
tMax = 80.00;

% system

kp = 1;

[A,B,C,D] = tf2ss(g1);


k = 2;

u(1) = 0; % para malha fechada

tplot(1) = 0;
X(:,1) = [0;0];
Xm(:,1) = [0; 0];

t = 0;

flag = false;

% controller parameters
lambda1 = .8; %for IO filters v1,v2
gamma = 0.1; %for g vector used for IO filter's calculation
g = gamma;
gamma1 = .1; % adaptative gain for \theta_v11
gamma2= .1; % adaptive gain for \theta_2
gamma3 = .5; % adaptative gain for \theta_v21
gamma4 = .1 ;% adaptive gain for \theta_4

%r = ones(1, tMax/dt + 2);
r(1) = 0;

e(1) = 0;

theta_v1(1) = 0;
theta_2(1) = 0;
theta_v2(1) = 0;
theta_4(1) = 0;

theta_v1_star = (0.8 - 2.588) / gamma;
theta_2_star = (0.834 - 1.6) / 67.29;
theta_v2_star = (1.993 - 0.8 - 67.29*theta_2_star*0.8) / (gamma*67.29);
theta_4_star = 1 / 67.29;

v1(1) = 0;
v2(1) = 0;

y(1) = 0;
ym(1) = 0;

tau = 0.1;

while (t < (tMax-dt))

   d = 3*sin(0.2*pi*t);
   dp(k) = d;

   dotX = A*X(:,k-1) + B*u(k-1);
   X(:,k) = X(:,k-1) + dt*(dotX);

   if (t < 15)
    y(:,k) = C*X(:,k);
   else
    y(:,k) = C*X(:,k) + d;
   endif

   %y(:,k) = C*X(:,k);
   dotXm = Am*Xm(:,k-1) + Bm*r(k-1); %sem is
   %r(k) = 1;
   r(k) = 3*sin(0.2 * pi * t);

   Xm(:,k) = Xm(:,k-1) + dt*(dotXm);
   ym(:,k) = Cm*Xm(:,k);

   dotv1 = -lambda1*v1(k-1) + g*u(k-1);
   v1(k) = v1(k-1) + dt*(dotv1);

   dotv2 = -lambda1*v2(k-1) + g*y(k);
   v2(k) = v2(k-1) + dt*(dotv2);

   e(k) = y(k) - ym(k);

   % Lyapunov-based laws - MRAC
   theta_v1(k) = theta_v1(k-1) + dt*(-gamma1 * e(k) * v1(k));
   theta_2(k) = theta_2(k-1) + dt*(-gamma2 * e(k) * y(k));
   theta_v2(k) = theta_v2(k-1) + dt*(-gamma3 * e(k) * v2(k));
   theta_4(k) = theta_4(k-1) + dt*(-gamma4 * e(k) * r(k));
################### FIM MRAC ###################

   % u é usado com MRAC
   u(k) = theta_v1(k)*v1(k) + theta_2(k)*y(k) + theta_v2(k)*v2(k) + theta_4(k)*r(k);

   t = t + dt;
   tplot(k) = t;
   k = k + 1;
endwhile


subplot(3,1,1);
plot(tplot, y, 'k-', 'LineWidth', 2, tplot, ym, 'b-', 'LineWidth', 2);
legend('y', 'ym');
title('�ngulo de rolagem (\beta) - MRAC sem Input Shaping');
ylabel('\beta (^o)');

subplot(3,1,2);
plot(tplot, e, 'k-', 'LineWidth', 2);
title('Erro de rastreamento ');
ylabel('^o');

subplot(3,1,3);
plot(tplot, u, 'k-', 'LineWidth', 2);
title('Sinal de controle - deslocamento do cilindro (m)');
xlabel('Tempo (s)');
ylabel('l_d (m)');

% Métricas

%IAE

iae = sum(abs(e));
itae = sum(t .* abs(e));
ise = sum(e.^2);
tu = sum(abs(u));
stdu = std(u);

printf('iae: %.2f\n', iae);
printf('itae: %.2f\n', itae);
printf('ise: %.2f\n', ise);
printf('total u: %.2f', tu);
printf('std u: %.2f', stdu);




