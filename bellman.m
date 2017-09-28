% Bellman equation v(k) = max u(f(k)-k')+Bv(k')
% log utility u=log(c) and f(k)=Ak^alpha+(a-delta)k
% Find the value and policy functions

% NOTE SEP 2017: WILL ADD MORE DOCUMENTATION LATER...

clear all;
clc;

% Set Variables
ALPHA    = (1/3);
BETA     = 0.75;
DELTA    = 1.00;
THETA    = 1.00;                
A        = 10;
KBAR     = ((A.*ALPHA)./((1/BETA)-(1-DELTA))).^(1/(1-ALPHA));

% Setting up the grid
NK       = 50;
ITER     = 0;
XMIN     = 2*KBAR/NK;          
XMAX     = 2*KBAR;
X        = (XMIN:(XMAX-XMIN)/(NK-1):XMAX)';

% Analytical values of policy and value functions
realV    = 4*((4/3)*log(10)+log(3/4)+(1/3)*log(1/4))+(4/9)*log(X);
realG    = 2.5*X.^ALPHA;

% Initialize a few vectors
v1       = zeros(NK,1);
g1       = zeros(NK,1);
f0       = zeros(NK,1);

% Convergence tolerance
DIFF     = 1;
EPS      = 1e-6;

while DIFF>EPS
    ITER = ITER + 1;
    
    for i=1:NK
        % Prevent wild values in gamma
        gamma = X <= A*X(i)^ALPHA+(1-DELTA)*X(i);  
        kp    = X(gamma);                           
        kt    = ones(length(kp),1)*X(i);
        
        c     = A*kt.^ALPHA+(1-DELTA)*kt-kp;        
        u     = log(c);                             
        
        f1    = u + BETA*f0(gamma);                 
        [v1(i),g1(i)] = max(f1);
    end
    gk = X(g1);
    DIFF = max(abs(v1-f0));
    meanresidv(ITER) = mean(abs(real(realV-v1)));
    meanresidg(ITER) = mean(abs(real(realG-gk)));
    f0 = v1;
end

h=figure(1);
subplot(2,1,1);
plot(X,f0,X,realV,'--'), hold on;
xlabel ('k'), ylabel ('Value Function at k');
grid on, grid minor;
title(sprintf('1.(b) Value Function v(k) (N=%d, delta=%.1f, mean resid=%.5f)', NK, DELTA, mean(meanresidv)));

subplot(2,1,2);
plot(X,gk,X,realG,'--'), hold on;
xlabel ('k'), ylabel ('Policy Function at k');
grid on, grid minor;
title(sprintf('1.(b) Policy Function g(k) (N=%d, delta=%.1f, mean resid=%.5f)', NK, DELTA, mean(meanresidg)));

saveas(h, '1_b_50.pdf');

% N = 500.

NK       = 500;                 
ITER     = 0;                   

XMIN     = 2*KBAR/NK;           
XMAX     = 2*KBAR;
X        = (XMIN:(XMAX-XMIN)/(NK-1):XMAX)';

realV    = 4*((4/3)*log(10)+log(3/4)+(1/3)*log(1/4))+(4/9)*log(X);
realG    = 2.5*X.^ALPHA;

v1       = zeros(NK,1);
g1       = zeros(NK,1);
f0       = zeros(NK,1);

DIFF     = 1;

while DIFF>EPS
    ITER = ITER + 1;
    
    for i=1:NK
        gamma = X <= A*X(i)^ALPHA+(1-DELTA)*X(i);   
        kp    = X(gamma);                           
        kt    = ones(length(kp),1)*X(i);
        
        c     = A*kt.^ALPHA+(1-DELTA)*kt-kp;        
        u     = log(c);                             
        
        f1    = u + BETA*f0(gamma);                 
        [v1(i),g1(i)] = max(f1);
    end
    gk = X(g1);
    DIFF = max(abs(v1-f0));
    meanresidv(ITER) = mean(abs(real(realV-v1)));
    meanresidg(ITER) = mean(abs(real(realG-gk)));
    f0 = v1;
end

h=figure(2);
subplot(2,1,1);
plot(X,f0,X,realV,'--'), hold on;
xlabel ('k'), ylabel ('Value Function at k');
grid on, grid minor;
title(sprintf('1.(b) Value Function v(k) (N=%d, delta=%.1f, mean resid=%.5f)', NK, DELTA, mean(meanresidv)));

subplot(2,1,2);
plot(X,gk,X,realG,'--'), hold on;
xlabel ('k'), ylabel ('Policy Function at k');
grid on, grid minor;
title(sprintf('1.(b) Policy Function g(k) (N=%d, delta=%.1f, mean resid=%.5f)', NK, DELTA, mean(meanresidg)));

saveas(h, '1_b_500.pdf');

% delta = 0.1.

DELTA    = 0.1;
DIFF     = 1;
v1       = zeros(NK,1);
g1       = zeros(NK,1);
f0       = zeros(NK,1);
ITER     = 0;                   
XMIN     = 2*KBAR/NK;           
XMAX     = 2*KBAR;
X        = (XMIN:(XMAX-XMIN)/(NK-1):XMAX)';

while DIFF>EPS
    ITER = ITER + 1;
    
    for i=1:NK
        gamma =X <= A*X(i)^ALPHA+(1-DELTA)*X(i);
        kp = X(gamma);
        kt = ones(length(kp),1)*X(i);
        
        c = A*kt.^ALPHA+(1-DELTA)*kt-kp;
        u = log(c);
        
        f1 = u+BETA*f0(gamma);
        [v1(i),g1(i)] = max(f1);
    end
    
    DIFF = max(abs(v1-f0));
    f0 = v1;
end

h=figure(3);
subplot(2,1,1);
plot(X,f0), hold on;
xlabel ('k'), ylabel ('Value Function at k');
grid on, grid minor;
title(sprintf('1.(c) Value Function v(k) (N=%d, delta=%.1f)', NK, DELTA));

subplot(2,1,2);
plot(X,gk), hold on;
xlabel ('k'), ylabel ('Policy Function at k');
grid on, grid minor;
title(sprintf('1.(c) Policy Function g(k) (N=%d, delta=%.1f)', NK, DELTA));

saveas(h, '1_c_500.pdf');