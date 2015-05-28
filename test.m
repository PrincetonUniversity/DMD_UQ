clear all
close all

Ac = [1 -2; 1 -1];
dt = 0.1;
A = expm(Ac*dt);
EigsTrue = eig(A);

x0 = [ 1; 0.1];
%%
Nsteps = 100;
x = zeros(2,Nsteps);

x(:,1) = x0;
%Attain noise-free data
for kk = 1:(Nsteps-1)
    x(:,kk+1) = A*x(:,kk);
end

r = 2;
[Phi, Lambda, U, S, V, Atilde] = DMDext(x,r);

%
Ntrials = 100;
s = 0.2; %noise level

ErrorTrueVec = zeros(4,Ntrials);
ErrorEstVec = zeros(4,Ntrials);

ErrorTrueNorm=zeros(1,Ntrials);

EigsNaive = zeros(2,8,Ntrials);
EigsTLS = zeros(2,8,Ntrials);
% run DMD and NCDMD with different seeds of noise

NoiseSampleVec = [1 0; 0 1; 1 1;2 0; 0 2; 1 2; 2 1; 2 2];



for qq = 1:8;
    
    for nn = 1:Ntrials
        nn
        Noise = s*randn(size(x));
        
        Noise(1,:) = NoiseSampleVec(qq,1).*Noise(1,:);
        Noise(2,:) = NoiseSampleVec(qq,2).*Noise(2,:);
        
        xn = x + Noise;
        
        [Phin, Lambdan, U, S, V, Atilden] = DMDext(xn,r);
        An = U*Atilden*U';
        EigsNaive(:,qq,nn) = sum(Lambdan);
        ErrorTrueVec(:,nn) = A(:)-An(:);
        
        ErrorTrueNorm(nn) = norm(A-An);
        
        Admdtls  = dmd_tls(xn);
        EigsTLS(:,qq,nn) = sort(eig(Admdtls));
        
        
    end
end

%% Time check
%
% Noise = s*randn(size(x));
%
% xn = x + Noise;
% tic
% Admdtls  = tlsdmd(xn);
% toc
%
% tic
% Admdtls  = tlsdmdFast(xn);
% toc
% Fast is a small amount faster


%%
figure
plot((eig(A)),'b^','MarkerFaceColor','b','MarkerSize',12)
hold on
%plot(-1,-1,'k.')

%plot(-1,-1,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',12)
set(gca,'FontSize',14)

for qq = 1:8
    subplot(4,2,qq)
    plot(real((squeeze(EigsNaive(:,qq,:)))),imag((squeeze(EigsNaive(:,qq,:)))),'k.')
    hold on
    plot(real((squeeze(EigsTLS(:,qq,:)))),imag((squeeze(EigsTLS(:,qq,:)))),'r.')
   % plot(mean((squeeze(EigsNaive(:,qq,:))),2),'o','MarkerSize',12)

ylabel('Im(\lambda_c)')%,'Interpreter','Latex')
xlabel('Re(\lambda_c)')%,'Interpreter','Latex')
end
% For cts time
%legend('DMD, all trials','DMD, mean','NCDMD, all trials',...
%    'ncDMD, mean','True eigenvalue','Unit Circle','Location','NorthWest')

legend('True eigenvalue','DMD, all trials','DMD, mean','Location','NorthWest')