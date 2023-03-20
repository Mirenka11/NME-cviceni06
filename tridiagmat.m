%%
% Aby skript fungoval spravne, tak je nejdrive potreba doplnit 2 radky oznacene takto: !INSTRUKCE!
% Az skript bude fungovat, tak muzete menit parametry v radkach oznacenych takto: ?INSTRUKCE?
% Na konci skriptu jsou otazky k zodpovezeni

%% procedura pro reseni tridiagonalni matice:

clc % vycistim si pracovni prostredi
clear all
close all

N = 10000;  % pocet kroku, ?VYZKOUSEJTE RUZNE POCTY KROKU?, sledujte jak se meni cas vypoctu pro obe pouzite metody
h = 1/N; % krok

x = 0+h:h:1-h; % rozdeleni intervalu na 10 podintervalu, hodnoty x jsou vnitrni uzly (neobsahuje krajni body x=0 a x=1)

%r = ones(1,length(x));
%r = sin(pi*x);
%r = (x-1/2).^2;
%r=log(x);
%r=rand(10,1)*x;
 % hodnoty nabojove hustoty v uzlech intervalu
 % ?ZKUSTE RUZNE FUNKCE?, napr.: sin(pi*x), kvadratickou funkci (x-1/2).^2, atd.
epsilon0 = 1; % konstantu muzem dat =1, skalovani nas ted nezajima

a = ones(N-1,1); % diagonala
b = (-1/2)*ones(N-1,1); % nad diagonalou
c = (-1/2)*ones(N-1,1); % pod diagonalou
f = (h.^2.*r)/(2.*epsilon0); % prava strana !!!DOPLNENO!!!

phi   = zeros(N+1,1); % vektor reseni, obsahuje i krajni body x=0 a x=1

tic % odtud merim cas vypoctu
% pomocne vektory
rho = zeros(N+1,1); 
mu  = zeros(N+1,1); 

for i = 1 : N-1   % primy beh
    mu(i+1)  = -b(i) / ( c(i)*mu(i) + a(i) ); % doplnim do vzorecku podle teorie
    rho(i+1) = ( f(i) - c(i)*rho(i) ) / ( c(i)*mu(i) + a(i) );
end

for i = N-1 : -1 : 1   % zpetny beh
    phi(i+1) = mu(i+1)*phi(i+2) + rho(i+1);
end
time_tridiag = toc % zapisu cas vypoctu

% vykresleni reseni
figure
plot([0, x, 1],phi)

%%  reseni standartnim matlabovskym '\'
tic % odtud merim cas vypoctu
b(end)=[]; % pro sestaveni tridiagonalni matice potrebuju vektory 'b' a 'c' o jeden clen kratsi
c(end)=[];
A = diag(b,1)+diag(c,-1)+diag(a); %!SESTAVTE TRIDIAGONALNI MATICI! % sestavim si tridiagonalni matici A
phi2 = A\f'; % provedu vypocet phi
time_backslash = toc % zapisu cas vypoctu

% vykresleni reseni
figure
plot([0, x, 1],phi)
%% 
% OTAZKY:
%% jake jsou casy vypoctu kazde metody pro 10 tis. kroku?
%   * pro r = ones
%       nase funkce:    0.0075
%       matlabovske \:  16.7951     
%   * pro r = sin(pi*x)
%       nase funkce:    0.0102
%       matlabovske \:  16.0390
%   * pro r = (x-1/2).^2
%       nase funkce:    0.0100
%       matlabovske \:  16.0813 
%   * pro r = log(x)
%       nase funkce:    0.0472
%       matlabovske \:  17.0689
%   * pro nahodnou distribuci
%       nase funkce:    0.0085
%       matlabovske \:  16.2115

%% na kterem radku jsme zadali okrajovou podminku phi(x=0)=0 a phi(x=1)=0?
%   phi   = zeros(N+1,1); radek 31