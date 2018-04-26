%Creo matrici per script fortran
close all; clear;

Re = 1e7;
alpha = 2; %In deg!!
%cr =1; %In m
%ct=1; %In m
b=1; %In m
AR = 10;
m=18;
n=50;
m_piu1 = m + 1;
n_piu1 = n + 1;
cr = 2*b/AR;
ct = cr;
NACA = 0012;

%Generazione geometria
Wing = build_wing(cr,ct,b,m_piu1,n_piu1, NACA, 0, 0);


%% code
%W_let = build_winglet(ct,0.505,0.252,0.1,m+1,25,deg2rad(1.992),0,deg2rad(1.35),NACA,Wing,0,1);

%Wing = assemble_wing(Wing,W_let);

%Wing2(:,:,1) = Wing(:,:,1);
%Wing2(:,:,2) = Wing(:,:,2);
%Wing2(:,:,3) = -Wing(:,:,3) + 2*max(max(Wing(:,:,3)));
%Wing2 = Wing2(:,end:-1:1,:);
%Wing = assemble_wing(Wing,Wing2);

m = size(Wing,1)-1;
n = size(Wing,2)-1;



%Plot
surf(Wing(:,:,1),Wing(:,:,2),Wing(:,:,3))
xlabel('chord')
ylabel('span')
drawnow
%axis equal
x=Wing(:,:,1);
y=Wing(:,:,2);
z=Wing(:,:,3);

%Salvataggio geometria
save X.dat x -ascii
save Y.dat y -ascii
save Z.dat z -ascii

% file data.dat
id = fopen('Data.dat','w');
form = '%s=%d\n%s=%d\n%s=%f\n%s=%f\n%s=%f\n%s=%d\n';
fprintf(id,form,'m',m,'n',n,'b',b,'AR',AR,'alpha',alpha,'Re',Re);
fclose(id);

%Compilo e linko il programma

Compilazione = system('f95 -c Subroutines.f90 InFiles.f90 Main.f90');
Link = system('f95 -o MAIN Subroutines.f90 InFiles.f90 Main.f90 -L/usr/local/lib -llapack -lblas');

%Lancio il programma
tic
if Compilazione==0 && Link ==0
Esecuzione = system('./MAIN');
toc
load Coeff
Coeff
% 
% CL=[];
% CDi=[];
% ALPHA = -5:5;
% for alpha = ALPHA
%    %file Data.dat
%   id = fopen('Data.dat','w');
%   form = '%s=%d\n%s=%d\n%s=%f\n%s=%f\n%s=%f\n%s=%d\n';
%   fprintf(id,form,'m',m,'n',n,'b',b,'AR',AR,'alpha',alpha,'Re',Re);
%   fclose(id);
%   system('./MAIN');
%   load Coeff
%   CL=[CL; Coeff(1)];
%   CDi=[CDi; Coeff(2)];
% end
% 
% %%
% 
% figure(2)
% subplot(2,2,1)
% plot(ALPHA, CL,'o--','linewidth',2)
% xlabel('\alpha [deg]')
% ylabel('C_L')
% grid on
% %axis([ALPHA(1) ALPHA(end) -5 5])
% 
% figure(2)
% subplot(2,2,2)
% plot(ALPHA, CDi,'o--','linewidth',2)
% xlabel('\alpha [deg]')
% ylabel('C_D_i')
% grid on
% %axis([ALPHA(1) ALPHA(end) -5 5])
% %load v1; load v2; load v3; load v4;
% figure(2)
% subplot(2,2,3)
% plot(CDi,CL,'o--','linewidth',2)
% xlabel('C_D_i')
% ylabel('C_L')
% grid on
% 
% figure(2)
% subplot(2,2,4)
% plot(ALPHA, CL./CDi,'o--','linewidth',2)
% xlabel('\alpha')
% ylabel('E')
% grid on
end
