function cd = wing_fitness(x)
% x(1) = sweep
% x(2) = twist

Re = 1e7;
alpha = 2; %In deg!!
%cr =1; %In m
%ct=1; %In m
b=6; %In m
AR = 15;
m=5;
n=4;
m_piu1 = m + 1;
n_piu1 = n + 1;
cr = 2*b/AR;
ct = cr;
NACA = 0012;

%Generazione geometria
Wing = build_wing(cr,ct,b,m_piu1,n_piu1, NACA, 0, 0);


%% code
W_let = build_winglet(ct,x(1),x(2),0.01,m+1,25,deg2rad(x(3)),0,deg2rad(x(4)),NACA,Wing,0,x(5));

Wing = assemble_wing(Wing,W_let);

%Wing2(:,:,1) = Wing(:,:,1);
%Wing2(:,:,2) = Wing(:,:,2);
%Wing2(:,:,3) = -Wing(:,:,3) + 2*max(max(Wing(:,:,3)));
%Wing2 = Wing2(:,end:-1:1,:);
%Wing = assemble_wing(Wing,Wing2);

m = size(Wing,1)-1;
n = size(Wing,2)-1;



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

Esecuzione = system('./MAIN');

load Coeff
cd = Coeff(2);



end