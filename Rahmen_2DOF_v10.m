% Skript zur dynamischen Berechnung eines Rahmens mit bis zu drei
% Geschossen infolge Anfangsverschiebungen und/oder -geschwindigkeiten
%
% Version v1.0, 25.04.2023
% E212-03, FoB Baumechanik und Baudynamik
% Grundlage ist das Script 'ML_29_1' aus:
% C. Petersen, H. Werkle, Dynamik der Baukonstruktionen
% 2. Auflage, Springer Vieweg, Wiesbaden, 2018
% 
% Ausgabe:
%   - Eigenfrequenzen und -formen inkl. vereinfachter Darstellung
%   - Verschiebung, Geschwindigkeit und Beschleunigung an den
%     Geschossdecken (numerische Integration mit Newmark-Verfahren)
%
% Hinweise:
%   - Höhe und Biegesteifigkeit aller Einzelstützen sind gleich
%   - Als Lagerungsart der Stützen kann zwischen gelenkig-eingespannt und
%     eingespannt-eingespannt gewählt werden
%   - Dämpfung ist in jedem Geschoss vorhanden (ANNAHME: modale Dämpfung)
%
clear; clc; close all
% ---------------------------- EINGABEBLOCK ------------------------------
% RAHMEN
% Masse je Geschoss in [kg]
m1 = 15000;         % 1. Geschoss
m2 = 10000;         % 2. Geschoss

% Massenmatrix 
M = [m1 0;
     0 m2];

% Stützen
H = 5;              % Höhe [m]
EI = 5*10^6;        % Biegesteifigkeit [Nm^2]
                       
% Äquivalente Federsteifigkeit einer Einzelstütze in [N/m]
kG = 3*EI/H^3;      % Lagerung: gelenkig-eigespannt 
kE = 12*EI/H^3;     % Lagerung: eingespannt-eingespannt

% Äquivalente Federsteifigkeit je Geschoss
k1 = 4*kG;          % 1. Geschoss
k2 = 2*kE;          % 2. Geschoss   

% Steifigkeitsmatrix
K = [k1+k2  -k2;
      -k2    k2];

% Dämpfung je Geschoss in [Ns/m]
zeta = 0.03;        % Lehr'sches Dämpfungsmaß [-]

c1 = 2*zeta*sqrt(m1*k1);
c2 = 2*zeta*sqrt(m2*k2);

% Dämpfungsmatrix 
C = [c1+c2  -c2;
      -c2    c2];

% ANFANGSWERTE
x_0 = [0.1 0];     % Anfangsverschiebungen [m]
v_0 = [0 0];       % Anfangsgeschwindigkeiten [m/s] 

% ZEIT
t_ber = 30;          % Berechnungszeit [s]
dt = 0.005;          % Berechnungszeitschritt [s]

% -------------------------- BERECHNUNGSBLOCK ----------------------------
% Lösung Eigenwertproblem
[EF,EW] = eig(K,M); % A-Eigenformmatrix, EW-Eigenwertmatrix

% Definition eines Vektors für die ermittelten Eigenwerte
D_EW = diag(EW);

% Ermittlung der Eigenkreisfrequenzen
Omega = sqrt(D_EW);

% Sortierung der Eigenwerte in aufsteigender Reihenfolge
[Omega,index] = sortrows(Omega);
EF = EF(:,index);

% Ermittlung der Eigenfrequenzen
Freq = Omega/(2*pi);

% Ermittlung der Eigenschwingzeiten
T = 1./Freq;

% Anzahl der Berechnungszeitschritte
nt = ceil(t_ber/dt)+1;
t_int = 0:dt:dt*(nt-1);

% Normierung der Eigenvektoren auf das betragsgrößte Element
n = size(K,1); % Matrixdimension
for j = 1:1:n
    if max(EF(:,j))>abs(min(EF(:,j)))
       z1 = max(EF(:,j));
       for i = 1:1:n
         EF(i,j) = (EF(i,j)/z1);
       end
    else
       z1 = min(EF(:,j));
       for i = 1:1:n
         EF(i,j) = (EF(i,j)/z1);
       end
    end
end

% Definition der Ergebnisvektoren
x = zeros(n,nt);  % Verschiebungsmatrix
v = zeros(n,nt);  % Geschwindigkeitsmatrix
a = zeros(n,nt);  % Beschleunigungsmatrix

% Berücksichtigung der Anfangsbedingungen
for i = 1:1:n
    x(i,1) = x_0(i);
    v(i,1) = v_0(i);
end

count_m0 = 0;  % Zähler der Nulleinträge auf der Diagonale der Massenmatrix
for ii = 1:1:size(M,1)    
    if M(ii,ii) == 0
        count_m0 = count_m0+1;
    end
end
if count_m0 == 0
    a(:,1) = (-M^-1)*(C*v(:,1)+K*x(:,1));
else
    a(:,1) = -pinv(M)*(C*v(:,1)+K*x(:,1));
end

% Integrationsparameter für das Newmark-Verfahren
alpha = 0.5;
beta = 0.25;

% Berechnung der Schwingreaktion mittels Newmark-Verfahren
for i = 2:1:nt
    a_h = ((1/beta)*M)+(alpha/beta)*C*dt+K*dt^2;
    b_h = ((1/beta)*M+(alpha/beta)*C*dt)*x(:,i-1)+((1/beta)*M+...
        (alpha/beta-1)*C*dt)*dt*v(:,i-1)+((1/(2*beta)-1)*M+...
        (alpha/(2*beta)-1)*C*dt)*dt^2*a(:,i-1);
    x(:,i) = a_h^-1*b_h;
    v(:,i) = (alpha/(beta*dt))*(x(:,i)-x(:,i-1))-((alpha/beta)-1)*...
        v(:,i-1)-(alpha/(2*beta)-1)*dt*a(:,i-1);
    a(:,i) = (1/(beta*dt^2))*(x(:,i)-x(:,i-1))-1/(beta*dt)*v(:,i-1)-...
        (1/(2*beta)-1)*a(:,i-1);
end

% Ausgabe der max. (Absolut-)Werte der Verschiebung je Geschoss
disp('Max. (Absolut-)Werte der Verschiebungen:')
disp(append('x1 = ',num2str(round(max(abs(x(1,:))),2)),' m'))
switch n
case 3
disp(append('x2 = ',num2str(round(max(abs(x(2,:))),2)),' m'))
disp(append('x3 = ',num2str(round(max(abs(x(3,:))),2)),' m'))
case 2
disp(append('x2 = ',num2str(round(max(abs(x(2,:))),2)),' m'))
end

% ------------------------- DARSTELLUNGSBLOCK ----------------------------
% Grafische Darstellung der Ergebnisse
name_fig1 = 'Eigenformen';
fig1 = figure('Name',name_fig1,'NumberTitle','off'); 
set(fig1,'Position',[200 500 1000 500]);
tl_EF = tiledlayout(1,n);
tl_EF.TileSpacing = 'compact'; tl_EF.Padding = 'loose';

% Koordinaten
X = [0 5];
switch n
case 3
    Y = [0 H 2*H 3*H];
case 2
    Y = [0 H 2*H];
case 1
    Y = [0 H];
end

for j = 1:1:n
nexttile
% DARSTELLUNG RAHMEN
% Stützen 1. Geschoss
line([X(1) X(1)],[Y(1) Y(2)],'linewidth',2)
line([X(2) X(2)],[Y(1) Y(2)],'linewidth',2)

% Geschossdecke 1. Geschoss
line([X(1) X(2)],[Y(2) Y(2)],'linewidth',10)

switch n
case 3
    % Stützen 2. Geschoss
    line([X(1) X(1)],[Y(2) Y(3)],'linewidth',2)
    line([X(2) X(2)],[Y(2) Y(3)],'linewidth',2)
    % Stützen 3. Geschoss
    line([X(1) X(1)],[Y(3) Y(4)],'linewidth',2)
    line([X(2) X(2)],[Y(3) Y(4)],'linewidth',2)
    
    % Geschossdecke 2. Geschoss
    line([X(1) X(2)],[Y(3) Y(3)],'linewidth',10)
    % Geschossdecke 3. Geschoss
    line([X(1) X(2)],[Y(4) Y(4)],'linewidth',10)
case 2
    % Stützen 2. Geschoss
    line([X(1) X(1)],[Y(2) Y(3)],'linewidth',2)
    line([X(2) X(2)],[Y(2) Y(3)],'linewidth',2)
    
    % Geschossdecke 2. Geschoss
    line([X(1) X(2)],[Y(3) Y(3)],'linewidth',10)
end

% DARSTELLUNG EIGENFORMEN
% Stützen 1. Geschoss
line([X(1) X(1)+EF(1,j)],[Y(1) Y(2)],'linewidth',2,'color','#D95319')
line([X(2) X(2)+EF(1,j)],[Y(1) Y(2)],'linewidth',2,'color','#D95319')

% Geschossdecke 1. Geschoss
line([X(1)+EF(1,j) X(2)+EF(1,j)],[Y(2) Y(2)],'linewidth',2,'color','#D95319')

switch n
case 3
% Stützen 2. Geschoss    
line([X(1)+EF(1,j) X(1)+EF(2,j)],[Y(2) Y(3)],'linewidth',2,'color','#D95319')
line([X(2)+EF(1,j) X(2)+EF(2,j)],[Y(2) Y(3)],'linewidth',2,'color','#D95319')
% Stützen 3. Geschoss
line([X(1)+EF(2,j) X(1)+EF(3,j)],[Y(3) Y(4)],'linewidth',2,'color','#D95319')
line([X(2)+EF(2,j) X(2)+EF(3,j)],[Y(3) Y(4)],'linewidth',2,'color','#D95319')

% Geschossdecke 2. Geschoss
line([X(1)+EF(2,j) X(2)+EF(2,j)],[Y(3) Y(3)],'linewidth',2,'color','#D95319')
% Geschossdecke 3. Geschoss
line([X(1)+EF(3,j) X(2)+EF(3,j)],[Y(4) Y(4)],'linewidth',2,'color','#D95319')

case 2
% Stützen 2. Geschoss 
line([X(1)+EF(1,j) X(1)+EF(2,j)],[Y(2) Y(3)],'linewidth',2,'color','#D95319')
line([X(2)+EF(1,j) X(2)+EF(2,j)],[Y(2) Y(3)],'linewidth',2,'color','#D95319')

% Geschossdecke 2. Geschoss
line([X(1)+EF(2,j) X(2)+EF(2,j)],[Y(3) Y(3)],'linewidth',2,'color','#D95319')
end

xlim([X(1)-2 X(2)+2]); ylim([0 Y(end)+1])
title(append(num2str(j),'. EF - \omega_',num2str(j),' = ',...
      num2str(round(Omega(j),2)),' rad/s'))
axis equal
end

name_fig2 = 'Schwingungsantworten';
fig2 = figure('Name',name_fig2,'NumberTitle','off'); 
set(fig2,'Position',[1300 300 1100 700]);

tl_SA = tiledlayout(3,n);
tl_SA.TileSpacing = 'compact'; tl_SA.Padding = 'loose';

for j = 1:1:n
nexttile
plot(t_int,x(j,:),'MarkerSize',3,'LineWidth',1)
title(append(num2str(j),'. Geschoss'))
xlabel('Zeit [s]')
ylabel('Verschiebung [m]')
ylim([-(max(abs(x),[],'all')+0.2) max(abs(x),[],'all')+0.2])
grid on
end

for j = 1:1:n
nexttile
plot(t_int,v(j,:),'MarkerSize',3,'LineWidth',1)
xlabel('Zeit [s]')
ylabel('Geschwindigkeit [m/s]')
ylim([-(max(abs(v),[],'all')+2) max(abs(v),[],'all')+2])
grid on
end

for j = 1:1:n
nexttile
plot(t_int,a(j,:),'MarkerSize',3,'LineWidth',1)
xlabel('Zeit [s]')
ylabel('Beschleunigung [m/s^2]')
ylim([-(max(abs(a),[],'all')+10) max(abs(a),[],'all')+10])
grid on
end


