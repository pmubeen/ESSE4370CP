clear all; clc;
%constants
E = 69e9;
A = 0.0225;

%input vectors
position = [1,0,0; 2,0,0.5; 3,0,1.5; 4,0,2.5; 5,0,3; 6,0.5,1; 7,0.5,2; 8,1,0.5; 9,1,1; 10,1,1.5;...
    11,1,2; 12,1,2.5; 13,2,1; 14,2,1.5; 15,2,2; 16,3,1.5];
L = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 1; 1; 1; 1; 1; 1; 1; 1; 1; sqrt(1.25); sqrt(1.25);...
    sqrt(1.25); sqrt(1.25); sqrt(1.25); sqrt(1.25); sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2;...
    sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(1.25); sqrt(1.25);];
intconnmat = [1,1,2; 2,4,5; 3,11,12; 4,10,11; 5,9,10; 6,8,9; 7,14,15; 8,13,14; 9,3,4; 10,2,3; 11,2,8;...
    12,3,10; 13,4,12; 14,9,13; 15,10,14; 16,11,15; 17,14,16; 18,1,8; 19,8,13; 20,13,16; 21,5,12;...
    22,12,15; 23,15,16; 24,2,6; 25,6,10; 26,3,7; 27,7,12; 28,4,7; 29,7,10; 30,3,6; 31,6,8; 32,11,14;...
    33,9,14];
theta = [90; 90; 90; 90; 90; 90; 90; 90; 90; 90; 0; 0; 0; 0; 0; 0; 0; 26.6; 26.6; 26.6; -26.6; -26.6;...
    26.6; 45; 45; 45; 45; -45; -45; -45; -45; -26.6; 26.6];
force = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; -3000; 0; 0; 0; 0; 0; 0; 0; -3000;...
    0; 0; 0; -3000];
glbstiffmat = zeros(size(position,1)*2);
% %creating the global stiffness matrix
 for i=1:length(L)
     l = cos(theta(i));
     m = sin(theta(i));
     n1 = 2*intconnmat(i,2)-1; %first node of element
     n2 = 2*intconnmat(i,3)-1; %second node of element
     locstiffmat = (1/L(i))*[l^2, l*m, -l^2, -l*m; l*m, m^2, -l*m, -m^2; -l^2, -l*m, l^2, l*m; -l*m,...
         -m^2, l*m, m^2]; %elemental stiffness matrix
     %adding the elemental stiffness matrix values to the global
     %stiffness matrix
     %row 1 of local stiffness matrix added
     glbstiffmat(n1,n1) = glbstiffmat(n1,n1)+locstiffmat(1,1); % u-disp first node
     glbstiffmat(n1,n1+1) = glbstiffmat(n1,n1+1)+locstiffmat(1,2); % v-disp first node
     glbstiffmat(n1,n2) = glbstiffmat(n1,n2)+locstiffmat(1,3); %u-disp second node
     glbstiffmat(n1,n2+1) = glbstiffmat(n1,n2+1)+locstiffmat(1,4); %v- disp second node
     
     %row 2 of local stiffness matrix added
     glbstiffmat(n1+1,n1) = glbstiffmat(n1+1,n1)+locstiffmat(2,1);
     glbstiffmat(n1+1,n1+1)= glbstiffmat(n1+1,n1+1)+locstiffmat(2,2);
     glbstiffmat(n1+1,n2) = glbstiffmat(n1+1,n2)+locstiffmat(2,3);
     glbstiffmat(n1+1,n2+1) = glbstiffmat(n1+1,n2+1)+locstiffmat(2,4);
     
     %row 3 of local stiffness matrix added
     glbstiffmat(n2,n1) = glbstiffmat(n2,n1)+locstiffmat(3,1);
     glbstiffmat(n2,n1+1)= glbstiffmat(n2,n1+1)+locstiffmat(3,2);
     glbstiffmat(n2,n2) = glbstiffmat(n2,n2)+locstiffmat(3,3);
     glbstiffmat(n2,n2+1) = glbstiffmat(n2,n2+1)+locstiffmat(3,4);
     
     %row 4 of local stiffness matrix added
     glbstiffmat(n2+1,n1) = glbstiffmat(n2+1,n1)+locstiffmat(4,1);
     glbstiffmat(n2+1,n1+1)= glbstiffmat(n2+1,n1+1)+locstiffmat(4,2);
     glbstiffmat(n2+1,n2) = glbstiffmat(n2+1,n2)+locstiffmat(4,3);
     glbstiffmat(n2+1,n2+1) = glbstiffmat(n2+1,n2+1)+locstiffmat(4,4);
 end
 
 glbstiffmat = E*A*glbstiffmat;
 displacement = zeros(size(position,1)*2,1);
 reducedgsm = glbstiffmat(11:end,11:end); %reducing the equation
 displacement(11:end,1) = reducedgsm^(-1)*force(11:end);
 force2 = glbstiffmat*displacement;
 
 for i=1:2:size(displacement,1)
     j=(i+1)/2;
    positionnew(j,2) = position(j,2)+displacement(i,1)*5e+4;
    positionnew(j,3) = position(j,3)+displacement(i+1,1)*5e+4;
 end
 
 figure(1); clf;
 hold on;
 for i=1:size(intconnmat,1)
     n1 = intconnmat(i,2); %first node of element
     n2 = intconnmat(i,3); %second node of element
     linevec = [position(n1,2),position(n1,3); position(n2,2), position(n2,3)];
     linevec2 = [positionnew(n1,2),positionnew(n1,3); positionnew(n2,2), positionnew(n2,3)];
     elemdisp(i) = sqrt((displacement(n1)*1e+5-displacement(n2)*1e+5)^2+(displacement(n1+1)*1e+5-displacement(n2+1)*1e+5)^2);
     elemstress(i) = sqrt((force2(n1)-force2(n2))^2+(force2(n1+1)-force2(n2+1))^2);
     colour = [1 0 elemdisp(i)*5];
     plot(linevec(:,1),linevec(:,2),'--k','LineWidth',2);
     p = plot(linevec2(:,1),linevec2(:,2),'-');
     p.Color = colour;
     p.LineWidth = 2;
 end
  xlim([-1 4]);
  ylim([-1 4]);
  title('Prototype 3: Undeformed and Deformed shapes, scaling factor 5e+4');
  xlabel('x');
  ylabel('y');
  axis equal
  
  fprintf('Nodal Displacements:\nMaximum nodal displacement in y: %.1d m\n', max(displacement(2:2:end)));
  fprintf('Minimum nodal displacement in y: %.3d m\n', min(displacement(2:2:end)));
  fprintf('Maximum nodal displacement in x: %.3d m\n', max(displacement(1:2:end)));
  fprintf('Minimum nodal displacement in x: %.3d m\n\n', min(displacement(1:2:end)));
  
  fprintf('Nodal Reaction Forces:\nNode 1: Reaction forces in x: %.d N\n', force2(1));
  fprintf('Node 1: Reaction forces in y: %.3d N\n', force2(2));
  fprintf('Node 2: Reaction forces in x: %.3d N\n', force2(3));
  fprintf('Node 2: Reaction forces in y: %.3d N\n', force2(4));
  fprintf('Node 3: Reaction forces in x: %.3d N\n', force2(5));
  fprintf('Node 3: Reaction forces in y: %.3d N\n', force2(6)); 
  fprintf('Node 4: Reaction forces in x: %.3d N\n', force2(7));
  fprintf('Node 4: Reaction forces in y: %.3d N\n', force2(8));
  fprintf('Node 5: Reaction forces in x: %.3d N\n', force2(9));
  fprintf('Node 5: Reaction forces in y: %.3d N\n', force2(10));
  fprintf('Elemental Stresses:\n');
for i=1:length(elemstress)
  fprintf('Element %i stress: %.3d N\n',i,elemstress(i));
end

  fprintf('Maximum elemental stress: %.3d N\n', max(elemstress));
  fprintf('Minimum elemental stress: %.3d N\n', min(elemstress));

  fprintf('\nElemental Strains:\n');
for i=1:length(elemstress)
  fprintf('Element %i strain: %.3d m/m\n',i,elemstress(i)/E);
end

fprintf('Maximum elemental strain: %.3d m\n', max(elemstress/E));
fprintf('Minimum elemental strain: %.3d m\n', min(elemstress/E));