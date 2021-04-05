clear all; clc;
%constants
E = 200e9;
A = 0.04;

%input vectors
position = [1,0,0; 2,0,1; 3,0.5,0.5; 4,1,0; 5,1,1; 6,1.5,0.5; 7,2,0; 8,2,1; 9,2.5,0.5; 10,3,0; 11,3,1];
L = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2;...
    sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2; sqrt(2)/2];
intconnmat = [1,1,2; 2,4,5; 3,7,8; 4,10,11; 5,1,4; 6,4,7; 7,7,10; 8,2,5; 9,5,8; 10,8,11; 11,1,3;...
        12,3,5; 13,4,6; 14,6,8; 15,7,9; 16,9,11; 17,2,3; 18,3,4; 19,5,6; 20,6,7; 21,8,9; 22,9,10;];
theta = [90; 90; 90; 90; 0; 0; 0; 0; 0; 0; 45; 45; 45; 45; 45; 45; -45; -45; -45; -45; -45; -45];
force = [1; 1; 1; 1; 0; 0; 0; -3000; 0; 0; 0; 0; 0; -3000; 0; 0; 0; 0; 0; -3000; 0; 0;];
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
 reducedgsm = glbstiffmat(5:end,5:end); %reducing the equation
 displacement(5:end,1) = reducedgsm^(-1)*force(5:end);
 force2 = glbstiffmat*displacement;
 
 for i=1:2:size(displacement,1)
     j=(i+1)/2;
    positionnew(j,2) = position(j,2)+displacement(i,1)*1e+5;
    positionnew(j,3) = position(j,3)+displacement(i+1,1)*1e+5;
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
     colour = [1 0 elemdisp(i)*2.87];
     plot(linevec(:,1),linevec(:,2),'--k','LineWidth',2);
     p = plot(linevec2(:,1),linevec2(:,2),'-');
     p.Color = colour;
     p.LineWidth = 5;
 end
  xlim([-0 4]);
  ylim([-1 2]);
  title('Prototype 1: Undeformed and Deformed shapes, scaling factor 1e+5');
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
  fprintf('Node 2: Reaction forces in y: %.3d N\n\n', force2(4));
  
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

  
     
 