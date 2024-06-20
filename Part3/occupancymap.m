% Criar um mapa ocupacional 3D com resolução de 0.1 metros
map = occupancyMap3D(0.1);

% Definir as zonas de exclusão como matrizes de pontos
[xNoFly1, yNoFly1, zNoFly1] = meshgrid(pNoFly1(1,1):0.1:pNoFly1(1,2), pNoFly1(2,1):0.1:pNoFly1(2,2), pNoFly1(3,1):0.1:pNoFly1(3,2));
pNoFly1_points = [xNoFly1(:), yNoFly1(:), zNoFly1(:)];

[xNoFly2, yNoFly2, zNoFly2] = meshgrid(pNoFly2(1,1):0.1:pNoFly2(1,2), pNoFly2(2,1):0.1:pNoFly2(2,2), pNoFly2(3,1):0.1:pNoFly2(3,2));
pNoFly2_points = [xNoFly2(:), yNoFly2(:), zNoFly2(:)];

% Adicionar zonas de exclusão ao mapa
setOccupancy(map, pNoFly1_points, 1);
setOccupancy(map, pNoFly2_points, 1);

% Visualizar o mapa 3D
figure;
show(map);

% Definir as posições iniciais e objetivos (para visualização)
pStart = [0, 0, 0]; % Ponto inicial
pGoal1 = [0, 1.1, 0.6]; % Objetivo

% Adicionar as posições iniciais e objetivos ao gráfico
hold on;
scatter3(pStart(1), pStart(2), pStart(3), 'bo', 'filled', 'DisplayName', 'Início');
scatter3(pGoal1(1), pGoal1(2), pGoal1(3), 'kx', 'LineWidth', 2, 'DisplayName', 'Objetivo');

% Configurar o gráfico
xlabel('X [meters]');
ylabel('Y [meters]');
zlabel('Z [meters]');
title('Mapa 3D com Zonas de Exclusão');
grid on;
axis equal;

legend;
hold off;
