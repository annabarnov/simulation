function display_setup(L, source_pos, mic_pos, fig, c)


floorX = [0 L(1) L(1) 0    0];
floorY = [0 0    L(2) L(2) 0];
floorZ = [0 0    0    0    0];

figure(fig)
hold on
% plot room:
plot3(floorX, floorY, floorZ,      'k');
plot3(floorX, floorY, floorZ+L(3), 'k');
plot3([0 0], [0 0], [0 L(3)], 'k');
plot3([L(1) L(1)], [0 0], [0 L(3)], 'k');
plot3([L(1) L(1)], [L(2) L(2)], [0 L(3)], 'k');
plot3([0 0], [L(2) L(2)], [0 L(3)], 'k');
axis equal; grid on

% plot source:
stem3(source_pos(:,1), source_pos(:,2), source_pos(:,3), 'r')
col = 'b';
if nargin == 5
    col = c;
end
stem3(mic_pos(:,1), mic_pos(:,2), mic_pos(:,3), col); hold off
title('Configuration')

end