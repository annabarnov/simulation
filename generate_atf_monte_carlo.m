% setup
dtheta      = 2;
radius      = 0.8;
Navg        = 100;
M           = 9;
spacing     = 0.05;

c           = 340;
f           = 1e3; FS = 16e3;
lambda      = c/f;
room        = [4,4,4];
t60         = 0.5;
ref_order   = -1;

% rir parameters
if t60 == 0 
    rir.rev_time          = 0.2;
    rir.ref_order         = 0;
else
    rir.rev_time          = t60;
    rir.ref_order         = ref_order;
end

rir.hp_filter         = 1;
rir.ns                = rir.rev_time*FS;
rir.type              = 'omnidirectional';
rir.orientation       = 0;

% init structure
atf = zeros(M, rir.ns, Navg);
hd5_filename = ['orig_atf_' num2str(t60*1000) '.h5'];




% for each theta
theta_arr   = (0:dtheta:360)*pi/180; theta_arr = theta_arr(1:end-1);

h5create(hd5_filename,'/dataset',[M, rir.ns, Navg, length(theta_arr)], 'ChunkSize',[M, rir.ns, Navg, 1])

for j = 1:length(theta_arr)
    tic
    fprintf('theta %d of %d\n', j, length(theta_arr))
    theta = theta_arr(j);
    speaker_rel_loc = [radius*cos(theta), radius*sin(theta), 0];
    % monte carlo
    parfor i = 1:Navg
        fprintf('itr#%d\n',i)
        rotation = 0 + (360)*rand(1,1);
        bound = lambda/2 + radius + spacing;
        a = 0 + bound; b = room(1) - bound; array_loc_x = a + (b-a)*rand(1,1);
        a = 0 + bound; b = room(2) - bound; array_loc_y = a + (b-a)*rand(1,1);
        a = 0 + lambda/2; b = room(2) - lambda/2; array_loc_z = a + (b-a)*rand(1,1);
        array_loc = [array_loc_x, array_loc_y, array_loc_z];


        [mic_pos, source_pos]  = p_absolute_position_from_relative(M, array_loc, spacing, speaker_rel_loc, rotation);
%         display_setup(room, source_pos, mic_pos, 10);
        atf(:,:,i)      = rir_generator(c, FS, mic_pos, source_pos(1,:), room, rir.rev_time, rir.ns, rir.type, rir.ref_order, 3, rir.orientation, rir.hp_filter);
%         atf(:,:,i)      = rchart_new(c, FS, mic_pos, source_pos(1,:), room, rir.rev_time, rir.ns, rir.type, rir.ref_order, 3, rir.orientation, rir.hp_filter); 
    end
    
    h5write(hd5_filename, '/dataset', atf, [1, 1, 1, j], [M, rir.ns, Navg 1]);

    toc
end