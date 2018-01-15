function image_shuffler(condition)
% shuffles up all the images into different folders at the start of the
% experiment in order to randomize the shapes that represent each opponent.
%
% input: condition (1:4)
%
% must be called from the directory where the function is, i.e. the
% experiment directory

%%
% randomize generator seed
%--------------------------
rng('shuffle')

%----------------------------------------------------------------------
%         assigning conditions 
%----------------------------------------------------------------------
% opponent images
if condition == 1
    socmat = [1 0 1 0];
    socpo  = 'filled';
    nonsocpo= 'frame';
elseif condition == 2
    socmat = [0 1 0 1];
    socpo  = 'filled';
    nonsocpo= 'frame';
elseif condition == 3
    socmat = [1 0 1 0];
    socpo  = 'frame';
    nonsocpo= 'filled';
elseif condition == 0
    socmat = [0 1 0 1];
    socpo  = 'frame';
    nonsocpo= 'filled';
end


% opponent images
soc_img_dir     = dir(['images' filesep '*',socpo,'*']);
nonsoc_img_dir  = dir(['images' filesep '*',nonsocpo,'*']);

% the image directories are set up such that each shape appears twice, once
% as a frame and once filled. Each shape will be assigned to represent one
% opponent, and will never be repeated. So if a framed triangle is selected
% to represent a nonsocial opponent, a filled triangle will never be used
% for a social opponent
all_idx = 1:length(soc_img_dir); % both soc and nonsoc should be the same length
all_idx = Shuffle(all_idx);      % shuffle all indices so that each shape is only ever associated with 1 opponent

% first, delete all the previous images
for ii = 1:4
 
    tmp_dir = dir(['block', num2str(ii), '_images' filesep '*png*']);
    for tt = 1:length(tmp_dir)
        delete(['block', num2str(ii), '_images', filesep, tmp_dir(tt).name])
    end

end


block = 1;
for ii = 1:length(all_idx)
    social = socmat(block);
    if social == 1
        imstruct = soc_img_dir;
    else
        imstruct = nonsoc_img_dir;
    end
    

    copyfile(['images', filesep, imstruct(all_idx(ii)).name], ['block', num2str(block), '_images', filesep, imstruct(all_idx(ii)).name]);
    
    if mod(ii, 3) == 0
        block = block + 1;
    end
    
end
