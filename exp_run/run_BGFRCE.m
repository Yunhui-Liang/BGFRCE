%
%
%

clear;
clc;

data_path = fullfile(pwd, '..',  filesep, "data_BPs", filesep);
addpath(data_path);
lib_path = fullfile(pwd, '..',  filesep, "lib", filesep);
addpath(lib_path);
code_path = genpath(fullfile(pwd, '..',  filesep, 'BGFRCE_ICASSP_2026'));
addpath(code_path);


dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};


exp_n = 'BGFRCE_RES';
for i1 =1 : length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    dir_name = [pwd, filesep, exp_n, filesep, data_name];
    create_dir(dir_name);
    clear BPs Y;
    load(strcat(data_path, datasetCandi{i1}));
    assert(size(BPs, 1) == size(Y, 1));
    nSmp = size(BPs, 1);
    nCluster = length(unique(Y));
    
    nBase = 20;

    Order = (2:1:6);
    Knn = (5:5:20);
    nOrder = length(Order);
    nK = length(Knn);
    %*********************************************************************
    % BGFRCE_RES
    %*********************************************************************
    fname2 = fullfile(dir_name, [data_name, '_', exp_n, '.mat']);
    if ~exist(fname2, 'file')
        nRepeat = 2;
        
        seed = 2024;
        rng(seed, 'twister')
        
        % Generate 50 random seeds
        random_seeds = randi([0, 1000000], 1, nRepeat * nRepeat );
        
        % Store the original state of the random number generator
        original_rng_state = rng;
        
        nMeasure = 15;
        
        BGFRCE_RES_result = zeros(nOrder, nK, nMeasure);
        BGFRCE_RES_result_std = zeros(nOrder, nK, nMeasure);
        BGFRCE_RES_time = zeros(nOrder, nK, 1);

        for iOrder = 1:nOrder
            for i2 = 1:nK
                res_1 = zeros(nRepeat, nMeasure);
                for iRepeat = 1:nRepeat
                    
                    idx = (iRepeat - 1) * nBase + 1 : iRepeat * nBase;
                    BPi = BPs(:, idx);
                    
                    t1_s = tic;
                    value_k = Knn(i2);
                    order = Order(iOrder);
                    [Hc, ~, Hcell] = compute_Hc_Hcell(BPi);
                    label_Y = kmeanspp(Hc', nCluster);
                    label = BGFRCE(Hc, Hcell, label_Y', nCluster, value_k, order);
                    
                    t1 = toc(t1_s);
                    result_10 = my_eval_y(label, Y);
                    res_1(iRepeat, :) = [result_10', t1];
                    
                end
                BGFRCE_RES_result(iOrder, i2, :) = mean(res_1, 1);
                BGFRCE_RES_result_std(iOrder, i2, :) = std(res_1, 1); 
            end
        end

        ACC = BGFRCE_RES_result(:,:,1);          
        [maxACC, linearIdx] = max(ACC(:));       
        [iacc, jacc] = ind2sub(size(ACC), linearIdx);     

        BGFRCE_RES_result_summary = squeeze(BGFRCE_RES_result(iacc,jacc, :))';
        BGFRCE_RES_result_summary_std = squeeze(BGFRCE_RES_result_std(iacc,jacc, :))';
        save(fname2, 'BGFRCE_RES_result', 'BGFRCE_RES_result_summary', 'BGFRCE_RES_result_std', 'BGFRCE_RES_result_summary_std');
        
        disp([data_name, ' has been completed!']);
    end
end
rmpath(data_path);
rmpath(lib_path);
rmpath(code_path);