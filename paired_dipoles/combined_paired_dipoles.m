%% Set/prepare folder structure
myelin_sim_dir = '<SIM_DIR>/myelin_sim_paired/';
myelin_sim_dir_single_dipoles = '<SIM_DIR>/myelin_sim/'; % to get the vertices for scaling the leadfields

if ~ isfolder(myelin_sim_dir)
    mkdir(myelin_sim_dir)
end

cd(myelin_sim_dir)
if ~ isfolder('results')
    mkdir('results')
end
cd('results')
if ~ isfolder('paired_dipoles')
    mkdir('paired_dipoles')
end
cd ..

% we use the single dipole simulations as a starting point (the leadfield matrix, patch locations,..)
copyfile([myelin_sim_dir_single_dipoles '/Sim_MEG_data_coreg_err_0.*'],myelin_sim_dir,'f')
copyfile([myelin_sim_dir_single_dipoles '/fixed_patches.mat'],myelin_sim_dir,'f')
copyfile([myelin_sim_dir_single_dipoles '/patch_loc_ind_60_bilateral.mat'],myelin_sim_dir,'f')
copyfile([myelin_sim_dir_single_dipoles '/U.mat'],myelin_sim_dir,'f')
copyfile([myelin_sim_dir_single_dipoles '/SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat'],[myelin_sim_dir 'SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat'],'f')

%% Set simulation parameters
dipole_mom = [20 6]; % dipole strength in nAm and dipole patch size
SNR = [0, -5, -10, -15, -20];% SNR in dB
source_recon_all = {'IID','EBB','MSP'}; % source priors Bayesian source reconstruction
rep_idx = 1:30; % simulated single dipole sources for each set

scale = [1 3/2 2/3 4/3 3/4 11/10 10/11 21/20 20/21]; % scaling factors
save('scaling','scale')
source_recon_idx = 1:3; % change if you want to perform simulations for specific source reconstruction appraoches only

%% Generate single dipole sources across 60 dipole locations
orig_data = [myelin_sim_dir '/Sim_MEG_data_coreg_err_0.mat'];
D = spm_eeg_load(orig_data);
D.inv{1}.gainmat = 'SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat';
D.save
load('patch_loc_ind_60_bilateral.mat')

for a = source_recon_idx
    source_recon_all_a = source_recon_all{a};
    for setNr = 1:2 % two sets
        if setNr == 1
            patch_loc_ind_sources = 1:30;
        elseif setNr == 2
            patch_loc_ind_sources = (1:30) + 30;
        end
        for sim_i = 1:7 % sim_scale
            for j = rep_idx
                % simulate paired dipoles
                vert_idx = Ip(patch_loc_ind_sources(j)); % current source location index
                loc_MNI = D.inv{1}.mesh.tess_mni.vert(vert_idx,:); % corresponding location
                if setNr == 1
                    vert_idx_2 = Ip(patch_loc_ind_sources(j)+30);
                elseif setNr == 2
                    vert_idx_2 = Ip(patch_loc_ind_sources(j)-30);
                end
                loc_MNI_2 = D.inv{1}.mesh.tess_mni.vert(vert_idx_2,:);

                save(sprintf('rep%d_params',j),'loc_MNI*', 'vert_idx*')

                for s = 1:length(SNR)
                    if SNR(s)==0
                        SNR_str = num2str(abs(SNR(s)));
                    else
                        SNR_str = ['minus' num2str(abs(SNR(s)))];
                    end

                    load('SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat','G','label')
                    load(sprintf('%srep%d_params_%d',myelin_sim_dir_single_dipoles,j,setNr),'weighted_x','x')
                    G(:,x) = G(:,x)+G(:,x).*(scale(sim_i)-scale(1)).*repmat(full(weighted_x)',size(G,1),1)*length(x);
                    load(sprintf('%srep%d_params_%d',myelin_sim_dir_single_dipoles,j,mod(setNr,2)+1),'weighted_x','x')
                    G(:,x) = G(:,x)+G(:,x).*(1./scale(sim_i)-scale(1)).*repmat(full(weighted_x)',size(G,1),1)*length(x);
                    save(sprintf('adapted_SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat'),'G','label')

                    D = spm_eeg_load(orig_data);
                    D.inv{1}.gainmat = 'adapted_SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat';
                    D.save

                    matlabbatch{1}.spm.meeg.source.simulate.D = {orig_data};
                    matlabbatch{1}.spm.meeg.source.simulate.val = 1;
                    matlabbatch{1}.spm.meeg.source.simulate.prefix = 'rep_';
                    matlabbatch{1}.spm.meeg.source.simulate.whatconditions.condlabel = {'Stim'};
                    matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = [100 400];
                    matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.fband = [10 30];
                    matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [dipole_mom; dipole_mom];
                    matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = [loc_MNI; loc_MNI_2];
                    matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR(s);
                    spm_jobman('run', matlabbatch);
                    clear matlabbatch
                    close all

                    for i = 1:length(scale) % do source reconstruction across scaling parameters
                        load('SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat')
                        load(sprintf('%srep%d_params_%d',myelin_sim_dir_single_dipoles,j,setNr),'weighted_x','x')
                        G(:,x) = G(:,x)+G(:,x).*(scale(i)-scale(1)).*repmat(full(weighted_x)',size(G,1),1)*length(x);
                        load(sprintf('%srep%d_params_%d',myelin_sim_dir_single_dipoles,j,mod(setNr,2)+1),'weighted_x','x')
                        G(:,x) = G(:,x)+G(:,x).*(1./scale(i)-scale(1)).*repmat(full(weighted_x)',size(G,1),1)*length(x);
                        save(sprintf('adapted_SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat'),'G','label')

                        rep_data = [myelin_sim_dir '/rep_Sim_MEG_data_coreg_err_0.mat']; % output file with simulated single dipole data from output from Simulation_Batch_singleDipole_job.m
                        D = spm_eeg_load(rep_data);
                        D.inv{1}.gainmat = 'adapted_SPMgainmatrix_Sim_MEG_data_coreg_err_0_1.mat';
                        D.save

                        matlabbatch{1}.spm.meeg.source.invertiter.D = {rep_data};
                        matlabbatch{1}.spm.meeg.source.invertiter.val = 1;
                        matlabbatch{1}.spm.meeg.source.invertiter.whatconditions.condlabel = {'Stim'};
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invtype = source_recon_all_a;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.woi = [100 400];
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.foi = [0 80];
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.hanning = 1;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {[myelin_sim_dir '/fixed_patches.mat']};
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = [1 Inf];
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = 0.6;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = 274;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.umodes = {[myelin_sim_dir '/U.mat']};
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 1;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
                        matlabbatch{1}.spm.meeg.source.invertiter.modality = {'All'};
                        spm_jobman('run', matlabbatch);
                        clear matlabbatch

                        D = spm_eeg_load(rep_data); % need to re-loaded

                        % store results
                        R2(s,j,i) = D.inv{1}.inverse.R2;
                        VE(s,j,i) = D.inv{1}.inverse.VE; % explained variance
                        F(s,j,i) = D.inv{1}.inverse.F; % negative Free Energy
                    end
                end
            end
            save([myelin_sim_dir '/results/paired_dipoles/DLE_model_ev_' source_recon_all_a '_sim_scale_' num2str(sim_i) '_' num2str(setNr)],'-regexp','^(?!(D)$).')
        end
    end
end
