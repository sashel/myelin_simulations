function [] = plot_model_ev_paired_across_SNRs(data_dir,source_recon_idx, sim_i)
% use as plot_model_ev_paired_across_SNRs('<SIM_DIR>/myelin_sim_paired/results/paired_dipoles',2,1)

cmap = brewermap(6,'Paired');
colormap(cmap)

if ~isfolder([data_dir 'figures'])
    mkdir([data_dir 'figures'])
end

SNR_label = {'0', '-5', '-10', '-15', '-20'};
source_recon_alg = {'IID','EBB','MSP'}; % can be 'IID' (Minimum Norm Estimate), 'EBB' (Empirical Bayesian Beamformer), 'MSP' (Multiple Sparse Priors)

for a = source_recon_idx
    source_recon_alg_a = source_recon_alg{a};
    % collect free energy across SNRs s, repetitions j and scalings ii
    load([data_dir '/DLE_model_ev_' source_recon_alg_a '_sim_scale_' num2str(sim_i) '_1'],'F')
    F_1 = F(:,:,[1,4:9]);

    load([data_dir '/DLE_model_ev_' source_recon_alg_a '_sim_scale_' num2str(sim_i) '_2'],'F')
    F_2 = F(:,:,[1,4:9]);

    F = cat(2,F_1,F_2);

    SNR_Scale_F = squeeze(mean(F(:,:,2:7)- F(:,:,1),2));

    b3 = bar3([SNR_Scale_F(:,1), SNR_Scale_F(:,2), SNR_Scale_F(:,3), SNR_Scale_F(:,4),...
        SNR_Scale_F(:,5), SNR_Scale_F(:,6)]);
   
    set(gca,'XTickLabel',{'133','75','110','91','105','95'})
    set(gca,'YTickLabel',SNR_label(1:end))
    for k = 1:length(b3)
        b3(k).FaceAlpha = 0.6;
    end
    colormap(cmap)
    set(gcf,'color','w')
    xlabel('Scaling in %')
    ylabel('SNR in dB')
    zlabel('\Delta F')
    title(source_recon_alg_a,'FontSize',16,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')
    set(gcf, 'Units', 'Normalized', 'Position', [0.1300 0.22 0.24 0.36]);
    print('-djpeg', '-r300', sprintf('%s/figures/F_barplot_%s_across_SNRs', data_dir, source_recon_alg_a))
end