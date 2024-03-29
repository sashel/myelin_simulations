function [] = plot_model_ev_paired(data_dir,source_recon_idx, sim_i)
% use as plot_model_ev_paired('<SIM_DIR>/myelin_sim_paired/results/paired_dipoles',2,1)

cmap = brewermap(6,'Paired');
colormap(cmap)

if ~isfolder([data_dir '/figures'])
    mkdir([data_dir '/figures'])
end

scaling_label = {'100-100%', '133-75%', '75-133%', '110-91%','91-110%','105-95%','95-105%'};
source_recon_alg = {'IID','EBB','MSP'}; % can be 'IID' (Minimum Norm Estimate), 'EBB' (Empirical Bayesian Beamformer) or'MSP' (Multiple Sparse Priors)

for a = source_recon_idx
    source_recon_alg_a = source_recon_alg{a};
    % collect free energy across SNRs s, repetitions j and scalings ii
    load([data_dir '/DLE_model_ev_' source_recon_alg_a '_sim_scale_' num2str(sim_i) '_1'],'F')
    F_1 = F(:,:,[1,4:9]);

    load([data_dir '/DLE_model_ev_' source_recon_alg_a '_sim_scale_' num2str(sim_i) '_2'],'F')
    F_2 = F(:,:,[1,4:9]);

    F = cat(2,F_1,F_2);
    SNR_Scale_F = squeeze(mean(F(1,:,2:7)- F(1,:,1),2));

    cm = colormap(brewermap(6, 'Paired'));

    h = bar(SNR_Scale_F,'linewidth',1.6);
    h.BaseLine.LineWidth = 1.6;
    hold on
    for i = 1:length(SNR_Scale_F)
        h = bar(i,SNR_Scale_F(i));
        if rem(i,2)
            set(h,'FaceColor',cm(i+1,:));
        else
            set(h,'FaceColor',cm(i-1,:));
        end
    end

    xlim([0 7])
    ylim([-4 6])

    box('off')
    set(gca,'XTickLabel',{'133-75', '75-133', '110-91','91-110','105-95','95-105'})
    colormap(cmap)
    set(gcf,'color','w')
    xlabel('Scaling in %')
    ylabel('\Delta F')
    set(gca,'linewidth',1.6)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',18)
    title(sprintf('Scaling: %s',scaling_label{sim_i}),'FontSize',24,'FontName','Helvetica','FontWeight', 'bold','FontAngle','italic','Interpreter','none')

    saveas(gcf,sprintf('%s/figures/infer_myelination_%s_scaling_%d', data_dir,source_recon_alg_a,sim_i),'svg')
    saveas(gcf,sprintf('%s/figures/infer_myelination_%s_scaling_%d', data_dir,source_recon_alg_a,sim_i),'jpg')
end