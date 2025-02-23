
#nichenet. Use Rama's DEGS
library("nichenetr"); library("tibble"); library("Seurat"); library("tidyr"); library("dplyr"); library("ggplot2"); library("purrr")
options(timeout=3600)
# lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
# ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
# weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
# save(list=c("lr_network", "ligand_target_matrix", "weighted_networks"), file="nichenet_data.RData")
load("rama_scs.RData")
load("nichenet_data.RData")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

rama_scs@meta.data %>% mutate(nichenet_labels = case_when(
    merged_cluster_scs ==  "DonorBCells-DLN" ~ "rama",
    merged_cluster_scs ==  "DonorBCells-NDLN" ~ "rama",
    TRUE ~ merged_cluster_scs)) -> rama_scs@meta.data

Idents(rama_scs) <- "merged_cluster_scs"
## 1 receiver
receiver = "rama"
DefaultAssay(rama_scs) <- "alra"
expressed_genes_receiver = get_expressed_genes("DonorBCells-NDLN", rama_scs, pct = 0.10, assay_oi = "alra")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## 2 sender
sender_celltypes = c("5_SSM")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, rama_scs, 0.10,  assay_oi = "SCT") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
# 3
Idents(rama_scs) <- "nichenet_labels"
seurat_obj_receiver = subset(rama_scs, idents = receiver)
Idents(seurat_obj_receiver) <- "merged_cluster_scs"
condition_reference = "DonorBCells-DLN"
condition_oi = "DonorBCells-NDLN" 

DE_table_receiver = readxl::read_excel("DEG_Bmem_HeatMap.xlsx") %>%
    filter(group==3 | group==4 | group == 7 | group == 2) 
#DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% tibble::rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% pull(genes) 
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
# 4
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
# 5
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
DefaultAssay(rama_scs) <- "SCT"
Idents(rama_scs) <- "merged_cluster_scs"
#DotPlot(rama_scs, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network
ggsave("lig-target-NDLN-up-n=1154.svg", p_ligand_target_network, height=15, width=30)
# 6 receptors
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to, -database, -source) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
order_ligands_receptor = order_ligands %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
ggsave("lig-receptor-NDLN-up-n=1154.svg", p_ligand_receptor_network)

#circos
Idents(rama_scs) <- "merged_cluster_scs"
ligand_expression_tbl = tibble(
    ligand = best_upstream_ligands, 
    SSM = FetchData(
        subset(rama_scs, idents = c("5_SSM")),
        vars = best_upstream_ligands,
        slot = "data") %>%
        apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)})
)

SSM_specific_ligands = ligand_expression_tbl %>% pull(ligand)

ligand_type_indication_df = tibble(
    ligand_type = c(rep("SSM_ligands", times = SSM_specific_ligands %>% length())),
    ligand = c(SSM_specific_ligands))
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "bcell") %>% inner_join(ligand_type_indication_df)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
grid_col_ligand =c("SSM_ligands" = "firebrick3")
grid_col_target =c(
    "bcell" = "darkgreen")
grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)
circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% dplyr::select(ligand,target, weight)
ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

target_order = circos_links$target %>% unique()
ligand_order = c(SSM_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5
gaps = c(
    rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "SSM_ligands") %>% distinct(ligand) %>% nrow() -1)), 
    width_ligand_target,
    rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "bcell") %>% distinct(target) %>% nrow() -1)),
    width_ligand_target
)
library("circlize")
svg("ligand_target_circos-NDLN-up-n=1154.svg", width = 15, height = 15)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()
dev.off()

# ligand activity heatmap
ligand_aupr_matrix = ligand_activities %>%
    select(aupr_corrected) %>% 
    as.matrix() %>%
    magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity",
                                                        color = "darkorange",
                                                        legend_position = "top",
                                                        x_axis_position = "top",
                                                        legend_title = "AUPR\n(target gene prediction ability)") +
    theme(legend.text = element_text(size = 9))
ggsave("lig-activity-NDLN-up-n=1154.svg", p_ligand_aupr, width=3, height=10)


save.image("nichenet_ndln_n=1154.RData")
