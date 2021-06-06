library(Seurat)
library(tidyverse)
library(readxl)
library(viridis)
library(ggrepel)
library(vegan)
library(ggpubr)

#create plotting subfolders
dir.create(file.path("plots"))
dir.create(file.path("plots", "umap"))
dir.create(file.path("plots", "heatmaps"))
dir.create(file.path("plots", "others"))

#load functions and colors
source(file.path("R", "functions.R"))

#load data
load(file.path("data", "seurat.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= all@meta.data[,"seurat_clusters"], row.names = rownames(all@meta.data)) %>%
  bind_cols(as.data.frame(t(all[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(all) <- order_clusters

#extract metadata
metadata <- all@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = order_clusters)
metadata$cell_type <- case_when(
  metadata$seurat_clusters %in% c("3")   ~ "CAMs",
  metadata$seurat_clusters %in% c("0", "1", "4", "7", "10") ~ "MG",
  metadata$seurat_clusters %in% c("9") ~ "Granulocytes",
  metadata$seurat_clusters %in% c("5", "8", "11") ~ "Monocytes",
  metadata$seurat_clusters %in% c("2") ~ "Lymphocytes",
  TRUE ~ "other"
)

#conditions
metadata$condition <- case_when(
  grepl("(Onset|EAE)", rownames(metadata)) ~ "EAE",
  grepl("(Acute)", rownames(metadata)) ~ "Acute",
  grepl("(Precl_)", rownames(metadata)) ~ "Precl",
  TRUE ~ "Ctrl"
)

#assign variables to seurat object
all$celltype <- metadata$cell_type
all$condition <- metadata$condition

# cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 3, object = all, .retain_cl = levels(all)) + labs(subtitle= paste0(i,' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap", paste0(i,"_signature.pdf")), useDingbats=F)
}

# clusters umap
all %>%
  DimPlot(label=T, pt.size = 3) +
  theme_void() +
  scale_color_manual(values=colors_many)
ggsave(file.path("plots", "umap", "clusters_umap.pdf"), useDingbats=F)

# cell types umap
all %>%
  DimPlot(label=T, pt.size = 3, group.by = "celltype") +
  theme_void() +
  scale_color_manual(values=colors_pat)
ggsave(file.path("plots", "umap", "celltypess_umap.pdf"), useDingbats=F)

# conditions umap
all %>%
  DimPlot(label=T, pt.size = 3, group.by = "condition") +
  theme_void() +
  scale_color_manual(values=rev(colors_pat))
ggsave(file.path("plots", "umap", "conditions_umap.pdf"), useDingbats=F)

#gene heatmap
if (!file.exists(file.path("data", "markers.RData"))) {
  markers <- FindAllMarkers(all)          
  save(markers, file = file.path("data", "markers.RData"))
  write_csv(markers, file.path("data", "markers.csv"))
} else {
  load(file.path("data", "markers.RData"))
}

top10 <- markers %>% 
  #remove non annotated genes
  filter(!gene %in% grep("(^Gm|^Rp|Rik|mt-|RP)", .$gene, value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(all,features = top10$gene, group.colors = colors_pat) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(file.path("plots", "heatmaps", "heatmap.pdf"), useDingbats=F)

DoHeatmap(all,features = top10$gene, group.bar = F) +
  theme_void() + 
  scale_fill_viridis(option = "B") +
  NoLegend()

ggsave(file.path("plots", "heatmaps", "heatmap_color_panel.png"))

#volcano plot
all2 <- all

Idents(all2) <- paste(metadata$condition, metadata$celltype, sep = "_")

#epiplexus
if (!file.exists(file.path("data", "diffgenes_ctrl_vs_aea.RData"))) {
  all2_genes <- FindMarkers(all2, 
                           ident.1 = "EAE_MG",
                           ident.2 = "Ctrl_MG",
                           logfc.threshold = 0,
                          min.pct = 0) %>%
    rownames_to_column(var="gene") %>%
    mutate(Condition_upregulated = case_when(
      avg_log2FC > 0 ~ "EAE",
      TRUE ~ "Ctrl"))
  
  save(all2_genes, file = file.path("data", "diffgenes_ctrl_vs_aea.RData"))
  write.csv(all2_genes, file = file.path("data", "diffgenes_ctrl_vs_aea.csv"))
} else {
  load(file.path("data", "diffgenes_ctrl_vs_aea.RData"))
}

#remove unannotated genes
all2_genes <- all2_genes %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP|ERCC)", .$gene))

all2_genes <- all2_genes %>%
  mutate(
    genes_sig = ifelse(p_val_adj < .05 & abs(avg_log2FC) > .2, "sig.", "not sig."),
    show_genes = ifelse(genes_sig == "sig.", gene, NA),
    
  ) 


volcano <- ggplot(all2_genes, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  annotate("rect", xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf,
            fill=colors_many[2], alpha =.4)+
  annotate("rect" ,xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf,
            fill=colors_many[3], alpha =.4) +
  geom_point(size=5) + 
  geom_text_repel(size=7, box.padding=1.15, max.overlaps = 30) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. log2FC", y="-log10 transf. adj. p-value") 

volcano 

ggsave(file.path("plots", "others", "eae_vs_ctrl_volcano.pdf"), useDingbats=F)


#ene expression umaps
genes <- c("Mrc1", "Stab1", "P2ry12", "Apoe", "Nr4a1", "Ccr2", "Enpp2", "H2-Aa", "Cd209a", "Hexb", "Fos", "Nkg7", "Tmem119", "Olfml3","Sall1")

for (i in genes) {
  plt <- plot_expmap_seurat(features=i, object=all,  point_size = 3,.retain_cl = levels(all))
  print(plt)
  ggsave(file.path("plots", "umap", paste0(i,".pdf")), useDingbats=F)
}  
