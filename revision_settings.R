mytheme <-     theme(
  # text = element_text(family=fontFamily),
  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.background = element_rect(fill = "transparent"),
  panel.grid.major.x =  element_blank(),
  panel.grid.minor.x =  element_blank(),
  axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
  legend.title = element_text(face="bold")
) 

all_pairs <- c(
  file.path("LI_40kb","GSE105381_HepG2_40kb", "TCGAlihc_norm_lihc"),
  file.path("LG1_40kb" ,"ENCSR444WCZ_A549_40kb", "TCGAluad_norm_luad"),
  file.path("LG2_40kb" ,"ENCSR444WCZ_A549_40kb" ,"TCGAluad_norm_luad"),
  file.path("LG1_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "TCGAprad_norm_prad"),
  file.path("GSE118514_RWPE1_40kb", "GSE118514_22Rv1_40kb", "TCGAprad_norm_prad")
)

all_normal_ds <- as.character(sapply(all_pairs, function(x) dirname(dirname(x))))
all_tumor_ds <-  as.character(sapply(all_pairs, function(x) basename(dirname(x))))


all_pairs_dt <- data.frame(
  hicds1 = dirname(dirname(all_pairs)),
  hicds2 = basename(dirname(all_pairs)),
  exprds = basename(all_pairs),
  stringsAsFactors = FALSE
)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


all_hicds <- c(
 "Barutcu_MCF-10A_40kb"="AWS_Barutcu_MCF-10A",
  "Barutcu_MCF-7_40kb"="AWS_Barutcu_MCF-7",
  "ENCSR079VIJ_G401_40kb" ="mega_ENCSR079VIJ_G401",
  "ENCSR312KHQ_SK-MEL-5_40kb"="mega_ENCSR312KHQ_SK-MEL-5",
  "ENCSR401TBQ_Caki2_40kb"="mega_ENCSR401TBQ_Caki2",
  "ENCSR504OTV_transverse_colon_40kb"="ENCSR504OTV_transverse_colon",
  "ENCSR549MGQ_T47D_40kb"="mega_ENCSR549MGQ_T47D",
  "ENCSR862OGI_RPMI-7951_40kb"="mega_ENCSR862OGI_RPMI-7951",
  "GSE105194_cerebellum_40kb"="mega_GSE105194_cerebellum",
  "GSE105194_spinal_cord_40kb"="mega_GSE105194_spinal_cord",
  "GSE105318_DLD1_40kb"="mega_GSE105318_DLD1", 
  "GSE109229_BT474_40kb"="GSE109229_BT474", 
  "GSE109229_SKBR3_40kb"="GSE109229_SKBR3",
  "GSE118588_Panc_beta_40kb"="GSE118588_Panc_beta",
  "GSE99051_786_O_40kb" = "GSE99051_786_O",
  "PA2_40kb"="Compendium_PA2",
  "PA3_40kb"="Compendium_PA3",
  "Panc1_rep12_40kb"="mega_Panc1_rep12",
  "Rao_HCT-116_2017_40kb"="AWS_Rao_HCT-116_2017",
  "K562_40kb"="AWS_K562",
  "HMEC_40kb"="AWS_HMEC",
    "LI_40kb"="Compendium_LI",
  "GSE105381_HepG2_40kb"="mega_GSE105381_HepG2",
  "LG1_40kb" ="Compendium_LG1",
  "ENCSR444WCZ_A549_40kb"="mega_ENCSR444WCZ_A549",
  "LG2_40kb"="Compendium_LG2",
  "ENCSR489OCU_NCI-H460_40kb"="mega_ENCSR489OCU_NCI-H460",
  "GSE118514_RWPE1_40kb"="mega_GSE118514_RWPE1",
  "ENCSR346DCU_LNCaP_40kb"="mega_ENCSR346DCU_LNCaP",
  "GSE118514_22Rv1_40kb"="GSE118514_22Rv1"
)


all_tissues <- c("Barutcu_MCF-10A_40kb"="breast",
              "Barutcu_MCF-7_40kb"="breast",               
"ENCSR079VIJ_G401_40kb"="kidney",
 "ENCSR312KHQ_SK-MEL-5_40kb"="skin",        
 "ENCSR346DCU_LNCaP_40kb"="prostate",            "ENCSR401TBQ_Caki2_40kb"="kidney",           
 "ENCSR444WCZ_A549_40kb"="lung",             "ENCSR489OCU_NCI-H460_40kb"="lung",        
 "ENCSR504OTV_transverse_colon_40kb"="colon", "ENCSR549MGQ_T47D_40kb"="breast",            
 "ENCSR862OGI_RPMI-7951_40kb"="skin",        "GSE105194_cerebellum_40kb"="brain",        
 "GSE105194_spinal_cord_40kb"="brain",        "GSE105318_DLD1_40kb"="colon",              
 "GSE105381_HepG2_40kb"="liver",              "GSE109229_BT474_40kb"="breast",             
 "GSE109229_SKBR3_40kb"="breast",              "GSE118514_22Rv1_40kb"="prostate",             
"GSE118514_RWPE1_40kb"="prostate",              "GSE118588_Panc_beta_40kb"="pancreas",         
 "GSE99051_786_O_40kb"="kidney",               "HMEC_40kb"="breast",                        
 "K562_40kb"="lymph",                         "LG1_40kb"="lung",                         
 "LG2_40kb"="lung",                          "LI_40kb"="liver",                          
 "PA2_40kb"="pancreas",                          "PA3_40kb"="pancreas",                         
 "Panc1_rep12_40kb"="pancreas",                  "Rao_HCT-116_2017_40kb"="colon"           
)

