species <- read.csv("RASF_model_components_2.csv", sep=";", header=F)

df = data.frame(species$V2)

# Componenents' name preprocessing

df$species.V2 = gsub("_simple_molecule", "",
                      gsub("_space_", "",
                           gsub("_Mitochondrion", "",
                                gsub("_Nucleus", "",
                                     gsub("_Secreted_space_Molecules","",
                                          gsub("_active", "",
                                               gsub("_Cytoplasm","",
                                                    gsub("_complex","",
                                                         gsub("_Extracellular_space_Space", "",
                                                              gsub("_ion","",
                                                                   gsub("_slash_","/",
                                                                        gsub("_empty","",
                                                                             gsub("_phenotype", "",
                                                                                  gsub("_sub_", "",
                                                                                       gsub("_endsub_", "",
                                                                                            gsub("_minus_","-",
                                                                                                 gsub ("_underscore_","_",
                                                                                                       gsub("_plus_", "+",
                                                                                                            gsub("_rna", "",
                                                                                                                gsub("_phosphorylated", "", df$species.V2))))))))))))))))))))

# separate the molecular complexes into unique components 

species <-tidyr::separate(
  data = df,
  col = species.V2,
  sep = "/",
  into = c("A","B","C","D","E","F","G"),
  remove = T
)

species_name_only=c(species$A,species$B,
                    species$C,species$D,
                    species$E,species$F,
                    species$G)

species_name_only=species_name_only[!is.na(species_name_only)]


# Remove duplicates

unique_species= unique(species_name_only)

# Remove the phenotypes (not considered in the cell-specificity calculations)


pheno_RA_map_V2=c("Angiogenesis","Bone Erosion","Apoptosis","Matrix Degradation","Osteoclastogenesis",
                  "Inflammation","Cell Chemotaxis","Recruitment","Infiltration",
                  "Cell growth","Survival","Proliferation", "Hypoxia")

unique_species=unique_species[!unique_species %in% pheno_RA_map_V2]


# Import the cell=speciic lists files

RASF_overlay = read.csv("RASF.csv", sep=";", header=TRUE) 
macro_overlay = read.csv("macrophages.csv", sep=";", header=TRUE)
PBMC_overlay = read.csv("PBMC.csv", sep=";", header=TRUE)
chondrocytes_overlay = read.csv("chondrocytes.csv", sep=";", header=TRUE) 
serum_overlay = read.csv("serum.csv", sep=";", header=TRUE)
synovial_fluid_overlay = read.csv("synovial_fluid.csv", sep=";", header=TRUE)
synovial_tissue_overlay = read.csv("synovial_tissue.csv", sep=";", header=TRUE)
blood_overlay = read.csv("blood_components.csv", sep=";", header=TRUE)


# Calculate the percentages of cell=specificity

RASF = round(((sum(tolower(unique_species) %in% tolower(RASF_overlay$name), na.rm = TRUE))/length(unique_species))*100)
Synovial_tissue = round(((sum(tolower(unique_species) %in% tolower(synovial_tissue_overlay$name), na.rm = TRUE))/length(unique_species))*100)
Synovial_fluid = round(((sum(tolower(unique_species) %in% tolower(synovial_fluid_overlay$name), na.rm = TRUE))/length(unique_species))*100)
PBMC = round(((sum(tolower(unique_species) %in% tolower(PBMC_overlay$name), na.rm = TRUE))/length(unique_species))*100)
Blood = round(((sum(tolower(unique_species) %in% tolower(blood_overlay$name), na.rm = TRUE))/length(unique_species))*100)
Chondrocytes = round(((sum(tolower(unique_species) %in% tolower(chondrocytes_overlay$name), na.rm = TRUE))/length(unique_species))*100)
Macrophages = round(((sum(tolower(unique_species) %in% tolower(macro_overlay$name), na.rm = TRUE))/length(unique_species))*100)
Serum = round(((sum(tolower(unique_species) %in% tolower(serum_overlay$name), na.rm = TRUE))/length(unique_species))*100)

# Plot the RASF=model's cell=specificity

perc=data.frame(c(RASF, Synovial_tissue, Synovial_fluid, Macrophages, Chondrocytes, 
                  PBMC, Blood, Serum))
row.names(perc)=c("RASF","Synovial Tissue", "Synovial Fluid", "Macrophages", "Chondrocytes", "PBMC", "Blood", "Serum")

library(data.table)
setDT(perc, keep.rownames=TRUE)[]
colnames(perc)=c("Cell_Type", "Percentage")

perc$Cell_Type = factor(perc$Cell_Type, levels = perc$Cell_Type[order(-perc$Percentage)])

p = ggplot(data=perc, aes(x=Cell_Type, y=Percentage)) + 
  geom_bar(stat="identity", fill="darkblue") +
  geom_text(aes(label=Percentage), vjust=-0.7, size=3.2)+
  theme_minimal()
p + ylim(0,100) + ggtitle("RASF-model cell-type specificity")


