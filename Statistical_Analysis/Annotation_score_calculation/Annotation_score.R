library("tidyverse")

### We start by calculating the annotation score of the compounds of the RA-map V2 from which we will remove the compounds that are not found in the RASF-model.

## RA-map V2's annotation score calculation

MIRIAM_export = read.csv(file="Extraction_RA-Map_V2-MIRIAM.csv", sep=";")

MIRIAM_tmp = distinct(select(MIRIAM_export,-starts_with("MIRIAM.relation")))
MIRIAM = MIRIAM_tmp[!MIRIAM_tmp$MIRIAM.datatype.=="",]

tmp <- MIRIAM %>% group_by(class, name)
tmp1 <- tmp %>% dplyr::select(MIRIAM.datatype.,MIRIAM.id.)
tmp2 <- tmp %>% dplyr::select(MIRIAM.datatype..1,MIRIAM.id..1) %>% dplyr::rename(MIRIAM.datatype. = MIRIAM.datatype..1 , MIRIAM.id. = MIRIAM.id..1)

# Rearranging the MIRIAM annotations

my_func <- function(x){
  initial <- x %>% mutate_all(as.character) %>% group_by(class, name)
  final=NULL
  for(i in seq(3,544,2)){
    colnames_initial <-  colnames(initial)
    tmp <- dplyr::select(initial,colnames_initial[c(i,i+1)]) 
    colnames(tmp) <- c("class","name","datatype","pubmed_ids")
    final <- bind_rows(final,tmp)
  }
  return(final)
}


MIRIAM_tmp <- my_func(MIRIAM)

MIRIAM_tmp$pubmed_ids[MIRIAM_tmp$datatype != "PubMed"] = "0"
MIRIAM_tmp$datatype[MIRIAM_tmp$datatype != "PubMed"] = "PubMed"

MIRIAM_final <- subset(MIRIAM_tmp, datatype == "PubMed") %>% 
  dplyr::select(-datatype) %>% 
  summarize(pubmed=paste(pubmed_ids,collapse = ", ")) %>% 
  arrange(name,class)

# Count of unique pmid
MIRIAM_final$pubmed %>% strsplit(split = ",") %>% unlist %>% unique %>% length

# Removing the class information and merge all the PMIds for species with the same names + keep only unique PMIDs
MIRIAM_finalV2 <- MIRIAM_final[,-1] %>% 
  group_by(name) %>% 
  mutate(pubmed = sapply(pubmed, function(x) strsplit(x,", "))) %>% 
  group_by(name) %>% 
  summarise(pubmed= paste0(unique(unlist(pubmed)),collapse = ", "))

# Counting the number of PMIDs for each component
MIRIAM_finalV2$count <- MIRIAM_finalV2$pubmed %>% 
  str_split(pattern = ", " ) %>% 
  lapply(FUN=length) %>% 
  unlist

## Importing the RASF-model

species <- read.csv("RASF_model_components_1.csv", sep=";", header=F)


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


# Removing duplicates

unique_species= unique(species_name_only)

pheno_RA_map_V2=c("Angiogenesis","Bone Erosion","Apoptosis","Matrix Degradation","Osteoclastogenesis",
                  "Inflammation","Cell Chemotaxis","Recruitment","Infiltration",
                  "Cell growth","Survival","Proliferation", "Hypoxia")

unique_species=unique_species[!unique_species %in% pheno_RA_map_V2]



MIRIAM_finalV2 <- MIRIAM_finalV2[MIRIAM_finalV2$name %in% unique_species, ]


### Calculating annotation score's distribuion

zero<-sum(MIRIAM_finalV2$count==1)
one=sum(MIRIAM_finalV2$count==2)
two_five=sum(MIRIAM_finalV2$count>2 & MIRIAM_finalV2$count <= 5)
six_ten=sum(MIRIAM_finalV2$count>5 & MIRIAM_finalV2$count <= 10)
eleven_fifteen=sum(MIRIAM_finalV2$count>10 & MIRIAM_finalV2$count <= 15)
sixteen_twenty=sum(MIRIAM_finalV2$count>15 & MIRIAM_finalV2$count <= 20)
twenty=sum(MIRIAM_finalV2$count>20)



count <- c(zero, one, two_five, six_ten, eleven_fifteen, sixteen_twenty, twenty)
pie(count, clockwise = TRUE, labels = c("0", "1", "2-5", "6-10", "11-15", "16-20", ">20"))

