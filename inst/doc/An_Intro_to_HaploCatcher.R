## ---- setup, include = FALSE--------------------------------------------------
#set universal options
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- loading, echo=FALSE, include=FALSE--------------------------------------
#library
suppressWarnings(library(HaploCatcher))
suppressWarnings(library(ggplot2))

#read data from package
data("geno_mat")
data("gene_comp")
data("marker_info")

## ---- head_gene_comp, echo=FALSE----------------------------------------------
#display head of gene_comp file
knitr::kable(head(gene_comp), 
             align = "c", 
             caption = "Gene Compendium Input File")

## ---- display_calls, echo=FALSE-----------------------------------------------
#display unique calls
knitr::kable(data.frame(Call=unique(gene_comp$Call)), 
             align = "c", 
             caption = "Haplotype Calls")

## ---- data_example, echo=FALSE, echo=FALSE------------------------------------
#make example data frame
a<-data.frame(FullSampleName=c("Geno_1", "Geno_2", "Geno_3"),
              Chromosome=rep("5D", 3),
              Gene=rep("pm34", 3),
              Call=c("pm34", "het_pm34", "non_pm34"))

#display
knitr::kable(a, 
             align = "c", 
             caption = "Pm34 Example Input File")

## ---- head_geno_mat, echo=FALSE-----------------------------------------------
#must define as a matrix
geno_mat<-as.matrix(geno_mat)

#manipulate a little to demonstrate
geno_mat[1,1]=1
geno_mat[5,5]=2

#see head
knitr::kable(geno_mat[1:5, 1:5], 
             align = "c", 
             caption = "Marker Matrix Input File")

## ---- head_marker_info, echo=FALSE--------------------------------------------
#display marker_info
knitr::kable(head(marker_info), 
             align = "c", 
             caption = "Marker Info Input File")

## ---- marker_info_example, echo=FALSE-----------------------------------------
#make example data set
a<-data.frame(Marker=c("TRACE_00102",
                       "TRACE_13112",
                       "TRACE_43821"),
              Chromosome=c(1, 
                           1, 
                           1),
              BP_Position=c(0.0,
                            0.5,
                            1.2))
#display
knitr::kable(a, 
             align = "c", 
             caption = "Marker Info Input File Example")

## ---- figure_1, echo=FALSE, fig.align='center', out.width="100%"--------------
#pull figure from internet
url<-"https://raw.githubusercontent.com/zjwinn/HaploCatcher/main/Figure.png"
knitr::include_graphics(url)

## ---- example_read_data-------------------------------------------------------
#library
library(HaploCatcher)

#read data from package
data("geno_mat")
data("gene_comp")
data("marker_info")

## ---- example_make_training_and_test------------------------------------------
#set seed (for reproducible results)
set.seed(022294)

#randomly partition the training data and test from total
training_genotypes=sample(rownames(geno_mat), size = round(nrow(geno_mat)*0.8, 0))
testing_genotypes=rownames(geno_mat)[!rownames(geno_mat) %in% training_genotypes]

#nullify the seed we set so we don't mess with cross validation
set.seed(NULL)

## ---- example_models_1, fig.width=12, fig.height=8, fig.align='center', out.width="100%"----
#run without heterozygous individuals sequentially
results1<-auto_locus(geno_mat = geno_mat,
                     gene_file = gene_comp,
                     gene_name = "sst1_solid_stem",
                     marker_info = marker_info,
                     chromosome = "3B",
                     training_genotypes = training_genotypes,
                     testing_genotypes = testing_genotypes,
                     set_seed = 022294,
                     n_perms = 10,
                     verbose = FALSE)

## ---- display_results---------------------------------------------------------
#show results for with and without hets
knitr::kable(results1$predictions[1:15,], 
             align = "c", 
             caption = "Heterozygous Excluded Results")

## ---- visualize_cross_validation, fig.width=12, fig.height=8, fig.align='center', out.width="100%"----
#plot out the cross validation results
plot_locus_perm_cv(results1$cross_validation_results)

## ---- look_at_summary_results-------------------------------------------------
#look at overall and by_class performance
knitr::kable(results1$cross_validation_results$Overall_Summary[,1:5], 
             align = "c", 
             caption = "Subsection of Overall Summary for 'results1' Object")
knitr::kable(results1$cross_validation_results$By_Class_Summary[,1:5],
             align = "c", 
             caption = "Subsection of By-Class Summary for 'results1' Object")

