## TP53 protein structure and function analysis
## PDB Structures Used: 20CJ,2ACO,2LZH,3DO5

## install packages
install.packages("bio3d")
install.packages("tibble")
install.packages("scales")

library(bio3d)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)

# 1.downloaded TP53 structures from RCSB PDB
wt_apo   = read.pdb("2OCJ")   # wild
wt_dna   = read.pdb("2AC0")   # wild + DNA
mut_r175 = read.pdb("2LZH")   # mutant
mut_r273 = read.pdb("3D05")   # mutant

print(wt_apo)
print(wt_dna)
print(mut_r175)
print(mut_r273)

# 2.amino acid composition
# i will look at chain A and CA atoms for all structures
chain_wt= wt_apo$atom[wt_apo$atom$chain == "A" & wt_apo$atom$elety == "CA", ]

# counting each amino acid type
aa_counts_wt = table(chain_wt$resid)
aa_counts_wt
aa_df_wt= as.data.frame(aa_counts_wt)
aa_df_wt
colnames(aa_df_wt)= c("residue", "freq")
aa_df_wt$percentage= aa_df_wt$freq / sum(aa_df_wt$freq) * 100
aa_df_wt$structure= "WT apo (2OCJ)"

# for R175H mutant
chain_mut= mut_r175$atom[mut_r175$atom$chain == "A" & mut_r175$atom$elety == "CA",]
aa_counts_mut = table(chain_mut$resid)
aa_df_mut= as.data.frame(aa_counts_mut)
aa_df_mut
colnames(aa_df_mut)= c("residue", "freq")
aa_df_mut$percentage =aa_df_mut$freq /sum(aa_df_mut$freq) * 100
aa_df_mut$structure="R175H (2LZH)"

## combining both to make plots
aa_combine =rbind(aa_df_wt, aa_df_mut)

aa_combine$category[aa_combine$residue %in% c("SER","THR","CYS","TYR","ASN","GLN")] ="polar"  
aa_combine$category[aa_combine$residue %in% c("ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO")] ="hydrophobic"
aa_combine$category[aa_combine$residue %in% c("ARG","LYS","HIS")] ="positively charged"
aa_combine$category[aa_combine$residue %in% c("ASP","GLU")] ="negatively charged"
aa_combine$category[aa_combine$residue == "GLY"] ="glycine"

# plots
install.packages("cowplot")
library(viridis)
library(cowplot)
library(ggplot2)
colors()
my_colors = c(
"hydrophobic" = "gold",
"positively charged" = "lightblue",
"negatively charged" = "coral",
"polar" = "lightgreen",
"glycine" = "lavenderblush3")

plot1 = ggplot(aa_combine, aes(x = reorder(residue, -percentage), y = percentage, fill = category)) +
geom_bar(stat = "identity") +
facet_wrap(~ structure) +
scale_fill_manual(values = my_colors) +
labs(title = "Amino acid composition of TP53 DBD",
x = "Amino acid", y = "% of residues",
fill = "Property") +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
strip.background = element_rect(fill = "gray"))

plot1

## B-factor analysis (flexibility analysis)
# l2 and l3 loops are flexible
get_bfact = function(pdb_obj, chain_id = "A") {
  ca_atoms = pdb_obj$atom[pdb_obj$atom$chain == chain_id & pdb_obj$atom$elety == "CA", ]
  return(ca_atoms[, c("resno", "resid", "b")]) 
}
bf_wt = get_bfact(wt_apo)
bf_wt
bf_dna = get_bfact(wt_dna)
bf_dna
bf_r273 = get_bfact(mut_r273)
bf_r273
bf_wt$structure   = "WTapo"
bf_dna$structure  = "WT+DNA"
bf_r273$structure = "R273H mutant"

## now i will merge them
bf_all = rbind(bf_wt, bf_dna, bf_r273)
bf_all$structure = factor(bf_all$structure, levels = c("WTapo", "WT+DNA", "R273H mutant"))

# L2 and L3 loop regions, approximate boundaries
# loop naming and general location from: Cho, Gorina, Jeffrey & Pavletich (1994) Science 265:346-355
# exact residue ranges approximated from PDB structure inspection
# loop l2 = 164-194(zinc binding) and l3 = 237-250(DNA contact)
plot2 = ggplot(bf_all, aes(x = resno, y = b, color = structure)) +
annotate("rect", xmin = 162, xmax = 195, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.15) +
annotate("rect", xmin = 235, xmax = 250, ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.15) +
annotate("text", x = 179, y = 80, label = "L2 loop", size = 3, color = "darkorange4") +
annotate("text", x = 243, y = 80, label = "L3 loop", size = 3, color = "steelblue4") +
geom_smooth(se = FALSE, span = 0.15) +   
labs(title = "Flexibility across TP53 DNA-binding domain", x = "Residue number", y = "B-factor", color = "Structure") +
theme_bw()

plot2

# to see where the flexibility changes between WT and R273H so i will do a simple difference plot
shared_res = intersect(bf_wt$resno, bf_r273$resno)
length(shared_res)
bf_wt_shared = bf_wt[bf_wt$resno %in% shared_res, ]
bf_mut_shared = bf_r273[bf_r273$resno %in% shared_res, ] 

bf_wt_shared = bf_wt_shared[order(bf_wt_shared$resno), ]
bf_mut_shared = bf_mut_shared[order(bf_mut_shared$resno), ]

all(bf_wt_shared$resno == bf_mut_shared$resno)

delta_b = data.frame(
  resno = bf_wt_shared$resno,
  delta = bf_mut_shared$b - bf_wt_shared$b
)
delta_b$direction = ifelse(delta_b$delta > 0, "More flexible in mutant", "More rigid in mutant")

p_delta = ggplot(delta_b, aes(x = resno, y = delta, fill = direction)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  labs(title = "B-factor difference: R273H mutant minus WT",
       x = "Residue", y = "Delta B-factor",
       fill = "") +
  theme_bw()

p_delta

# step 4 keeps breaking turns out fit.xyz() needs the coordinates in some different format than what i have
# found this Kabsch algorithm online
kabsch_align = function(P, Q) {
cP = colMeans(P)
cQ = colMeans(Q)
P0 = sweep(P, 2, cP)
Q0 = sweep(Q, 2, cQ)
H = t(Q0) %*% P0
s = svd(H)
d= sign(det(s$v %*% t(s$u)))
D = diag(c(1, 1, d))
R = s$v %*% D %*% t(s$u)
Q_aligned = t(R %*% t(Q0))
Q_aligned = sweep(Q_aligned, 2, cP, "+")
return(Q_aligned)
}

# STEP 4: structural superposition and per-residue RMSD
get_ca = function(pdb_obj, chain_id = "A") {
atoms = pdb_obj$atom
ca = atoms[atoms$chain == chain_id & atoms$elety == "CA", ]
xyz = as.matrix(ca[, c("x","y","z")])
rownames(xyz) = ca$resno
return(list(resno = ca$resno, xyz = xyz))
}
calc_rmsd_perres = function(ref_pdb, mob_pdb, chain = "A") {
ref = get_ca(ref_pdb, chain)
mob = get_ca(mob_pdb, chain)
# only keeping residues that show up in both
common = intersect(ref$resno, mob$resno)
  ri = ref$xyz[as.character(common), , drop = F]
  mi = mob$xyz[as.character(common), , drop = F]
# was using fit.xyz() here before 
  mob_fit = kabsch_align(ri, mi)
  per_res = sqrt(rowSums((ri - mob_fit)^2))
 return(data.frame(resno = common, rmsd = per_res))
}

# WT vs WT+DNA - does DNA binding change the shape much
rmsd_wt_dna = calc_rmsd_perres(wt_apo, wt_dna)
rmsd_wt_dna$comparison = "WT vs WT+DNA"

# WT vs R273H mutant
rmsd_wt_r273 = calc_rmsd_perres(wt_apo, mut_r273)
rmsd_wt_r273$comparison = "WT vs R273H"

rmsd_all = rbind(rmsd_wt_dna, rmsd_wt_r273)

plot3 = ggplot(rmsd_all, aes(x = resno, y = rmsd, color = comparison)) +
geom_smooth(method = "loess", span = 0.12, se = F, linewidth = 1.5) +
geom_vline(xintercept = 273, linetype = "dashed", color = "gray") +
annotate("text", x = 276, y = 4.5, label = "R273", size = 3) +
scale_color_manual(values = c("WT vs WT+DNA" = "forestgreen","WT vs R273H"  = "orange1")) +
labs(title = "Per-residue RMSD relative to WT p53 (2OCJ)",
x = "Residue number", y = "RMSD(A)", color = "") +
theme_bw() 

plot3



# STEP 5: pairwise RMSD matrix (heatmap)
# comparing the 3 crystal ones against each other again not using 2LZH bc its NMR

structs = list(
  "WT apo"     = wt_apo,
  "WT + DNA"   = wt_dna,
  "R273H mut"  = mut_r273
)
structs
n = length(structs)
n
struct_names = names(structs)

rmsd_matrix = matrix(NA, nrow = n, ncol = n)
rmsd_matrix
rownames(rmsd_matrix) = struct_names
colnames(rmsd_matrix) = struct_names
rmsd_matrix
for (i in 1:n) {
for (j in 1:n) {
if (i == j) {
rmsd_matrix[i, j] = 0
} else {
ref = get_ca(structs[[i]])
mob = get_ca(structs[[j]])
# this part was wrong before, i was just cutting both to the same length instead of actually matching residue numbers,
 common = intersect(ref$resno, mob$resno)
      
ri = ref$xyz[as.character(common), , drop = F]
mi = mob$xyz[as.character(common), , drop = F]
fitted = kabsch_align(ri, mi)
rmsd_matrix[i, j] = round(sqrt(mean(rowSums((ri - fitted)^2))), 2)
}
}
}
rmsd_matrix
rmsd_df = as.data.frame(rmsd_matrix)
rmsd_df$ref_struct = rownames(rmsd_df)

rmsd_long = pivot_longer(rmsd_df, cols = -ref_struct,
                         names_to = "comp_struct", values_to = "rmsd")

plot4 = ggplot(rmsd_long, aes(x = ref_struct, y = comp_struct, fill = rmsd)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = ifelse(rmsd == 0, "-", paste0(rmsd, " A"))), size = 4) +
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "RMSD (A)") +
  labs(title = "Pairwise RMSD between TP53 structures", x = NULL, y = NULL) +
  theme_minimal()
plot4

##STEP 6: mutation frequency data
mutations = c("R175H", "G245S", "R248W", "R248Q", "R249S", "R273H", "R273C", "R282W")
freq = c(6.3, 2.8, 4.1, 3.2, 2.1, 5.9, 3.7, 2.0)

# structural vs DNA-contact
mut_class = c("Structural", "DNA-contact", "DNA-contact", "DNA-contact", "DNA-contact", "DNA-contact", "DNA-contact", "Structural")
mut_class
hotspot_df = data.frame(mutation = mutations, frequency = freq, class = mut_class)
hotspot_df

plot5 = ggplot(hotspot_df, aes(x = reorder(mutation, -frequency),y = frequency, fill = class)) +
geom_bar(stat = "identity", width = 0.6) +
geom_text(aes(label = paste0(frequency, "%")), vjust = -0.5, size = 3.2) +
scale_fill_manual(values = c("Structural" = "maroon", "DNA-contact" = "navy")) +
labs(title = "TP53 hotspot mutations in human cancer",x = "Mutation", y = "Frequency (%)", fill = "Type") +
theme_classic()
plot5
## SAVING ALL PLOTS
ggsave("TP53_PLOT1.png", plot1, width = 10, height = 6, dpi = 300)
ggsave("TP53_PLOT2.png", plot2, width = 10, height = 6, dpi = 300)
ggsave("TP53_p_delta.png", p_delta, width = 10, height = 6, dpi = 300)
ggsave("TP53_PLOT3.png", plot3, width = 10, height = 6, dpi = 300)
ggsave("TP53_plot4.png", plot4, width = 10, height = 6, dpi = 300)
ggsave("TP53_plot5.png", plot5, width = 10, height = 5, dpi = 300)

