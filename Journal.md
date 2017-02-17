# CUVIER'S TEAM
# Heurteau
# SINGLE CELL PROJECT

# Configuration : https://github.com/u-n-i-v-e-r-z/single_Cell01 (Clone with ssh if on Genotoul)

The aim of this study is to analyze the effect of Mes4 depletion (siRNA) on Heterochromatin spreading and its effect upon genes expression.


# 02/02/17
Among the 19 single celle exp, all have been contaminated resulting in so few reads that will be used (max 40% of Uniquely mapped)
Starting to play with RNAseq script on Rstudio localy.




#03/02/17
We're doing the Alignement of Drosophila on dm6 genome with STAR (see STAR_alignment_noGFF.R)



#07/02/17
Script on Diff analysis and Variability analysis already done by Steph.


#14/02/17
The Diff Exp analysis has already been led by Steph, so we're focusing on analyzing the Diff variance.
To do so we have multiple option:
-Pascal's solution : NSD (Normalized Standard Dev).
-Raphael's solution : estimates the mean- variance relationship robustly and generates a precision weight for each individual with VOOM package.
-Steph's solution : EIQ (Ecart Inter Quartile), a way of analyzing dispersion on median based measure.
 
