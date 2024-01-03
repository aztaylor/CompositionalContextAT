# Recreation and modification of models used in the paper "Biophysical Constraints Arising from the Compositional Context of Synthetic Gene Networks".
## Alec Taylor, Bahareh Dehghan Banadaki, Enoch Yeung.

This code arose from a final project which attempted to recreate and modify the models of, 
described in the titular paper from Yeung et al. which describes dynamics of cell free 
expression of convergent genes in a synthetic network.
Our aim was to modify the code to include the effect of DNA binding proteins (specifically
targeted CRISPR-Cas proteins and protein the potential protein fusions). 
This is still a work in progress.

Currently, we are working towards a Julia Implementation of the model that is modularized
such that it can add more genes, different context, and targeted binding proteins with 
potential modifications (such as the fusion of topoisomerases.)

The goal is to provide a Julia Package with the ability to incorporate these modifications
for use in circuit design, with minimal effort from the end user.

All work is inspired by Enoch Yeung (director of the biological control lab (BCL) at
UCSB), current and former members of BCL, and the many people who inspired this work 
beforehand.

Please refer to the latex document to explicitly see changes to syntax and model design.
