wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/README.md

pandoc --filter pandoc-citeproc --bibliography=nascent_rna_methods.bib -s README.md -o Methods_Ms.pdf
