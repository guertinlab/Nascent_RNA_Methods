wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/README.md
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/nascent_rna_methods.bib
#pandoc --citeproc --number-sections --bibliography=nascent_rna_methods.bib -s README.md -o Methods_Ms.pdf
pandoc --number-sections --bibliography=nascent_rna_methods.bib -s README.md -o Methods_Ms.pdf
