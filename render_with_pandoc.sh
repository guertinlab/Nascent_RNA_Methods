wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/README.md
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/nascent_rna_methods.bib
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/springer-basic-parentheses-bold-italics.csl
pandoc --number-sections --citeproc  --bibliography nascent_rna_methods.bib -s README.md -o Methods_Ms.pdf --csl springer-basic-brackets.csl
