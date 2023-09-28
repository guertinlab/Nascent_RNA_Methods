wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/README.md
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/nascent_rna_methods.bib
#adapted from https://raw.githubusercontent.com/citation-style-language/styles/master/springer-basic-brackets.csl
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/springer-basic-parentheses-bold-italics.csl
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/Figure_1.pdf
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/Figure_2.pdf
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/Figure_3.pdf
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/Figure_4.pdf
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/Figure_5.pdf
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/Figure_6.pdf

pandoc --number-sections --citeproc --bibliography nascent_rna_methods.bib -s README.md -o Methods_Ms.pdf --csl springer-basic-parentheses-bold-italics.csl

pandoc --pdf-engine=/Library/TeX/texbin/pdflatex --number-sections --citeproc --bibliography nascent_rna_methods.bib -s README.md -o Methods_Ms.pdf --csl springer-basic-parentheses-bold-italics.csl
