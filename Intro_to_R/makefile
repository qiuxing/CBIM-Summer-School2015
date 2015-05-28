TEX = pdflatex
# LATEXFILES = presentation article
# TEXFLAGS = 
SRC	:= $(shell egrep -l '^[^%]*\\begin\{document\}' *.tex)
TRG	= $(SRC:%.tex=%.dvi)
PSF	= $(SRC:%.tex=%.ps)
PDF	= $(SRC:%.tex=%.pdf)

all : presentation.pdf article.pdf

presentation.pdf : presentation.tex
	$(TEX) $(TEXFLAGS) presentation.tex
article.pdf : article.tex
	$(TEX) $(TEXFLAGS) article.tex

article.pdf presentation.pdf : main.tex 

clean	:
	rm *.aux *.log *.snm *.toc *.nav *.out