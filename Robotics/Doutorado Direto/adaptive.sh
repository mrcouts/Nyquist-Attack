# adapta.sh
latex --src -interaction=nonstopmode Artigo
bibtex Artigo
makeindex Artigo.nlo -s nomencl.ist -o Artigo.nls
latex --src -interaction=nonstopmode Artigo
# latex --src -interaction=nonstopmode Artigo
# pdflatex -synctex=1 -interaction=nonstopmode Artigo.tex
pdflatex -synctex=1 Artigo.tex
# rm -f *.aux *.log 
rm -f *.blg *.dvi *.ps *.lot *.lof *.idx *.nlo *.nls *.ilg *.out *.toc *.bbl
