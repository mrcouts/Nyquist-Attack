#relrenato
latex --src -interaction=nonstopmode relrenato 
bibtex relrenato
latex --src -interaction=nonstopmode relrenato
latex --src -interaction=nonstopmode relrenato 
pdflatex -synctex=1 relrenato.tex
rm -f *.aux *.log *.blg *.dvi *.ps *.toc *.lot *.lof *.idx *.nlo *.nls *.ilg *.out *.fls
#
#relatorio_renato_maia_matarazzo_orsino
#pdftk relrenato2.pdf paperijmrs.pdf cat output relatorio_renato_maia_matarazzo_orsino_133192_2012-1.pdf
