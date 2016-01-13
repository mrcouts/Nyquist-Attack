#relrenato
latex --src -interaction=nonstopmode relrenato 
#bibtex relrenato
makeindex relrenato
#makeindex relrenato.glo -s relrenato.ist -o relrenato.gls 
makeindex relrenato.cno -s relrenato.ist -o relrenato.cns 
makeindex relrenato.syo -s relrenato.ist -o relrenato.sys

latex --src -interaction=nonstopmode relrenato
#latex --src -interaction=nonstopmode relrenato 
pdflatex -synctex=1 relrenato.tex
#xelatex -synctex=1 -interaction=nonstopmode relrenato
#xelatex -synctex=1 relrenato
rm -f *.aux *.log *.blg *.dvi *.ps *.toc *.lot *.lof *.idx *.nlo *.nls *.ilg *.out *.fls
#
#relatorio_renato_maia_matarazzo_orsino
#pdftk relrenato2.pdf paperijmrs.pdf cat output relatorio_renato_maia_matarazzo_orsino_133192_2012-1.pdf
