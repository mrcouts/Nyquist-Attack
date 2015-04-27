# adapta.sh
latex --src -interaction=nonstopmode [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms
bibtex [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms
makeindex [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms.nlo -s nomencl.ist -o [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms.nls
latex --src -interaction=nonstopmode [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms
# latex --src -interaction=nonstopmode [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms
# pdflatex -synctex=1 -interaction=nonstopmode [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms.tex
pdflatex -synctex=1 [ORSINO,\ COUTINHO,\ HESS-COELHO]\ Dynamic\ modelling\ and\ control\ of\ balanced\ parallel\ mechanisms.tex
# rm -f *.aux *.log 
rm -f *.blg *.dvi *.ps *.lot *.lof *.idx *.nlo *.nls *.ilg *.out *.toc *.bbl
