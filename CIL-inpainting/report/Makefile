OUTDIR = build
LATEX_FLAGS=-interaction=nonstopmode -output-directory=$(OUTDIR)

all: report.pdf

report.pdf: report.tex report.bib $(shell find graphics)
	mkdir -p build
	pdflatex $(LATEX_FLAGS) report.tex > /dev/null
	cp -f report.bib build && cp -f IEEEtran.bst $(OUTDIR) && cd build && bibtex report > /dev/null
	pdflatex $(LATEX_FLAGS) report.tex > /dev/null
	pdflatex $(LATEX_FLAGS) report.tex > /dev/null
	
debug: report.tex $(shell find graphics)
	pdflatex $(LATEX_FLAGS) report.tex
	cp -f report.bib build && cp -f IEEEtran.bst $(OUTDIR) && cd build && bibtex report
	pdflatex $(LATEX_FLAGS) report.tex
	pdflatex $(LATEX_FLAGS) report.tex

.PHONY: all clean cleanall spell
	
all: report.pdf

clean:
	cd build && rm -rf *.bbl *.blg *.aux *.log *~ *.bak *.ps *.dvi *.log *.out *.tmp 

cleanall:
	cd build && rm -rf *.bbl *.blg *.aux *.log *~ *.bak *.ps *.dvi *.log *.pdf svnver.tex *.out *.tmp *.bst *.bib

spell:
	ispell -f ispell.words -t *.tex
