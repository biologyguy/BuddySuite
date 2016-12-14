#!/usr/bin/env sh
pdflatex MBE_article.tex;
bibtex MBE_article;
pdflatex MBE_article.tex;
pdflatex MBE_article.tex;
mv MBE_article.log MBE_article.aux MBE_article.pdf MBE_article.bbl MBE_article.blg MBE_article.idx compile_dir;
