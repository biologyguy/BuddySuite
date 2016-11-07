#!/usr/bin/env sh
pdflatex bmc_article.tex;
bibtex bmc_article;
pdflatex bmc_article.tex;
pdflatex bmc_article.tex;
mv bmc_article.log bmc_article.aux bmc_article.pdf bmc_article.bbl bmc_article.blg compile_dir;
