docker run -i --rm --name texlive -v "$PWD":/home -w /home texlive/texlive:TL2020-historic pdflatex elsi_manual.tex
docker run -i --rm --name texlive -v "$PWD":/home -w /home texlive/texlive:TL2020-historic bibtex elsi_manual.aux
docker run -i --rm --name texlive -v "$PWD":/home -w /home texlive/texlive:TL2020-historic pdflatex elsi_manual.tex
docker run -i --rm --name texlive -v "$PWD":/home -w /home texlive/texlive:TL2020-historic pdflatex elsi_manual.tex
rm -f elsi_manual.aux elsi_manual.bbl elsi_manual.blg elsi_manual.log elsi_manual.out elsi_manual.toc
