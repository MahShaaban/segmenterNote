# segmenterNote

A technical note to describe the `segmenter` package.
- This note is goint to be presented as an elevator pitch during [useR 2021](https://user2021.r-project.org/)
- This note is based on the template [technical note](https://github.com/useRconf/templates/tree/main/technical_note)

## Create the technical note

```bash
# build a docker image with the required packages
docker build -t segmenternote .

# render the rmd files while keeping md and files
docker run --rm -v $(pwd):/home segmenternote Rscript -e "rmarkdown::render('/home/technical_note.Rmd', clean = FALSE)"
```
