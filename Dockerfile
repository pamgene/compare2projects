FROM rocker/r-ver:4.3.0

LABEL maintainer="Dora Schuller <dschuller@pamgene.com>"

RUN apt-get update && apt-get install -y --no-install-recommends \
    # sudo \
    # libcurl4-gnutls-dev \
    libcurl4-openssl-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libproj-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxml2-dev \
    libmagick++-dev \
    pandoc \
    pandoc-citeproc \
    && rm -rf /var/lib/apt/lists/*

ENV RENV_VERSION 1.0.3

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

RUN addgroup --system app \
    && adduser --system --ingroup app app

WORKDIR /report

COPY . .

RUN chown app:app -R /home/app
USER app

RUN R -e 'renv::restore()'

EXPOSE 5151

CMD ["R", "-e", "shiny::runApp('/report', port = 5151, host = '0.0.0.0')"]
