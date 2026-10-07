FROM python:3.15.0rc2-trixie

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir numpy pandas biopython vina

WORKDIR /data
