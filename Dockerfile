FROM mambaorg/micromamba:2.0.5

WORKDIR /opt/rnaseq
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml pyproject.toml README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src

RUN micromamba create --yes --name rnaseq --file environment.yml && \
    micromamba run --name rnaseq sh -c 'command -v python && command -v Rscript && Rscript -e "packages <- c(\"DESeq2\", \"tximport\", \"ggplot2\", \"pheatmap\", \"yaml\", \"jsonlite\", \"clusterProfiler\", \"AnnotationDbi\", \"org.Hs.eg.db\", \"org.Mm.eg.db\"); quit(status=if (all(vapply(packages, requireNamespace, logical(1), quietly=TRUE))) 0 else 1)"' && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/envs/rnaseq/bin:$PATH
ENV PYTHONUNBUFFERED=1

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
