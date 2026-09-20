FROM mambaorg/micromamba:2.0.5

ARG NF_RNA_SOURCE_REVISION=unlabeled-development-image
LABEL org.opencontainers.image.title="nf-rna first-party execution image" \
      org.opencontainers.image.source="https://github.com/nf-rna/nf-rna" \
      org.opencontainers.image.revision="${NF_RNA_SOURCE_REVISION}"

WORKDIR /opt/rnaseq
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml environment.docker.yml pyproject.toml README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow ./workflow

USER root
RUN mkdir -p /opt/nf-rna/r && chown "$MAMBA_USER:$MAMBA_USER" /opt/nf-rna /opt/nf-rna/r
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER src/rnaseq/r/ /opt/nf-rna/r/

ENV PATH="/opt/conda/envs/rnaseq/bin:${PATH}"
ENV PYTHONUNBUFFERED=1

RUN micromamba create --yes --name rnaseq --file environment.yml && \
    micromamba install --yes --name rnaseq --file environment.docker.yml && \
    sh -c 'python -m pip install --no-deps . && test -f /opt/nf-rna/r/l1_analysis.R && command -v ps && ps --version && command -v python && command -v Rscript && Rscript -e "packages <- c(\"DESeq2\", \"tximport\", \"ggplot2\", \"pheatmap\", \"yaml\", \"jsonlite\", \"clusterProfiler\", \"AnnotationDbi\", \"org.Hs.eg.db\", \"org.Mm.eg.db\"); quit(status=if (all(vapply(packages, requireNamespace, logical(1), quietly=TRUE))) 0 else 1)" && python -c "import rnaseq.models, rnaseq.workflow_support; from rnaseq.workflow_assets import required_workflow_assets; assert all(path.is_file() for path in required_workflow_assets().values())" && python -m rnaseq.workflow_support report --help >/dev/null' && \
    micromamba clean --all --yes

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
