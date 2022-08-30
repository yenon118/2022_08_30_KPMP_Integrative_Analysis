for i in seurat_clusters HumanPrimaryCellAtlasData BlueprintEncodeData MouseRNAseqData ImmGenData DatabaseImmuneCellExpressionData NovershternHematopoieticData MonacoImmuneData Experiment; do \
sbatch --time=0-02:00 --mem-per-cpu=30G -n 2 -N 1 \
--job-name=2022_08_30_markers_${i} \
--output=log_2022_08_30_markers_${i}_%j.out \
--wrap="cd /class/infoinst8005/chanye/2022_08_30_KPMP_Integrative_Analysis/scripts && \
source activate /class/infoinst8005/scRNA_seq_environment && \
Rscript 2022_08_30_markers.R -l ${i}"; \
done;