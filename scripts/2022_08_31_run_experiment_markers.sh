
for i in HumanPrimaryCellAtlasData_Experiment BlueprintEncodeData_Experiment MouseRNAseqData_Experiment ImmGenData_Experiment DatabaseImmuneCellExpressionData_Experiment NovershternHematopoieticData_Experiment MonacoImmuneData_Experiment; do \
sbatch --time=0-02:00 --mem-per-cpu=30G -n 2 -N 1 \
--job-name=2022_08_31_experiment_markers_${i} \
--output=log_2022_08_31_experiment_markers_${i}_%j.out \
--wrap="cd /class/infoinst8005/chanye/2022_08_30_KPMP_Integrative_Analysis/scripts && \
source activate /class/infoinst8005/scRNA_seq_environment && \
Rscript 2022_08_31_experiment_markers.R -l ${i}"; \
done;
