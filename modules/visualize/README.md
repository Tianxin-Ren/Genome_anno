# Visualize Module - Genome_anno Snakemake Pipeline

The `visualize` module provides automated visualizations and generates an interactive, high-quality HTML dashboard report summarizing abundance patterns, genome quality metrics, and taxonomic profiles across your samples.

## Configuration

The module is highly customizable. You can control which plots are generated and configure global report styling.

### 1. Default Configuration
The standard configuration file is located at:
- [config/visualize.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/config/visualize.yaml)

### 2. Standalone config_vi Folder
For dedicated visualization runs, use the configurations in the standalone folder:
- [config_vi/pipeline.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/config_vi/pipeline.yaml)
- [config_vi/visualize.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/config_vi/visualize.yaml)
- [config_vi/samples.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/config_vi/samples.yaml)

## Parameters

| Parameter | Default | Options | Description |
| :--- | :--- | :--- | :--- |
| `visualize.enable` | `false` | `true`, `false` | Set to `true` to activate the module |
| `visualize.results_dir` | `"results/visualize"` | Any valid path | Output directory for all visualization files |
| `visualize.plots.abundance_heatmap.enable` | `true` | `true`, `false` | Enable abundance heatmap generation |
| `visualize.plots.abundance_heatmap.palette` | `"viridis"` | `"viridis"`, `"plasma"`, `"inferno"`, `"coolwarm"` | Color scheme for abundance mapping |
| `visualize.plots.quality_assessment.enable` | `true` | `true`, `false` | Enable genome completeness/contamination scatter plot |
| `visualize.plots.taxonomic_profile.enable` | `true` | `true`, `false` | Enable taxonomic classification profiles |
| `visualize.report.theme` | `"dark"` | `"dark"`, `"light"` | Base color scheme for the interactive report |

## Outputs

All outputs will be saved inside the configured `results_dir` (default: `results/visualize`):

- **Summary Dashboard**: `summary_report.html` (Comprehensive interactive HTML report with charts and metrics)
- **Abundance Heatmap**: `abundance/abundance_heatmap.pdf` (Heatmap of genome/gene abundance)
- **Quality Scatter Plot**: `quality/quality_scatter_plot.png` (CheckM completeness vs contamination scatter plot)
- **Taxonomic Distribution**: `taxonomy/taxonomic_distribution.html` (GTDB-Tk taxonomic profile)

## Defensive Fallback Mechanism

If the upstream analysis rules (e.g. abundance quantification or binning) have not been run or are disabled, the visualization engine will **automatically and gracefully fallback to generating highly-realistic demonstration data**. This allows you to compile and test the visualization pipeline standalone or inspect layout formats before launching full genomic pipelines.

## Execution

To run the visualization module standalone using the dedicated `config_vi` directory:

```bash
# Dry-run validation
snakemake --config conf=config_vi -n

# Standard execution
snakemake --config conf=config_vi -j 8
```
