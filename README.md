UnBlender allows respiratory scientists to perform cell type deconvolution in a validated manner. Not all cell type deconvolution approaches are feasible, especially when using highly granular (i.e. high resolution, very specific) cell type labels. UnBlender leverages the Human Lung Cell Atlas reference data to allow the user to specify a cell type label granularity that suits the needs of their research question, and then validates whether that approach is feasible.

(Why is this important? Because otherwise you run a very realistic risk of generating meaningless deconvolution results. For reliable results, it is essential to test whether the chosen combination of cell types *can* be deconvoluted in the chosen sample type.)

### Webtool
UnBlender is available as an online webtool at www.unblender.io, which can be used for deconvolution of nasal brush, bronchial brush, bronchial biopsy and parenchymal resection samples. No installation of software required, with a clickable user interface to configure analysis parameters.

### Command line interface
However, if you wish to deconvolute your bulk RNA-seq dataset on your own computer, this GitHub repository provides a pipeline for the command line interface (CLI).

Output of the CLI pipeline:
- A signature matrix that you can use to perform cell type deconvolution
- A per-cell type evaluation whether deconvolution results using this signature matrix are reliable in your sample type
- The reference data used to generate the signature matrix

## How to use:
1) Have CIBERSORTx (docker) installed, and activate it.
2) Install & activate the conda environment supplied with the code.
3) Configure your pipeline run using a config.yaml file (example file is provided, see below for parameters)
4) [Download the .h5ad HLCA file (core).](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)
4) In your command line interface, go to the directory that contains the Snakemake file.
5) Run the pipeline using the following command: 
	snakemake -c1 --configfile [config_file_name.yaml]

## Configuration of parameters:
- config_filename: the absolute path + file name to the config.yaml file, e.g. /home/user/Documents/my_deconvolution_analysis/config.yaml
- sample_type: "parenchyma", "bronchial_brush", "nasal_brush", or "bronchial_biopsy"
- output_dir: absolute path to the output directory, e.g.  /home/user/Documents/my_deconvolution_analysis/output
- DEG_file: optional, the absolute path + file name to a file with any gene names you want to remove from the analysis (e.g. because they are differentially expressed in your bulk RNA-seq cohort). If you don't want to use this feature, give it value "no_filter". It is *not* recommended to remove large numbers of genes, as this will reduce deconvolution accuracy. However, in specific cases it may be worthwhile to remove a few. Format: .csv/txt/tsv etc. file with a single gene identifier (HGNC gene names) on each new line.
- HLCA_file: path to the .h5ad HLCA file.
- email: the email address with which you run CIBERSORTx
- token: your private CIBERSORTx token
- cell_types: the cell types you wish to deconvolute, and which HLCA annotation level they correspond to. (See [HLCA supplementary table 4](https://www.nature.com/articles/s41591-023-02327-2#Sec57) for available cell types & levels.) Format: a nested YAML list. Sometimes it is useful to merge two cell types with similar gene expression profiles into one category, which is possible as shown in the example below for basal and secretory cells:

### Partial example cell type selection of a deconvolution setup into three categories:
```cell_types:
  - - "ann_level_2_clean"
    - "Lymphatic EC"

  - - "ann_level_3_clean"
    - - "Basal"
      - "Secretory"

  - - "ann_level_3_clean"
    - "Multiciliated lineage"
```

# UnBlender is in beta version
UnBlender pipeline code is in beta version. If you encounter any bugs, please contact unblender.info@gmail.com.
(A known bug is that currently, use of certain cell types with uncommon non-alphanumeric characters in their name may cause problems running the pipeline.)
