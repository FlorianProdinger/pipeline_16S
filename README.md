# 18S rRNA gene amplicon analysis pipeline

This pipeline can be used as an example of how a pipeline can be built. It's not aparticularly good pipeline, but it runs.
I found it difficult to setup the metadata sheet and the qiime manifest. Feel free to use the reference of this pipeline.

It can be used with and without the qsub system, but QIIME2 and R need to be installed.

If you want to use the pipeline yourself keep in mind that the following variabels need to be changed:
- ```OUT_DIR```
- ```MANIFEST```
- ```METADATA_FILE```
- ```SILVA_FILE_SEQ_97```
- ```SILVA_FILE_97```
- path to the R script "command_rename_first_line_taxonomy_file.R"



