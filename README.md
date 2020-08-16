This repository contains a discriminatory modeling framework that is utilized downstream of the analysis found in [apis-methylation](https://github.com/littleblackfish/apis-methylation). 
As in the case of the other repository, it is not polished for public use but mainly documents the workflow for sake of reproducibility.

In conclusion of upstream analysis we observe that a substantial fraction of DNA methylation signal in the honey bee (and some other Hymenopterans) is robust to caste, developmental stage, tissue-of-origin and environmental history. 
Lacking dynamic response, we remark that these static epigenetic patterns can be efficiently carried by the underlying genetic sequence, 
and predict that they are strongly coupled to it. 
We show this to be the case through effectiveness of models that predict robustly methylated sites based on surrounding sequence.
We term this approach of deferring mechanistic explanation in favor of no-frills empirical support "*in-silico* assaying".

The primary task here is to predict a robustly methylated CpG site given a sequence context of a few hundred bases.
We set up logistic regression models trained on kmer compositions as a linear baseline, 
and show that shallow convolutional neural networks offer significantly improved performance, yielding over 0.93 area under precision-recall curves (AUPRC) for the honey bee.

A secondary task is to predict non-robustly methylated CpG sites, which still make up a very small subset of the genome.
We show that the same modeling approach yields reduced performance where AUPRC generally remains below 0.6, 
suggesting decreased coupling to sequence context for dynamically methylated CpG sites.

A tertiary task is to predict robustly methylated sites within the methylome, reinforcing distinct cis-regulation of robustly methylated sites.

A number of scripts are provided to encode output from [apis-methylation](https://github.com/littleblackfish/apis-methylation) in preparation for modeling.
In general, we start with pandas DataFrames that are indexed by seqid and position of all CpG sites in the genome, further annotated in terms of genomic features and methylation levels.
In summary: 
  * [extract_windows.py](scripts/extract_windows.py) extracts windows of a given radius for each CpG site from a reference genome into a multifasta file, optionally excluding overlapping windows. Methylation levels are json-encoded in fasta descriptions.
  * [encode_kmer.py](scripts/encode_kmer.py) prepares a kmer composition representation of a given multifasta file. 
  * [encode_onehot.py](scripts/encode_onehot.py) prepares a one-hot encoded representation of a given multifasta file. 

With sequences in place, those [scripts](scripts/) with the "train_" prefix execute training of various models for either task, with cross validation over chromosome based folds. 
Logistic models, being cheap, are utilized in largely disposable fashion while convolutional models are retained.
Initial hyperparameter sampling for convolutional models are stored in a [json file](models/models.json), while a [single set](models/models_one.json) was chosen for wider application.

Finally, performance of the resulting models are summarized in some Jupyter [notebooks](notebooks/).


