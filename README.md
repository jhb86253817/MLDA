Master thesis: Training Algorithms for Multilingual Latent Dirichlet Allocation

Two algorithms are implemented for training MLDA under the standard and the proposed approximate training framework:
* collapsed Gibbs sampling
* variational inference

The Gibbs sampling part is developed from [GibbsLDA++] (http://gibbslda.sourceforge.net/).

Programming languages:
* python, for scripts
* c++, for main algorithms

Extra tools required:
* Stanford segmenter, for Chinese phrase segmentation
* opencc, for traditional-to-simplified Chinese transformation
