This folder contains the working files for the PhD.

These functions will all come together in the PhD to create:
1. a classifier to identify hypermutated tumour exome sequence samples
2. A classifier to assign tumour samples to subtype populations using a network processed network of genes*samples.



network_informed_clustering_function.R
This file contains functions for processing a gene*sample network of TRUE FALSE.
TRUE=at least 1 non-silent mutation in gene i in sample j.


create_geneAdj_and_mutM.R
This file creates the geneAdj and mutM matrices needed for the network clustering function above.
The program takes as input a TCGA .maf file and a .simple network file. 
The ".maf" is a vcf like exome sequence file. 
The ".simple" file is a network file in the following format:

vertex1	vertex2
vertex1	vertex3
vertex2	vertex3

Each line represents an edge in the network and each column represents the vertices incident upon the edge.


identify_hypermutated_samples.R
This file is the basis of a function to identify hypermutated tumour samples. This will be developed in to a more complex classifier function.


create_metadata_structure.R
This function is designed to create a unified datatable created from the "Biotab" clinical data tables. It creates a samples*variables table.
It is important to makes sure the tables included have unique ids for each sample rather than ids for numerous tissue samples for each sample.

