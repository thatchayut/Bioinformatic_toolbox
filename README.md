# Toolbox for Bioinformatic Tools

Since the microarray technology has continuously advanced, many techniques used to identifydiagnostic markers by analyzing the genome-wide expression profiles have emerged. Gene expressionis shown as colors on a microarray plate. These colors will be converted to numerical data representing gene expressions, and these numerical data will be used in identifying diagnostic markers. 
  
Many research studies show that using **_pathways_** (or groups) of genes that interact with each other in molecular level to identify genes that are used to classify an occurence of a specific disease yields more *reliability* and *reproducibility* compared with using genes themselves to identify diagnostic markers.

This toolbox offers multiple existing pathway-based classification methods including *Mean*, *[Condition-Responsive Genes (CORGs)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2563693/)* , and *[Log-Likelihood Ratio (LLR)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2781165/)*. Moreover, a gene-based classification method is also provided in this toolbox.

Other basic tools are also provided in this toolbox. These tools can be used to assist this toolbox including
- Gene Probe IDs to Gene Entrez IDs Converter
  - This program is used to convert gene probe IDs to gene entrez IDs.
- Feature Tester
  - This program is used to test generated feature set.
- Feature Extractor
  - This program is used to extract features into member genes.
- Linear Discrimination Analysis
  - This program is created as a prototype to be used as one of important functions in a toolbox.

## Prerequisite ##
These libraries are required for using this toolbox.
- pandas
- random
- math
- time
- scipy
- numpy
- copy
- sklearn

These libraries are required for using other basic tools.
- xlsxwriter
- collections
- xlwt

# Get Started

To use the **Bioinformatic Toolbox**

1. Clone this project
```
git clone https://github.com/thatchayut/Bioinformatic_toolbox.git
```
2. Move into this project's direcory
3. Install all required libraries as shown in the prerequisite.
4. Execute the toolbox
```
python toolbox.py
```
**_Note_ :** *All files asked by this toolbox **must** be moved into the project's directory.* 
5. Answer all questions asked by the toolbox on your terminal
6. Wating for the result. This step might take a lot of time depending on the classification method.

  
