# Toolbox for Bioinformatic Tools

Since the microarray technology has continuously advanced, many techniques used to identify diagnostic markers by analyzing the genome-wide expression profiles have emerged. Gene expression is shown as colors on a microarray plate. These colors will be converted to numerical data representing gene expressions, and these numerical data will be used in identifying diagnostic markers. 
  
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
- GSE3494 Sample to Class Mapper
  - This program is used to match samples of dataset *[GSE3494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3494)* with their class.
- calculate.py
  - This is a file containing a bundle of functions used in calaulation processes.
- add_ons.py
  - This is a file containing functions which are used to deal with input data. 

## Prerequisite 
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

All required files must be prepared as shown in *file_format.pdf*

## Get Started 
1. Clone this project
```
git clone https://github.com/thatchayut/Bioinformatic_toolbox.git
```
2. Move into this project's direcory
```
cd Bioinformatic_toolbox
```
3. Install all required libraries as shown in the prerequisite
4. Place all required files in the **same** directory as the project's directory

### To use the *Bioinformatic Toolbox*
1. Execute this program
```
python toolbox.py
```
2. Answer all questions asked by the toolbox on your terminal
3. Waiting for the result. This step might take various ranges of time depending on the classification method.

### To use the *Gene Probe IDs to Gene Entrez IDs Converter*
1. Execute this program
```
python convert_to_entrez.py
```
2. Answer all questions asked by the toolbox on your terminal
3. Waiting for the result

### To use the *Feature Tester*
1. Execute this program
```
python feature_tester.py
```
2. This program requires manually configuration in feature_tester.py. Follow instructions provided by this program 
on your terminal.
3. Execute this program again
4. Answer all questions asked by the toolbox on your terminal
5. Waiting for the result

### To use the *Feature Extractor*
1. Execute this program
```
python feature_extractor.py
```
2. Answer all questions asked by the toolbox on your terminal
3. Waiting for the result

### To use the *Linear Discrimination Analysis*
1. Execute this program
```
python lda.py
```
2. This program requires manually configuration in feature_tester.py. Follow instructions provided by this program 
on your terminal.
3. Execute this program again
4. Waiting for the result

### To use the *GSE3494 Sample to Class Mapper*
1. Execute this program
```
python gse3494_sample_to_class_mapper.py
```
2. This program requires manually configuration in gse3494_sample_to_class_mapper.py. Follow instructions provided by this program 
on your terminal.
3. Execute this program again
4. Answer all questions asked by the toolbox on your terminal
5. Waiting for the result






  
