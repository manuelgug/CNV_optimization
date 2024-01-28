# Copy Number Variant Optimization

## Objective

The objective of this script is to optimize the selection of amplicons (from ampseq) for the calculation of copy number variants (CNVs) in _Plasmodium falciparum_ through the use of a Generalized Additive Model (GAM). To achieve this, a genetic algorithm (GA) is employed to systematically identify the most informative set of amplicons from a larger pool.

## Functionality Overview

1. **Input Data:**
   - **Loci of Interest for Reference:** Load reference loci data from an RDS file.
   - **Amplicon Coverage Data:** Provide the path to the amplicon coverage file.

2. **Data Formatting:**
   - Use the `formating_ampCov` function to format the amplicon coverage data, filtering out irrelevant controls and creating additional columns for analysis.

3. **CNV Estimation:**
   - Employ the `estCNV` function to estimate copy number variants for each sample. The GAM model takes into account amplicon length and specified loci of interest.

4. **Genetic Algorithm:**
   - Utilize a genetic algorithm (`GA`) to optimize the selection of amplicons for each sample. The GA dynamically evolves a population of potential solutions to maximize the accuracy of CNV predictions known CNV fold changes in control strains of _P. falciparum_.

5. **Fitness Function:**
   - Define a fitness function that evaluates the performance of each potential solution (set of amplicons) based on the expected and observed fold changes of loci of interest.

6. **GA Parameters:**
   - Set parameters such as population size, generations, mutation probability, and elitism to guide the genetic algorithm's search for optimal solutions.

7. **Run GA:**
   - Execute the GA with real-time plotting to visualize the evolution of fitness scores over generations.

8. **Extract Results:**
   - Extract the optimal set of amplicons identified by the genetic algorithm.

9. **Fold Change Calculation:**
   - Calculate fold changes for each locus of interest using the optimal amplicon set.
