# Q-value calculation

Q-values were calculated using the `qvalue` R package from nominal P values obtained from the statistical comparisons described in the Methods.

For each omics dataset, a CSV file containing one nominal P value per tested feature was prepared in advance. The same R script was then applied to the `p` column of each input file to calculate q-values.

The following datasets and statistical comparisons were included:

## Bulk RNA-seq

- Human PBMC:
  TPM values were used.
  HV1-5 versus MBC1-10 was tested using a two-tailed unpaired Welch's t-test.

- Mouse PBMC:
  TPM values were used.
  Sham1-4 versus 4T1-1-4 was tested using a two-tailed unpaired Welch's t-test.

- Mouse neutrophils:
  TPM values were used.
  Ly6C-low1-5 versus Ly6C-high1-5 was tested using a two-tailed paired t-test.

- Mouse liver:
  log2(TPM + 1) values were used.
  Sham-PBS1-4 versus 4T1-PBS1-7 and 4T1-PBS1-7 versus 4T1-Alb1-7 were each tested using a two-tailed unpaired Welch's t-test.

- Mouse skeletal muscle:
  log2(TPM + 1) values were used.
  Sham-PBS1-4 versus 4T1-PBS1-7 and 4T1-PBS1-7 versus 4T1-Alb1-7 were each tested using a two-tailed unpaired Welch's t-test.

## Plasma metabolomics

- Human plasma metabolome:
  HV1-5 versus MBC1-10 was tested using a two-tailed unpaired Welch's t-test.

- Mouse plasma metabolome:
  Sham1-3 versus 4T1-1-3 was tested using a two-tailed unpaired Welch's t-test.

## Proteomics

- Mouse neutrophil proteome:
  Ly6C-low1-5 versus Ly6C-high1-5 was tested using a two-tailed paired t-test.

## Input format

The script expects a CSV file containing a column named:

`p`

Each row should correspond to one tested feature, such as a gene, metabolite, or protein. Additional annotation columns may be included and will be retained in the output.

The input CSV files used in the study are not included in this repository. Nominal P values were obtained from the statistical comparisons described above and in the Methods. Please refer to the corresponding Dataset EV files for the nominal P values used in each analysis.

## Output

The script appends the following columns to the input table:

- `qvalue`
- `lfdr`

It also reports the estimated `pi0` value and saves a CSV file containing the original input together with the calculated q-values.

The q-value calculation script is intended to reproduce the multiple-testing correction step only. The nominal P values used as input were generated from the statistical comparisons described above and in the Methods; the corresponding nominal P values are provided in the relevant Dataset EV files.
