# Microarray Gene Expression Analysis with Python
Ahmed Moustafa
October 5, 2023

- [<span class="toc-section-number">1</span> Importing required
  libraries](#importing-required-libraries)
- [<span class="toc-section-number">2</span> Loading the
  data](#loading-the-data)
- [<span class="toc-section-number">3</span> Checking the data
  behavior](#checking-the-data-behavior)
  - [<span class="toc-section-number">3.1</span>
    Transforming](#transforming)
- [<span class="toc-section-number">4</span> Exploring the
  data](#exploring-the-data)
  - [<span class="toc-section-number">4.1</span> Boxplot](#boxplot)
  - [<span class="toc-section-number">4.2</span> Hierarchical
    clustering](#hierarchical-clustering)
- [<span class="toc-section-number">5</span> Biological
  signifiance](#biological-signifiance)
  - [<span class="toc-section-number">5.1</span> Slicing the dataset by
    condition](#slicing-the-dataset-by-condition)
  - [<span class="toc-section-number">5.2</span> Gene-wise mean
    expression](#gene-wise-mean-expression)
  - [<span class="toc-section-number">5.3</span> Scatter plot using
    Matplotlib](#scatter-plot-using-matplotlib)
  - [<span class="toc-section-number">5.4</span> Calculating the
    fold-change](#calculating-the-fold-change)
  - [<span class="toc-section-number">5.5</span> Scatter plot using
    Seaborn](#scatter-plot-using-seaborn)
  - [<span class="toc-section-number">5.6</span> Historgram of the
    fold-change](#historgram-of-the-fold-change)
- [<span class="toc-section-number">6</span> Statistical
  significance](#statistical-significance)
  - [<span class="toc-section-number">6.1</span> Calculating t-test
    *p*-value](#calculating-t-test-p-value)
  - [<span class="toc-section-number">6.2</span> Historgram of the
    *p*-value](#historgram-of-the-p-value)
- [<span class="toc-section-number">7</span> Biological & statistical
  signifiance](#biological-statistical-signifiance)
  - [<span class="toc-section-number">7.1</span> Volcano Plot (*p*-value
    vs. fold-change)](#volcano-plot-p-value-vs.-fold-change)
- [<span class="toc-section-number">8</span> Differentially expressed
  genes (DEGs)](#differentially-expressed-genes-degs)
  - [<span class="toc-section-number">8.1</span> Genes with significant
    fold-change](#genes-with-significant-fold-change)
  - [<span class="toc-section-number">8.2</span> Genes with significant
    *p*-value](#genes-with-significant-p-value)
  - [<span class="toc-section-number">8.3</span> Genes with significant
    fold-change & significant
    *p*-value](#genes-with-significant-fold-change-significant-p-value)
  - [<span class="toc-section-number">8.4</span> Heatmap](#heatmap)

## Importing required libraries

``` python
import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import seaborn as sns
```

## Loading the data

``` python
data = pd.read_table("https://media.githubusercontent.com/media/ahmedmoustafa/gene-expression-datasets/main/datasets/irf6/irf6.tsv", index_col='ID')
data.head()
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
&#10;    .dataframe tbody tr th {
        vertical-align: top;
    }
&#10;    .dataframe thead th {
        text-align: right;
    }
</style>

|              | KO1     | KO2     | KO3     | WT1     | WT2     | WT3     |
|--------------|---------|---------|---------|---------|---------|---------|
| ID           |         |         |         |         |         |         |
| 1415670_at   | 6531.0  | 5562.8  | 6822.4  | 7732.1  | 7191.2  | 7551.9  |
| 1415671_at   | 11486.3 | 10542.7 | 10641.4 | 10408.2 | 9484.5  | 7650.2  |
| 1415672_at   | 14339.2 | 13526.1 | 14444.7 | 12936.6 | 13841.7 | 13285.7 |
| 1415673_at   | 3156.8  | 2219.5  | 3264.4  | 2374.2  | 2201.8  | 2525.3  |
| 1415674_a_at | 4002.0  | 3306.9  | 3777.0  | 3760.6  | 3137.0  | 2911.5  |

</div>

``` python
number_of_genes = data.shape[0]
number_of_genes
```

    45101

``` python
ids = data.index # The ids of the genes are the names of the rows
ids
```

    Index(['1415670_at', '1415671_at', '1415672_at', '1415673_at', '1415674_a_at',
           '1415675_at', '1415676_a_at', '1415677_at', '1415678_at', '1415679_at',
           ...
           'AFFX-r2-P1-cre-5_at', 'AFFX-ThrX-3_at', 'AFFX-ThrX-5_at',
           'AFFX-ThrX-M_at', 'AFFX-TransRecMur/X57349_3_at',
           'AFFX-TransRecMur/X57349_5_at', 'AFFX-TransRecMur/X57349_M_at',
           'AFFX-TrpnX-3_at', 'AFFX-TrpnX-5_at', 'AFFX-TrpnX-M_at'],
          dtype='object', name='ID', length=45101)

## Checking the data behavior

Check the behavior of the data (e.g., normal?, skewed?)

``` python
data.hist()
```

    array([[<Axes: title={'center': 'KO1'}>, <Axes: title={'center': 'KO2'}>],
           [<Axes: title={'center': 'KO3'}>, <Axes: title={'center': 'WT1'}>],
           [<Axes: title={'center': 'WT2'}>, <Axes: title={'center': 'WT3'}>]],
          dtype=object)

![](Irf6_py_files/figure-commonmark/cell-6-output-2.png)

### Transforming

$log2$ transformation

``` python
data2 = np.log2(data)
data2.head()
data2.hist()
```

    array([[<Axes: title={'center': 'KO1'}>, <Axes: title={'center': 'KO2'}>],
           [<Axes: title={'center': 'KO3'}>, <Axes: title={'center': 'WT1'}>],
           [<Axes: title={'center': 'WT2'}>, <Axes: title={'center': 'WT3'}>]],
          dtype=object)

![](Irf6_py_files/figure-commonmark/cell-7-output-2.png)

## Exploring the data

### Boxplot

#### Using Matplotlib

``` python
plt.boxplot(data2)
```

    {'whiskers': [<matplotlib.lines.Line2D at 0x13fa01690>,
      <matplotlib.lines.Line2D at 0x13fa00450>,
      <matplotlib.lines.Line2D at 0x14c2d4410>,
      <matplotlib.lines.Line2D at 0x14c2d7f90>,
      <matplotlib.lines.Line2D at 0x13fe45250>,
      <matplotlib.lines.Line2D at 0x13fe46990>,
      <matplotlib.lines.Line2D at 0x13fa35990>,
      <matplotlib.lines.Line2D at 0x13fa36850>,
      <matplotlib.lines.Line2D at 0x175147610>,
      <matplotlib.lines.Line2D at 0x175145750>,
      <matplotlib.lines.Line2D at 0x13ffb97d0>,
      <matplotlib.lines.Line2D at 0x13ffba2d0>],
     'caps': [<matplotlib.lines.Line2D at 0x13fa003d0>,
      <matplotlib.lines.Line2D at 0x13fa02f90>,
      <matplotlib.lines.Line2D at 0x14c2d4310>,
      <matplotlib.lines.Line2D at 0x14c2d45d0>,
      <matplotlib.lines.Line2D at 0x13fe47090>,
      <matplotlib.lines.Line2D at 0x13fe44fd0>,
      <matplotlib.lines.Line2D at 0x13fa37f10>,
      <matplotlib.lines.Line2D at 0x13fa34ed0>,
      <matplotlib.lines.Line2D at 0x175147d90>,
      <matplotlib.lines.Line2D at 0x175144090>,
      <matplotlib.lines.Line2D at 0x13ffb9c10>,
      <matplotlib.lines.Line2D at 0x13ffba9d0>],
     'boxes': [<matplotlib.lines.Line2D at 0x14c0d2dd0>,
      <matplotlib.lines.Line2D at 0x14c2d4a90>,
      <matplotlib.lines.Line2D at 0x14c2d5cd0>,
      <matplotlib.lines.Line2D at 0x13fa34290>,
      <matplotlib.lines.Line2D at 0x175147390>,
      <matplotlib.lines.Line2D at 0x13ffb9950>],
     'medians': [<matplotlib.lines.Line2D at 0x13fa01c90>,
      <matplotlib.lines.Line2D at 0x14c2d4990>,
      <matplotlib.lines.Line2D at 0x13fe47590>,
      <matplotlib.lines.Line2D at 0x13fa36550>,
      <matplotlib.lines.Line2D at 0x107d95390>,
      <matplotlib.lines.Line2D at 0x13ffbba50>],
     'fliers': [<matplotlib.lines.Line2D at 0x13fc20390>,
      <matplotlib.lines.Line2D at 0x14cd3a7d0>,
      <matplotlib.lines.Line2D at 0x13fe45610>,
      <matplotlib.lines.Line2D at 0x175144910>,
      <matplotlib.lines.Line2D at 0x13fcbe950>,
      <matplotlib.lines.Line2D at 0x13fe45950>],
     'means': []}

![](Irf6_py_files/figure-commonmark/cell-8-output-2.png)

#### Using Pandas

``` python
data2.boxplot()
```

    <Axes: >

![](Irf6_py_files/figure-commonmark/cell-9-output-2.png)

#### Using Seaborn

``` python
sns.boxplot(data2)
```

    <Axes: >

![](Irf6_py_files/figure-commonmark/cell-10-output-2.png)

### Hierarchical clustering

``` python
linkage_matrix = hierarchy.linkage(data2.T, method='ward')  # Transpose data with .T
hierarchy.dendrogram(linkage_matrix, labels = data2.columns)
```

    {'icoord': [[15.0, 15.0, 25.0, 25.0],
      [5.0, 5.0, 20.0, 20.0],
      [45.0, 45.0, 55.0, 55.0],
      [35.0, 35.0, 50.0, 50.0],
      [12.5, 12.5, 42.5, 42.5]],
     'dcoord': [[0.0, 219.95424812831652, 219.95424812831652, 0.0],
      [0.0, 238.48396280958053, 238.48396280958053, 219.95424812831652],
      [0.0, 231.7745678766896, 231.7745678766896, 0.0],
      [0.0, 243.84222507735646, 243.84222507735646, 231.7745678766896],
      [238.48396280958053,
       281.5274405240738,
       281.5274405240738,
       243.84222507735646]],
     'ivl': ['KO2', 'KO1', 'KO3', 'WT1', 'WT2', 'WT3'],
     'leaves': [1, 0, 2, 3, 4, 5],
     'color_list': ['C0', 'C0', 'C0', 'C0', 'C0'],
     'leaves_color_list': ['C0', 'C0', 'C0', 'C0', 'C0', 'C0']}

![](Irf6_py_files/figure-commonmark/cell-11-output-2.png)

## Biological signifiance

### Slicing the dataset by condition

``` python
ko = data2[['KO1', 'KO2', 'KO3']] # KO dataframe (KO1,KO2,KO3)
ko.head()
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
&#10;    .dataframe tbody tr th {
        vertical-align: top;
    }
&#10;    .dataframe thead th {
        text-align: right;
    }
</style>

|              | KO1       | KO2       | KO3       |
|--------------|-----------|-----------|-----------|
| ID           |           |           |           |
| 1415670_at   | 12.673088 | 12.441596 | 12.736064 |
| 1415671_at   | 13.487627 | 13.363957 | 13.377400 |
| 1415672_at   | 13.807677 | 13.723458 | 13.818253 |
| 1415673_at   | 11.624247 | 11.116019 | 11.672602 |
| 1415674_a_at | 11.966505 | 11.691264 | 11.883025 |

</div>

``` python
wt = data2[['WT1', 'WT2', 'WT3']] # WT dataframe (WT1,WT2,WT3)
wt.head()
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
&#10;    .dataframe tbody tr th {
        vertical-align: top;
    }
&#10;    .dataframe thead th {
        text-align: right;
    }
</style>

|              | WT1       | WT2       | WT3       |
|--------------|-----------|-----------|-----------|
| ID           |           |           |           |
| 1415670_at   | 12.916645 | 12.812017 | 12.882624 |
| 1415671_at   | 13.345433 | 13.211356 | 12.901282 |
| 1415672_at   | 13.659171 | 13.756734 | 13.697587 |
| 1415673_at   | 11.213226 | 11.104468 | 11.302239 |
| 1415674_a_at | 11.876747 | 11.615170 | 11.507547 |

</div>

### Gene-wise mean expression

Note: the mean function can take the `axis` parameter to determine the
direction of computing the mean, where:

- `axis=0` → vertical (by column), the default direction
- `axis=1` → horizontal (by row), the direction that we want in this
  case.

``` python
ko_means = ko.mean(axis=1) # Compute the means of the KO samples
ko_means.head()
```

    ID
    1415670_at      12.616916
    1415671_at      13.409661
    1415672_at      13.783129
    1415673_at      11.470956
    1415674_a_at    11.846931
    dtype: float64

``` python
wt_means = wt.mean(axis=1) # Compute the means of the WT samples
wt_means.head()
```

    ID
    1415670_at      12.870428
    1415671_at      13.152690
    1415672_at      13.704497
    1415673_at      11.206644
    1415674_a_at    11.666488
    dtype: float64

### Scatter plot using Matplotlib

``` python
plt.scatter(x = wt_means, y = ko_means)
```

    <matplotlib.collections.PathCollection at 0x13cf26110>

![](Irf6_py_files/figure-commonmark/cell-16-output-2.png)

### Calculating the fold-change

``` python
fold_change = ko_means - wt_means # The difference between means
fold_change.head()
```

    ID
    1415670_at     -0.253513
    1415671_at      0.256971
    1415672_at      0.078632
    1415673_at      0.264312
    1415674_a_at    0.180443
    dtype: float64

Note:

- **+ve** fold-change → **Up**-regulation ↑
- **-ve** fold-change → **Down**-regulation ↓

### Scatter plot using Seaborn

``` python
sns.scatterplot(x = wt_means, y = ko_means, hue = fold_change, palette='RdBu_r')
plt.xlabel('WT')
plt.ylabel('KO')
plt.title('Irf6 KO vs. WT Gene Expression')
```

    Text(0.5, 1.0, 'Irf6 KO vs. WT Gene Expression')

![](Irf6_py_files/figure-commonmark/cell-18-output-2.png)

### Historgram of the fold-change

#### Using Matplotlib

``` python
plt.hist (fold_change)
```

    (array([3.0000e+00, 8.0000e+00, 6.0000e+00, 3.2000e+01, 4.0100e+02,
            1.1552e+04, 3.1788e+04, 1.2640e+03, 3.4000e+01, 1.3000e+01]),
     array([-11.17862753,  -9.36607657,  -7.55352561,  -5.74097464,
             -3.92842368,  -2.11587271,  -0.30332175,   1.50922922,
              3.32178018,   5.13433114,   6.94688211]),
     <BarContainer object of 10 artists>)

![](Irf6_py_files/figure-commonmark/cell-19-output-2.png)

#### Using Seaborn

``` python
sns.histplot (fold_change)
```

    <Axes: ylabel='Count'>

![](Irf6_py_files/figure-commonmark/cell-20-output-2.png)

## Statistical significance

### Calculating t-test *p*-value

(a *p*-value for each gene i.e., horizontally)

``` python
t_stat, p_value = stats.ttest_ind(ko, wt, axis=1)
t_stat_df = pd.DataFrame({'t_stat': t_stat, 'p_value': p_value})
t_stat_df.head()
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
&#10;    .dataframe tbody tr th {
        vertical-align: top;
    }
&#10;    .dataframe thead th {
        text-align: right;
    }
</style>

|     | t_stat    | p_value  |
|-----|-----------|----------|
| 0   | -2.677586 | 0.055367 |
| 1   | 1.872445  | 0.134450 |
| 2   | 1.904526  | 0.129561 |
| 3   | 1.413610  | 0.230364 |
| 4   | 1.321077  | 0.256980 |

</div>

### Historgram of the *p*-value

#### Using Matplotlib

``` python
plt.hist (-np.log10(t_stat_df['p_value']))
```

    (array([3.0659e+04, 9.2820e+03, 3.2870e+03, 1.2330e+03, 4.3300e+02,
            1.3700e+02, 5.4000e+01, 1.1000e+01, 3.0000e+00, 2.0000e+00]),
     array([2.47714579e-06, 7.35796647e-01, 1.47159082e+00, 2.20738499e+00,
            2.94317915e+00, 3.67897332e+00, 4.41476749e+00, 5.15056166e+00,
            5.88635583e+00, 6.62215000e+00, 7.35794417e+00]),
     <BarContainer object of 10 artists>)

![](Irf6_py_files/figure-commonmark/cell-22-output-2.png)

#### Using Seaborn

``` python
sns.histplot (-np.log10(t_stat_df['p_value']))
```

    <Axes: xlabel='p_value', ylabel='Count'>

![](Irf6_py_files/figure-commonmark/cell-23-output-2.png)

## Biological & statistical signifiance

### Volcano Plot (*p*-value vs. fold-change)

#### Using Matplotlib

``` python
plt.scatter (x = fold_change, y = -np.log10(t_stat_df['p_value']))
```

    <matplotlib.collections.PathCollection at 0x13c951190>

![](Irf6_py_files/figure-commonmark/cell-24-output-2.png)

#### Using Seaborn

``` python
sns.scatterplot (x = fold_change.values, y = -np.log10(t_stat_df['p_value']), hue = fold_change.values, palette='RdBu_r')
```

    <Axes: ylabel='p_value'>

![](Irf6_py_files/figure-commonmark/cell-25-output-2.png)

## Differentially expressed genes (DEGs)

### Genes with significant fold-change

``` python
fold_change_cutoff = 2
np.sum(abs(fold_change.values) >= fold_change_cutoff)
```

    1051

### Genes with significant *p*-value

``` python
pvalue_cutoff = 0.001
np.sum(t_stat_df['p_value'] <= pvalue_cutoff)
```

    575

### Genes with significant fold-change & significant *p*-value

``` python
np.sum((abs(fold_change.values) >= fold_change_cutoff) & (t_stat_df['p_value'] <= pvalue_cutoff))
```

    163

``` python
filtered = data2.reset_index().loc[(abs(fold_change.values) >= fold_change_cutoff) & (t_stat_df['p_value'] <= pvalue_cutoff)]
filtered.set_index('ID', inplace=True)
filtered.shape
filtered.head()
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }
&#10;    .dataframe tbody tr th {
        vertical-align: top;
    }
&#10;    .dataframe thead th {
        text-align: right;
    }
</style>

|              | KO1       | KO2       | KO3       | WT1       | WT2       | WT3       |
|--------------|-----------|-----------|-----------|-----------|-----------|-----------|
| ID           |           |           |           |           |           |           |
| 1416200_at   | 13.312004 | 12.973357 | 12.868456 | 7.404290  | 8.558803  | 8.683696  |
| 1416236_a_at | 14.148397 | 14.039236 | 14.130007 | 12.236044 | 12.022402 | 11.495056 |
| 1417808_at   | 5.321928  | 5.442943  | 4.053111  | 15.169780 | 15.070087 | 14.753274 |
| 1417932_at   | 10.602884 | 10.257152 | 10.496055 | 13.984454 | 14.203294 | 13.720960 |
| 1418050_at   | 10.622052 | 10.975490 | 10.795066 | 12.865134 | 13.012048 | 12.658122 |

</div>

### Heatmap

``` python
sns.clustermap(filtered, cmap='RdYlGn_r', standard_scale = 0)
```

![](Irf6_py_files/figure-commonmark/cell-30-output-1.png)
