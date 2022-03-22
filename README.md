# properTy
Naive, memory, activation, effector, exhaustion properties of T cells


## Install
```
devtools::install_github('zgyaru/properTy')
```



## Example

### Step 1. Loading example data
```
library(properTy)

## It is a SingleCellExperiment toy dataset
data = loadExample()

```
### Step 2. Calculating immune property score
```
score = properTy(scMat = data@assays@data$counts, 
                 clusterNames = data.frame(data$clusterNames, row.names = colnames(data)))
```
### Step 3. Visulizing immune properties of interested T subSets
```
DensityProp(score, 
            c('CD8+ Tex-SPRY1', 'CD8+ Tex-XAF1', 'CD8+ Tex-MKI67')   ## interested cell cluster names
            )
```

<div align=center style="border:5px solid #000"><img  src="https://github.com/zgyaru/properTy/blob/main/pic/density.png"/> </div>


```
PieProp(score, 
        c('CD8+ Tex-SPRY1', 'CD8+ Tex-XAF1', 'CD8+ Tex-MKI67')   ## interested cell cluster names
        )
```

<div align=center style="border:5px solid #000"><img  
  src="https://github.com/zgyaru/properTy/blob/main/pic/pie.png"/> </div>
