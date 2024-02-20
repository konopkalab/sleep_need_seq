## Remove Ambient RNA using CellBender

```
## load module and activate conda
module load python=3.7
conda activate cellbender

## running cellbender
ls | grep -P "^[A|X][1-8]$" | xargs -I % -n 1 -P 24 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000 --total-droplets-included 50000'


## deactivate conda
conda deactivate
```
