# Syze

An Unequal Caching Strategy for Graph Analytics

Syze is developed from [GPOP](https://github.com/souravpati/GPOP), by adding two unique features:
1. Unequally partitioning subgraphs for workload balance (e.g., `dynamicSplit`)
2. Encoding vertex ID to facilitate data propagation (e.g., `locateSub`)

The overhead introduced by Syze is minor, but the acceleration on imbalanced power-law graphs is remarkable.

The graph [converter](https://github.com/yuang-chen/Graph-Converter) contains the tools for generating the needed CSR graphs.

## How to run:
```
make

./app/app <filename> -v <root vertex (optional)>
                     -s <size(default 256KB)>
                     -i <#iterations(default 20)>
                     -r <#rounds (default 3)> 
                     -t <#threads (default 20)>
                     -y <is_dynamic [1|0] (default 0)>
```

## ToDo
I am working on rewritting Syze based on modern C++ design. Thus, the codes in this repository might be flushed in (hopefully near) future. 