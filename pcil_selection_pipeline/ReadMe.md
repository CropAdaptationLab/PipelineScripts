# PCIL Selection Pipeline

## Overview
This pipeline identifies:

- **PCIL (+)** → lines carrying introgressions overlapping target regions  
- **PCIL (–)** → matched non-carrier lines for comparison  

It supports two modes:
- **Full (genomics)** → best possible biological matches  
- **Available (operational)** → only lines with seed available  

This pipeline currently has three functions:
1) load_pcil_data()
2) select_pcil_positive()
3) select_pcil_negative()


## Files in this folder are:
1) load_pcil_data()_ReadMe.md, description of the load_pcil_data() function located here:
"https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/27fa2df7d1ed23ccc6faae951ac95e6426520fcf/load_pcil_data.R"

2) select_pcil_positive()_ReadMe.md, description of the select_pcil_positive() function located here:
"https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/e360eaa231b57a31077358d2766dc423dac93b3a/select_pcil_positive.R"

3) select_pcil_negative()_ReadMe.md, description of the select_pcil_negative() function located here:
"https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/51d1952f4936e0f55a817474b3d8da8e83ce2e56/select_pcil_negative.R"

4) Development_files, this folder contains readme for the additional files generated for the pipeline included in "load_pcil_data"
