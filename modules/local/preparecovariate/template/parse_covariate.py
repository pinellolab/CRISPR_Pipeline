#!/usr/bin/env python

import json
import sys 
import pandas as pd

print(json.loads(sys.argv[1]))
dict = json.loads(sys.argv[1])
df = pd.DataFrame(dict)
df.to_csv('parse_covariate.csv', index=False)
