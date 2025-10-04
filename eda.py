import pandas as pd

df = pd.read_csv('neo_2025_10_04.csv')

for col in df.columns:
    print(col)