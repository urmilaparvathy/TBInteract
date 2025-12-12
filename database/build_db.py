import pandas as pd
from sqlalchemy import create_engine
edges = pd.read_csv(snakemake.input[0])
engine = create_engine("sqlite:///{}".format(snakemake.output[0]))
edges.to_sql("edges", engine, if_exists="replace", index=False)
print("Wrote DB:", snakemake.output[0])
