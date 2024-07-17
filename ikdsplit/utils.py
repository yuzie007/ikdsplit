import pandas as pd


def format_df(df: pd.DataFrame) -> pd.DataFrame:
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(12)
    return df
