import sqlite3
import pandas as pd


def get_all_bad_regions_list() -> None:
    """Merges all the bad regions into 1 file"""
    df_bad = pd.read_csv("data/centromeres.csv")
    df_bad = df_bad[["chrom", "chromStart", "chromEnd"]]
    # telomeres, scaffold, contig
    df_bad = pd.concat(
        [df_bad, pd.read_csv("data/gaps_full.csv")[["chrom", "chromStart", "chromEnd"]]]
    )
    # ucsc hg38 -> mapping and sequencing -> problematic regions -> ENCODE Blicklist V2
    df = pd.read_csv("data/encode_blacklist_v2.csv")
    df = df.rename(columns={'#"chrom"': "chrom"})
    df_bad = pd.concat([df_bad, df[["chrom", "chromStart", "chromEnd"]]])
    df_bad["chromStart"] = df_bad["chromStart"].astype(int)
    df_bad["chromEnd"] = df_bad["chromEnd"].astype(int)
    df_bad = df_bad.drop_duplicates()
    df_bad = df_bad.sort_values(["chrom", "chromStart", "chromEnd"])
    df_bad.to_csv("data/all_excluded_regions.csv")


def get_intersected_rows(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """Finds the intersection between coordinates set in df1 and df2

    Args:
        df1 (pd.DataFrame): First Dataframe with columns "chr", "start", "end"
        df2 (pd.DataFrame): Second DataFrame with columns "chr", "start", "end"

    Returns:
        pd.DataFrame: DataFrame of intersected rows
    """
    df1 = df1.reset_index(drop=True, inplace=False).reset_index(
        drop=False, inplace=False
    )
    df2 = (
        df2.reset_index(drop=True, inplace=False)
        .reset_index(drop=False, inplace=False)
        .rename(
            columns={
                "chr": "chr_1",
                "start": "start_1",
                "end": "end_1",
                "index": "index_1",
            }
        )
    )
    chr_list = df1["chr"].unique()

    conn = sqlite3.connect(":memory:")
    df1.to_sql("df1", conn, index=False)
    df2.to_sql("df2", conn, index=False)

    all_intersections = []

    for chr_n in chr_list:
        # start_1 between start and end
        query = f"""
            select * 
            from df1 join df2 on start_1 between start and end
            where df1.chr = {chr_n} and df2.chr_1 = {chr_n}
            """
        df_intersect = pd.read_sql_query(query, conn)[["index", "index_1"]]
        all_intersections.append(df_intersect)

        # start between start_1 and end_1
        query = f"""
            select * 
            from df1 join df2 on start between start_1 and end_1
            where df1.chr = {chr_n} and df2.chr_1 = {chr_n}
            """
        df_intersect = pd.read_sql_query(query, conn)[["index", "index_1"]]
        all_intersections.append(df_intersect)

        # end_1 between start and end
        query = f"""
            select * 
            from df1 join df2 on end_1 between start and end
            where df1.chr = {chr_n} and df2.chr_1 = {chr_n}
            """
        df_intersect = pd.read_sql_query(query, conn)[["index", "index_1"]]
        all_intersections.append(df_intersect)

    df_intersected = pd.concat(all_intersections)
    df_intersected = df_intersected.drop_duplicates()
    df_intersected = pd.merge(df_intersected, df1, on="index", how="inner")
    df_intersected = pd.merge(df_intersected, df2, on="index_1", how="inner")
    df_intersected.drop(["index", "index_1"], axis=1, inplace=True)
    print("Number of rows in first DF:", df1.shape[0])
    print("Number of rows in second DF:", df2.shape[0])
    print("Number of intersected rows:", df_intersected.shape[0])
    return df_intersected


def filter_bad_regions(
    breakpoints_path: str, bad_regions_path: str, out_path: str
) -> None:
    """Removes Y chromosome and bad regions from breakpoints data

    Args:
        breakpoints_path (str): path to breakpoints
        bad_regions_path (str): path to bad regions data
        out_path (str): path to output data with excluded regions
    """
    df_bkpt = (
        pd.read_csv(breakpoints_path)[["hg38_chr", "hg38_coord"]]
        .rename(columns={"hg38_chr": "chr", "hg38_coord": "start"})
        .assign(end=lambda x: x['start'] + 1)
    )
    df_bkpt['chr'] = df_bkpt['chr'].astype(str)

    # remove Y chromosome
    print(df_bkpt["chr"].unique())
    df_bkpt = df_bkpt[df_bkpt["chr"] != "Y"]
    
    # remove bad regions
    df_bad_regions = pd.read_csv(bad_regions_path).rename(
        columns={"chrom": "chr", "chromStart": "start", "chromEnd": "end"}
    )
    df_bad_regions['chr'] = df_bad_regions['chr'].map(lambda x: x.replace("chr", ""))
    df_bad_regions['chr'] = df_bad_regions['chr'].astype(str)
    df_intersected = get_intersected_rows(df_bkpt, df_bad_regions)
    df_intersected.to_csv("data/tmp.csv")
    df_bkpt_all = pd.merge(df_bkpt, df_intersected, how="left")
    df_bkpt_all = df_bkpt_all[df_bkpt_all["chr_1"].isnull()]
    df_bkpt_all[["chr", "start", "end"]].to_csv(out_path)
    # 468 472
    print("Number of resulting breakpoints:", df_bkpt_all.shape[0])


if __name__ == "__main__":
    # get_all_bad_regions_list()
    filter_bad_regions(
        breakpoints_path="data/hg38_breakpoints_wo_err.csv",
        bad_regions_path="data/all_excluded_regions.csv",
        out_path="data/breakpoints_wo_bad_regions.csv",
    )
