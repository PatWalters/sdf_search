#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import PandasTools
from tqdm.auto import tqdm
import csv
from contextlib import ExitStack
import os
import subprocess
import duckdb
import time
from functools import wraps
from FPSim2 import FPSim2Engine
import pandas as pd
import click
from typing import List

def timed(func):
    """
    Decorator to measure and store the runtime of a function in seconds.

    The runtime is stored in the `last_runtime` attribute of the decorated function.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        wrapper.last_runtime = elapsed
        return result
    wrapper.last_runtime = 0.0
    return wrapper

@timed
def process_sdf(sdf_name: str, base_name: str) -> None:
    """
    Converts an SDF file to CSV and SMILES formats, extracting all properties.

    Args:
        sdf_name: Path to the input SDF file.
        base_name: Base name for output files (CSV and SMILES).
    """
    csv_name = f"{base_name}.csv"
    smi_name = f"{base_name}.smi"

    try:
        suppl = Chem.SDMolSupplier(sdf_name)
        with ExitStack() as stack:
            csv_fs = stack.enter_context(open(csv_name, "w", newline=""))
            smi_fs = stack.enter_context(open(smi_name, "w"))
            writer = csv.writer(csv_fs)
            prop_cols = []
            for idx, mol in enumerate(tqdm(suppl)):
                if mol is None:
                    continue
                smiles = Chem.MolToSmiles(mol)
                prop_dict = mol.GetPropsAsDict()
                prop_dict["Title"] = mol.GetProp("_Name")
                if idx == 0:
                    prop_cols = [k for k in prop_dict if k.upper() not in ["SMILES"]]
                    out_cols = ["SMILES", "Index"] + prop_cols
                    writer.writerow(out_cols)
                row = [smiles, idx] + [prop_dict.get(col, "") for col in prop_cols]
                writer.writerow(row)
                print(f"{smiles} {idx}", file=smi_fs)
    except Exception as e:
        print(f"Error processing {sdf_name}: {e}")

@timed
def build_fpsim2_database(smi_name: str, h5_name: str, db_type: str) -> int:
    """
    Builds an FpSim2 database from a SMILES file.

    Args:
        smi_name: Path to the input SMILES file.
        h5_name: Path to the output HDF5 database file.
        db_type: Type of fingerprint ('sim' or 'sub').

    Returns:
        0 if successful, 1 otherwise.
    """
    param_dict = {
        "sim": ['Morgan', '{"radius": 2, "fpSize": 256}'],
        "sub": ['RDKitPattern','{"fpSize": 2048, "tautomerFingerprints": true}']
    }
    if db_type not in param_dict:
        print(f"Invalid db_type: {db_type}. Must be 'sim' or 'sub'.")
        return 1

    fp_type, fp_params = param_dict[db_type]

    try:
        if os.path.isfile(h5_name):
            os.remove(h5_name)
        subprocess.run(
            ['fpsim2-create-db', "--fp_type", fp_type, "--fp_params", fp_params, smi_name, h5_name],
            capture_output=True,
            text=True,
            check=True
        )
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Command failed. Stderr: {e.stderr}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1

@timed
def build_duckdb_database(db_name: str, csv_name: str) -> None:
    """
    Builds a DuckDB database from a CSV file.

    Args:
        db_name: Path to the output DuckDB database file.
        csv_name: Path to the input CSV file.
    """
    if os.path.isfile(db_name):
        os.remove(db_name)
    try:
        with duckdb.connect(database=db_name) as conn:
            conn.execute(f"""
                CREATE TABLE sdf_data AS
                SELECT *
                FROM read_csv_auto('{csv_name}');
            """)
    except Exception as e:
        print(f"Error building DuckDB database: {e}")


def build_sdf_database(sdf_name: str, outfile_prefix: str) -> None:
    """
    Orchestrates the conversion of an SDF to FPSim2, and DuckDB databases,
    printing the runtime for each step.

    Args:
        sdf_name: Path to the input SDF file.
        outfile_prefix: Prefix for all output files.
    """
    smi_name = f"{outfile_prefix}.smi"
    h5_sim_name = f"{outfile_prefix}_sim.h5"
    h5_sub_name = f"{outfile_prefix}_sub.h5"
    ddb_name = f"{outfile_prefix}.ddb"
    csv_name = f"{outfile_prefix}.csv"

    try:
        print(f"Reading {sdf_name} to create {csv_name} and {smi_name}")
        process_sdf(sdf_name, outfile_prefix)
        print(f"Done. Runtime: {process_sdf.last_runtime:.2f} seconds")

        print(f"Building FpSim2 substructure database {h5_sub_name} from {smi_name}")
        build_fpsim2_database(smi_name, h5_sub_name, "sub")
        print(f"Done. Runtime: {build_fpsim2_database.last_runtime:.2f} seconds")

        print(f"Building FpSim2 similarity database {h5_sim_name} from {smi_name}")
        build_fpsim2_database(smi_name, h5_sim_name, "sim")
        print(f"Done. Runtime: {build_fpsim2_database.last_runtime:.2f} seconds")

        print(f"Building DuckDB database {ddb_name} from {csv_name}")
        build_duckdb_database(ddb_name, csv_name)
        print(f"Done. Runtime: {build_duckdb_database.last_runtime:.2f} seconds")

        print(f"Removing temporary files {smi_name} and {csv_name}")
        os.unlink(smi_name)
        os.unlink(csv_name)
    except Exception as e:
        print(f"Error in build_sdf_database: {e}")

@timed
def similarity_search(
    h5_name: str,
    ddb_name: str,
    query_smiles: str,
    threshold: float = 0.35,
    limit_size: int = 10000
) -> pd.DataFrame | None:
    """
    Performs a similarity search using FpSim2 and retrieves matching records from DuckDB.

    Args:
        h5_name: Path to the FpSim2 HDF5 database.
        ddb_name: Path to the DuckDB database.
        query_smiles: Query molecule in SMILES format.
        threshold: Similarity threshold for hits.
        limit_size: Maximum number of hits to return.

    Returns:
        DataFrame of matching records, or None if no hits found.
    """
    try:
        engine = FPSim2Engine(h5_name)
        sim_results = engine.similarity(query_smiles, threshold=threshold)
        if len(sim_results) == 0:
            print("No similar molecules found.")
            return None

        sim_df = pd.DataFrame(sim_results)
        if len(sim_df) > limit_size:
            print(f"Your search returned more than {limit_size} hits. Truncating the hit list to {limit_size}.")
            sim_df = sim_df.head(limit_size)

        idx_list = sim_df['mol_id'].tolist()
        if not idx_list:
            print("No valid indices found in similarity results.")
            return None

        placeholders = ','.join(['?'] * len(idx_list))
        query = f'SELECT * FROM sdf_data WHERE "Index" IN ({placeholders})'

        with duckdb.connect(database=ddb_name) as conn:
            res_df = conn.execute(query, idx_list).df()
        # Add the similarity scores to the result DataFrame
        res_df = res_df.merge(sim_df, left_on='Index', right_on='mol_id', how='left')
        res_df.rename(columns={'coeff': 'Tanimoto'}, inplace=True)
        return res_df
    except Exception as e:
        print(f"Error during similarity search: {e}")
        return


def smarts_in_smiles(smiles_list: List[str], smarts: str) -> List[bool]:
    """
    Returns a list of booleans indicating if the SMARTS is contained in each SMILES.

    Args:
        smiles_list: List of SMILES strings.
        smarts: SMARTS pattern to search for.

    Returns:
        List of booleans.
    """
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        raise ValueError("Invalid SMARTS pattern")

    results = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results.append(False)
        else:
            results.append(mol.HasSubstructMatch(query))
    return results

@timed
def substructure_search(
    h5_name: str,
    ddb_name: str,
    query_smiles: str,
    limit_size: int = 10000
) -> pd.DataFrame | None:
    """
    Performs a substructure search using FpSim2 and retrieves matching records from DuckDB.

    Args:
        h5_name: Path to the FpSim2 HDF5 database.
        ddb_name: Path to the DuckDB database.
        query_smiles: Query molecule in SMILES format.
        limit_size: Maximum number of hits to return.

    Returns:
        DataFrame of matching records, or None if no hits found.
    """
    try:
        engine = FPSim2Engine(h5_name)
        hit_idx = engine.substructure(query_smiles)
        if len(hit_idx) == 0:
            print("No substructure matches found.")
            return None

        if len(hit_idx) > limit_size:
            print(f"Your search returned more than {limit_size} hits. Truncating the hit list to {limit_size}.")
            hit_idx = hit_idx[:limit_size]
        hit_idx = hit_idx.tolist()
        placeholders = ','.join(['?'] * len(hit_idx))
        query = f'SELECT * FROM sdf_data WHERE "Index" IN ({placeholders})'

        with duckdb.connect(database=ddb_name) as conn:
            res_df = conn.execute(query, hit_idx).df()
        # Filter results to ensure substructure match (in case of false positives)
        res_df['match'] = smarts_in_smiles(res_df['SMILES'].tolist(), query_smiles)
        hit_df = res_df[res_df['match']]
        print(f"Filtered {len(res_df) - len(hit_df)} false positives; {len(hit_df)} true matches remain.")
        return hit_df.drop(columns=['match'])
    except Exception as e:
        print(f"Error during substructure search: {e}")
        return None


def process_search(
    search_type: str,
    prefix: str,
    query_smiles: str,
    outfile_name: str = None,
    threshold: float = 0.35,
    limit_size: int = 10000
) -> None:
    """
    Processes a similarity or substructure search and writes results to an output file.

    Args:
        search_type: 'similarity' or 'substructure'.
        prefix: Prefix for the input database files.
        query_smiles: Query molecule in SMILES format.
        outfile_name: Output file name (.sdf or .csv). If None, returns a DataFrame.
        threshold: Similarity threshold (for similarity search only).
        limit_size: Maximum number of hits to return.
    """
    sub_h5_name = f"{prefix}_sub.h5"
    sim_h5_name = f"{prefix}_sim.h5"
    ddb_name = f"{prefix}.ddb"

    search_res = None
    try:
        if search_type == "sim":
            search_res = similarity_search(sim_h5_name, ddb_name, query_smiles, limit_size=limit_size, threshold=threshold)
            search_res.sort_values(by=['Tanimoto'], ascending=False, inplace=True)
            print(f"Runtime for similarity search: {similarity_search.last_runtime:.2f} seconds")
        elif search_type == "sub":
            search_res = substructure_search(sub_h5_name, ddb_name, query_smiles, limit_size=limit_size)
            print(f"Runtime for substructure search: {substructure_search.last_runtime:.2f} seconds")
        else:
            print(f"Unknown search type: {search_type}")
            return

        if outfile_name is None:
            print(f"Search returned {len(search_res)} results.")
            return search_res

        if search_res is not None and not search_res.empty:
            search_res = search_res.fillna("")
            print(f"Found {len(search_res)} hits.")
            if outfile_name.lower().endswith(".sdf"):
                PandasTools.AddMoleculeColumnToFrame(search_res, smilesCol="SMILES")
                if "Title" in search_res.columns:
                    [x.SetProp("_Name",str(y)) for x,y in search_res[["ROMol","Title"]].values]
                PandasTools.WriteSDF(search_res, out=outfile_name, properties=list(search_res.columns))
            elif outfile_name.lower().endswith(".csv"):
                search_res.to_csv(outfile_name, index=False)
            else:
                print(f"Unsupported output file format: {outfile_name}")
                return
            print(f"Search results written to {outfile_name}")
        else:
            print("No hits found.")
    except Exception as e:
        print(f"Error in process_search: {e}")
    return None


@click.group()
def cli():
    """sdf_search command line interface."""
    pass

@cli.command("build")
@click.argument("sdf_name", type=click.Path(exists=True))
@click.argument("outfile_prefix", type=click.Path())
def build_db_cmd(sdf_name, outfile_prefix):
    """Build FpSim2 and DuckDB databases from an SDF."""
    build_sdf_database(sdf_name, outfile_prefix)

@cli.command("search")
@click.argument("search_type", type=click.Choice(["sim", "sub"]))
@click.argument("prefix", type=click.Path())
@click.argument("query_smiles")
@click.argument("outfile_name", type=click.Path())
@click.option("--limit", default=10000, show_default=True, help="Maximum number of hits to return.")
@click.option("--threshold", default=0.35, show_default=True, help="Similarity threshold (for similarity search only).")
def search_cmd(search_type, prefix, query_smiles, outfile_name, limit, threshold):
    """Run a similarity or substructure search"""
    process_search(search_type, prefix, query_smiles, outfile_name, limit_size=limit, threshold=threshold)

if __name__ == "__main__":
    cli()
