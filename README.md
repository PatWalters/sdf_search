# sdf_search

## Installation

```
git clone git@github.com:PatWalters/sdf_search.git  
cd sdf_search  
pip install .
```

Then put the script sdf_search.py somewhere in your PATH.

## TLDR

```
sdf_search.py build infile.sdf db_name

sdf_search.py search [sub|sim] db_name query output.[csv|sdf] [--limit N] [--threshold T]
```

## Introduction

The script `sdf_search.py` is a straightforward command-line Python program designed for performing substructure and similarity searches on chemical structures and data extracted from SD files. I often encounter situations where I need to search an SDF for structures containing a specific substructure or similar to a compound of interest. While there are various graphical interfaces available to load structures into, I prefer working from the [command line](https://www.amazon.com/Beginning-Was-Command-Line-ebook/dp/B0011GA08E/). I wanted a simple utility that could convert an SDF into a database and then perform substructure or similarity searches. The script operates in two modes, which are specified by the first command line argument.

- `build` - transform an SDF into a searchable databaase  
- `search` - perform a **sub**structure or **sim**ilarity search on the database created in the `build` phase

## Building a database

Before performing a search, we need to build a database. To do this we use the `sdf_search.py build` command, like this.

`sdf_search.py build infile.sdf db_name`

where infile.sdf is the input SDF and db_name is the base name of the database.

**Please note: If the database already exists, it will be overwritten without warning.**

As an example, let's assume we want to create a database from the 11 million compounds in the Aldrich Mark Select collection, which is contained in the SDF `MarketSelectMain.sdf`.

`sdf_search.py build MarketSelectMain.sdf ams`

This process will create two files.

- `ams.h5`  an HDF5 file containing a molecule index and the fingerprints for the molecules.  
- `ams.ddb`  a DuckDB database containing all the data fields in the SDF.

## Searching a database

To perform a substructure or similarity search, we use the `simsearch.py search` command. The syntax for searching is  
   
`sdf_search.py search [sub|sim] db_name query output`  
   
In a hopefully intuitive way, you can specify `'sub'` to perform a substructure search or `'sim'` to perform a similarity search. The `db_name` parameter is the same as the one used during the `build` phase. For example, if we created `ams.h5` and `ams.bdb`, the db_name would be `'ams'`. The argument `'query'` is a SMILES string for a similarity search or a SMARTS string for a substructure search. To avoid errors, it's best to always quote this argument. The final argument, `'output',` specifies the output file name. The output file type, which can be `csv` or `sdf`, is determined by the file extension.

**Please note: If the output file already exists, it will be overwritten without warning.**

The `search` command also supports two optional flags.

- --limit - sets an upper limiit on the number of results returned (default 10000\)  
- --threshold - sets a lower threshold for similarity searches.

## Examples

### Perform a substructure search for molecules containing 3-aminopyridine

`sdf_search.py search ams '[NH2]c1cccnc1' out.csv`

This command perform a substructure search on database with more than 11 million records in less than a second. The output file type is taken from the file extension of the output file.

- csv - creates an output csv file  
- sdf - creates an output SDF

If we wanted the command above to output an SDF, we could simply change the file extension on the output file to `.sdf`

`sdf_search.py search ams '[NH2]c1cccnc1' out.sdf`

### Perform a similarity search for molecules similar to 3-aminopyridine.

`sdf_search.py search sim  ams  'Nc1cccnc1'  out.csv`

By default, the Tanimoto similarity search threshold is set to 0.35. This threshold can be adjusted using the `--threshold` flag. If we wish to change this threshold to 0.5, we can add the flag.

`sdf_search.py search sim  ams  'Nc1cccnc1'  out.csv --threshold 0.5`

### Search limits

By default, searches return a maximum of 10000 results. If this limit is reached, you will see a message like this.

`Your search returned more than 10000 hits. Truncating the hit list to 10000.`

You can use the `--limit` flag to change this limit.

sdf_search.py search sub mydb  'c1ccccc1' out.csv --limit 20000

## Caveats

This script assumes that you have a properly formatted SDF with consistent data fields for each record. I can't guarantee the results if this is not the case.

## Acknowledgements

This script is a simple wrapper around the [FPSim2](https://github.com/chembl/FPSim2) and [DuckDB](https://duckdb.org/) libraries, both of which are excellent open-source projects created by people much smarter than I am. Please consider supporting these projects if you find this script helpful.  
