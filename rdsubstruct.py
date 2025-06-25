from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary
import pickle
from tqdm.auto import tqdm
from sdf_search import timed

@timed
def create_and_save_substruct_library(smi_file: str, lib_file: str) -> None:
    """
    Reads a SMILES file and creates a SubstructLibrary, then saves it to disk.
    """
    mols = []
    names = []
    with open(smi_file) as f:
        num_lines = sum(1 for _ in f)
    with open(smi_file) as f:
        for line in tqdm(f,total=num_lines):
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            smi, name = parts[0], parts[1]
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mols.append(mol)
                names.append(name)
    lib = rdSubstructLibrary.SubstructLibrary()
    for mol in mols:
        lib.AddMol(mol)
    # Save both the library and names list
    with open(lib_file, "wb") as out:
        pickle.dump((lib, names), out)

@timed
def search_substruct_library(lib_file: str, smarts: str):
    """
    Loads a SubstructLibrary from disk and searches with a SMARTS query.
    Returns a list of (index, name) for matches.
    """
    with open(lib_file, "rb") as f:
        lib, names = pickle.load(f)
    query = Chem.MolFromSmarts(smarts)
    if not query:
        raise ValueError("Invalid SMARTS pattern")
    hits = lib.GetMatches(query, maxResults=10000)
    return [(idx, names[idx]) for idx in hits]

if __name__ == "__main__":
    # Example usage
    smi_file = "mydb.smi"  # Input SMILES file
    lib_file = "mydb.pkl"  # Output library file
    create_and_save_substruct_library(smi_file, lib_file)
    print("Done")
    print(create_and_save_substruct_library.last_runtime)

    smarts_query = "c1ccccc1N"  # Example SMARTS query (benzene ring)
    results = search_substruct_library(lib_file, smarts_query)
    print(search_substruct_library.last_runtime)

    print(len(results))
#    for idx, name in results:
#        print(f"Match found: Index {idx}, Name {name}")