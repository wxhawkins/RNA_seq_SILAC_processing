import argparse
import re
import requests
import time
from concurrent import futures
from pathlib import Path

import pandas as pd

# Establish argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", action="store", type=str, dest="in_path", default=Path(Path.cwd().parent / "in_files" / "rna_seq_input.xlsx"), help="path to input file containing RNA seq data.")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_file_name", default="id_dump.csv", help="Name of output file to be generated.")
parser.add_argument("-e", "--error", action="store", type=str, dest="err_file_name", default="error_log.csv", help="Name of error log file to be generated.")
parser.add_argument("-s", "--sample", action="store", type=str, dest="sample", default="None", help="Take random subsample of DataFrame.")
parser.add_argument("-t", "--threads", action="store", type=int, dest="num_theads", default=1, help="Run queries accross multiple threads.")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False, help="Print ID and error information")
parser.add_argument("-f", "--force", action="store_true", dest="force", default=False, help="Suppress non-critical input requests.")

args = parser.parse_args()

ERR_PATH = Path(Path.cwd() / "out_files" / args.err_file_name)
OUT_PATH = Path(Path.cwd() / "out_files" / args.out_file_name)



class Protein:
    def __init__(self, rna_id_="None", prot_id_="None", gene_name_="None"):
        self.rna_id = rna_id_
        self.prot_id = prot_id_
        self.gene_name = gene_name_
    def __str__(self):
        return "hello?" 


def extract_rna_id(str_):
    """
        Extract RNA ID (GenBank accession number) from Nr Description information.
    """

    regex = re.compile("(.+)//")
    hit = re.search(regex, str_)
    if hit is not None:
        return hit.group(1)

    # Return "N/A" if no accession number found
    return None


def get_cleaned_df(in_file_path, delim="\t"):
    """
        Perform minor data cleaning on starting DataFrame in preparation for downstream processes.
    """
    try:
        df = pd.read_excel(in_file_path, sep=delim)
    except PermissionError:
        print("Could open input file. Ensure file is not open.")
        return False

    print(df.shape)

    # Fill empty cells with placeholder
    df.fillna("N/A", inplace=True)

    # Delete rows with no nr_description
    filt = df["Nr Description"] == "N/A"
    df.drop(df[filt].index, inplace=True)

    # Take sample of input data if requested
    if args.sample != "None":
        df = df.sample(min(df.shape[0], int(args.sample)))

    # Add accession column with info extracted form "nr_description"
    df.loc[:, "accession"] = df["Nr Description"].apply(extract_rna_id)
    print(df.shape)
    return df


def dump_failures(err_accs):
    """
        Write accessions for which there was no hit in Uniprot to a CSV file.
    """

    with open(ERR_PATH, "w") as err_out_file:
        err_out_file.write("\n".join(err_accs))


def dump_data(data):
    """
        Write RNA IDs, Uniprot IDs and gene names to Excel file.
    """

    with open(OUT_PATH, "w") as out_file:
        for row in data:
            out_file.write(",".join(row) + "\n")


def make_request(rna_id):
    """
        Queries Uniprot database with Genbank accession number and extracts Uniprot accession number and gene name
    """

    url = "http://www.uniprot.org/uniprot/"
    payload = {
                "format": "tab",
                "query": (rna_id),
                "columns": "id,genes",
              }

    try:
        r = requests.get(url, params=payload)

        # Parse data to extract Uniprot ID and gene name
        lines = r.text.split("\n")
        info = lines[1].split("\t")
        u_id = info[0]
        gene_name = info[1].split(" ")[0]

        # Set gene name to "None" if none found
        if gene_name.strip() == "":
            gene_name = "None"
    except IndexError:
        # If an error is encountered, display error message and return filler values for Uniprot ID and gene name
        if args.verbose:
            print(f"\n-------------ERROR ENCOUNTERED AT ACC {rna_id}--------------\n")
        
        u_id = "None"
        gene_name = "None"

    return u_id, gene_name


def get_info(acc):

    """
        Submit Genbank accession numbers for querying against Uniprot database.
    """
    
    u_id, gene_name = make_request(acc)

    # If query fails, retry with version number suffix removed
    if "None" in (acc, u_id, gene_name):
        re_1 = re.compile(r"(.*)\.")
        mod_acc = re.search(re_1, acc).group(1)
        u_id, gene_name = make_request(mod_acc)

    return (acc, u_id, gene_name)


def handle_paths():
    """
        Check with user before overriding existing files with destination path.
    """

    # Skip override checking if requested
    if args.force:
        return True

    try:
        if OUT_PATH.exists():
            response = input(f"{args.out_file_name} already exists. Delete this file? (y/n)\t").lower()
            if response not in ("y", "yes"):
                return False

        if ERR_PATH.exists():
            response = input(f"{args.err_file_name} already exists. Delete this file? (y/n)\t").lower()
            if response not in ("y", "yes"):
                return False
    except PermissionError:
        print("Could not overriede existing output file. Ensure output files are not open.")
        return False

    return True


def spawn_threads(num_workers):
    """
        Spawn specified number of threads and execute Uniprot ID queries.
    """

    data_list = []
    with futures.ThreadPoolExecutor(max_workers=len(acc_pool)) as executor:
        new_rows = (list(executor.map(get_info, acc_pool)))
        data_list += [row for row in new_rows if "None" not in row]
    
        if args.verbose:
            print("\n".join([",  ".join(row) for row in new_rows]))

    return data_list


if __name__ == "__main__":
    # Abort if user does not approve file override
    if handle_paths():
        start = time.time()

        # Read in and clean starting RNA seq dataset
        df = get_cleaned_df(args.in_path)
        accs = df["accession"]
        acc_pool  = []
        data_list = []
        err_accs  = []

        # Submit query for every accession number
        for index, acc in enumerate(accs):
            # Add accessions to pool until number equals number of threads
            if len(acc_pool) < args.num_theads:
                acc_pool.append(acc)
            else:
                # Have each thread submit Uniprot query for one accession from pool
                data_list += spawn_threads(len(acc_pool))

                # Empty pool
                acc_pool = []
                acc_pool.append(acc)
            
                # Print status of run
                print(f" {index} ".center(100, "-"))

        # Clear out queue at end of the run
        num_workers = len(acc_pool)
        data_list += spawn_threads(num_workers)
        
        # Dump information to CSV files
        dump_data(data_list)
        dump_failures(err_accs)

        print(f"Finished in {round((time.time() - start), 1)} seconds")
        print(f"Information successfully collected for {len(data_list)} RNA IDs of {df.shape[0]}")
