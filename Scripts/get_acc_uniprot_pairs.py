import argparse
import re
import requests
import time
from concurrent import futures
from pathlib import Path

import pandas as pd

# Establish argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", action="store", type=str, dest="in_path", default=Path(Path.cwd().parent / "in_files" / "bulk_rna_seq.tsv"), help="path to input file containing RNA seq data.")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_file_name", default="id_dump.csv", help="Name of output file to be generated.")
parser.add_argument("-e", "--error", action="store", type=str, dest="err_file_name", default="error_log.csv", help="Name of error log file to be generated.")
parser.add_argument("-s", "--sample", action="store", type=str, dest="sample", default="None", help="Take random subsample of DataFrame.")
parser.add_argument("-t", "--threads", action="store", type=int, dest="num_theads", default=1, help="Run queries accross multiple threads.")

args = parser.parse_args()

ERR_PATH = Path(Path.cwd() / "out_files" / args.err_file_name)
OUT_PATH = Path(Path.cwd() / "out_files" / args.out_file_name)


def clean_df(df):
    """
        Perform minor changes to starting DataFrame in preparation for downstream processes.
    """

    # Fill empty cells with placeholder
    df.fillna("N/A", inplace=True)

    # Delete rows with no nr_description
    filt = df["Nr Description"] == "N/A"
    df.drop(df[filt].index, inplace=True)

    # Take sample of input data if requested
    if args.sample != "None":
        df = df.sample(min(df.shape[0], int(args.sample)))

    # Add accession column with info extracted form "nr_description"
    df.loc[:, "accession"] = df["Nr Description"].apply(get_accession)

    return df


def get_accession(str_):
    """
        Extract GenBank accession number from Nr Description information.
    """

    re_1 = re.compile("(.+)//")
    re_hit_1 = re.search(re_1, str_)
    if re_hit_1 is not None:
        return re_hit_1.group(1)

    # Return "N/A" if no accession number found
    return "N/A"

def dump_errs(err_accs):
    """
        Write accessions for which there was no hit in Uniprot to a CSV file.
    """

    with open(ERR_PATH, "w") as err_out_file:
        err_out_file.write("\n".join(err_accs))

def dump_info(data):
    """
        Write accessions, Uniprot IDs and gene names to CSV file.
    """

    with open(OUT_PATH, "w") as out_file:
        for row in data:
            out_file.write(",".join(row) + "\n")


def get_info(acc):

    """
        Query Uniprot database for given accession number, returning matched accession number, Uniprot ID and gene name.
    """

    url = "http://www.uniprot.org/uniprot/"
    payload = {
                "format": "tab",
                "query": (acc),
                "columns": "id,genes",
              }

    try:
        r = requests.get(url, params=payload)
        lines = r.text.split("\n")
        info = lines[1].split("\t")
        u_id = info[0]
        gene_name = info[1].split(" ")[0]

        # Set gene name to "N/A" if none found
        if gene_name.strip() == "":
            gene_name = "N/A"
    except IndexError:
        # If an error is encountered, display error message and return filler values for Uniprot ID and gene name
        print(f"\n-------------ERROR ENCOUNTERED AT ACC {acc}--------------\n")
        
        u_id = "None"
        gene_name = "None"

    return (acc, u_id, gene_name)


def handle_paths():
    """
        Check with user before overriding existing files with destination path.
    """

    if OUT_PATH.exists():
        response = input(f"{args.out_file_name} already exists. Delete this file? (y/n)\t").lower()
        if response not in ("y", "yes"):
            return False

    if ERR_PATH.exists():
        response = input(f"{args.err_file_name} already exists. Delete this file? (y/n)\t").lower()
        if response not in ("y", "yes"):
            return False

    return True


def retry_accs(err_accs, data_list):
    """
        For accessions that had no hit to Uniprot database, try to query again without version number (.X) tail.
    """

    err_accs_final = []
    for acc in err_accs:
        # Try without extension
        re_1 = re.compile(r"(.*)\.")
        mod_acc = re.search(re_1, acc).group(1)
        print(f"mod_acc = {mod_acc}")
        data = get_info(mod_acc)

        # Append to data to data_list if hit found
        if "None" not in data:
            data_list += data
        else:
            err_accs_final.append(acc)

    return data_list, err_accs_final


if __name__ == "__main__":
    # Abort if user does not approve file override
    if handle_paths():
        start = time.time()

        # Read in and clean starting RNA seq dataset
        df = pd.read_csv(args.in_path, sep="\t")
        df = clean_df(df)

        accs = df["accession"]
        acc_pool = []
        data_list = []  # Will be list of tuples (accession, u_id, gene_name)
        err_accs = []

        # Submit query for every accession number
        for acc in accs:
            # Add accessions to pool until number equals number of threads
            if len(acc_pool) < args.num_theads:
                acc_pool.append(acc)
            else:
                # Have each thread submit Uniprot query for one accession from pool
                temp_dict = {}
                with futures.ThreadPoolExecutor(max_workers=args.num_theads) as executor:
                    new_rows = (list(executor.map(get_info, acc_pool)))
                    print("\n".join([",  ".join(row) for row in new_rows]))
                    data_list += [row for row in new_rows if "None" not in row]
                    err_accs += [row[0] for row in new_rows if "None" in row]

                # Empty pool
                acc_pool = []
                acc_pool.append(acc)

        # Clear out queue at end of the run
        for acc in acc_pool:
            new_row = (get_info(acc))
            if "None" not in new_row:
                print(",  ".join(new_row))
                data_list.append(new_row)
            else:
                err_accs.append(acc)
            
        # Retry accs without version number (.X) tails
        data_list, err_accs = retry_accs(err_accs, data_list)

        # Dump information to CSV files
        dump_info(data_list)
        dump_errs(err_accs)

        print(f"Finished in {round((time.time() - start), 1)} seconds")
