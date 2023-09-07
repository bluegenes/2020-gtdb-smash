import os
import sys
import glob
import argparse
import pandas as pd


def summarize_gather_at(rank, tax_assign, gather_results):
    # from https://github.com/dib-lab/sourmash/issues/1011
    # collect!

    sum_uniq_weighted = defaultdict(float)
    for row in gather_results:
        match_ident = row['name']
        match_ident = get_ident(match_ident)
        lineage = tax_assign[match_ident]
        lineage = pop_to_rank(lineage, rank)
        assert lineage[-1].rank == rank, lineage[-1]

        f_uniq_weighted = row['f_unique_weighted']
        f_uniq_weighted = float(f_uniq_weighted)
        sum_uniq_weighted[lineage] += f_uniq_weighted

    items = list(sum_uniq_weighted.items())
    items.sort(key = lambda x: -x[1])

    for k, v in items:
        print(rank, f'{v:.3f}', sourmash.lca.display_lineage(k))













def find_filename(row, glob_dir, glob_col, full_paths=False):
    identifier = row[glob_col]
    filepattern = os.path.join(glob_dir, f"*{identifier}*")
    found_files = glob.glob(filepattern)
    if len(found_files) >= 1:
        if full_paths:
            row["filename"] = os.path.abspath(found_files[0])
        else:
            row["filename"] = os.path.basename(found_files[0])
        if len(found_files) > 1:
            fileinfo = "\n".join(found_files)
            sys.stderr.write(f"Warning: {identifier} found more than one match in directory {glob_dir}. " \
                             f"Only first match returned. All matches are:\n {fileinfo}\n")
    elif len(found_files) == 0:
        row["filename"] = ""
        sys.stderr.write(f"Warning: {identifier} found no matches in directory {glob_dir}.\n")
    return row


def main(args):
    #read in csv
    queryDF= try_reading_csv(args.query_csv)
    # glob for filenames
    queryDF = queryDF.apply(find_filename, axis=1, args=(args.sigfile_directory,args.identifier_column, args.full_paths))
    # write full csv
    queryDF.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("query_csv")
    p.add_argument("--output-csv")
    p.add_argument("--sigfile-directory", default=os.getcwd(), help="directory with sig or fasta files that will be used to build sbt")
    p.add_argument("--identifier-column", default="accession", help="column in query_csv that can be used to identify corresponding filename")
    p.add_argument("--full-paths", action="store_true", help="return full filepath instead of basename")
    args = p.parse_args()
    if not args.output_csv:
        args.output_csv=(args.query_csv).rsplit(".", 1)[0] + ".with_filenames.csv"
    sys.exit(main(args))
