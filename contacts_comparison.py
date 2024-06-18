#!/usr/bin/env python3

"""
Created on 05 Jun. 2024
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
import logging
import os
import re
import shutil
import sys

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from pymsaviz import MsaViz, get_msa_testdata
import pandas as pd


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def match_aligned_positions_with_original(aln):
    """
    From an alignment, match the position alignment with the position of the residue in the original sequence.
    A dictionary is created with the following structure:
    <SEQ ID>:
        <POSITION 1 IN ORIGINAL SEQUENCE>:
            aligned position: <POSITION>
            residue: <RESIDUE>
        <POSITION 2 IN ORIGINAL SEQUENCE>:
            aligned position: <POSITION>
            residue: <RESIDUE>
        ...

    :param aln: the alignment.
    :type aln: Bio.Align.MultipleSeqAlignment
    :return: the match between the original position and the aligned position.
    :rtype: dict
    """
    logging.info("get the corresponding positions between the sequences of the alignment and the original ones.")
    data = {}
    for record in aln:
        data[record.id] = {}
        aln_seq = record.seq
        for index in range(len(aln_seq)):
            gaps = aln_seq[:index].count('-')
            original_index = len(aln_seq[:index]) - gaps
            data[record.id][original_index + 1] = {"aligned position": index + 1, "residue": aln_seq[index]}
    return data


def get_domains(domains_file_path, aln, data_whole_positions):
    """
    Extract the domains and update the positions with the ones of the alignment.

    :param domains_file_path: the path to the protein domains CSV file.
    :type domains_file_path: str
    :param aln: the alignment.
    :type aln: Bio.Align.MultipleSeqAlignment
    :param data_whole_positions: the match between the original position and the aligned position.
    :type data_whole_positions: dict
    :return: the domains' coordinates.
    :rtype: panddas.DataFrame
    """
    logging.info(f"extracting the domains positions from the domains file:")
    data = {"domain": [], "start": [], "end": []}

    # getting the length of the original sequence from the sequence in the alignment
    seq_length = None
    pattern = re.compile("(.+)_domains.csv")
    match = pattern.search(os.path.basename(domains_file_path))
    if match:
        smp = match.group(1)
    else:
        logging.error(f"\tno match in the file name \"{os.path.basename(domains_file_path)}\" and the pattern "
                      f"\"{pattern.pattern}\".")
        sys.exit(1)
    for record in aln:
        if record.id == smp:
            seq_length = len(record.seq.replace("-", ""))
            break
    if not seq_length:
        logging.error(f"\tthe sequence \"{smp}\" from the domains file was not found in the sequences of the "
                      f"alignment.")
        sys.exit(1)

    df = pd.read_csv(domains_file_path, sep=",")
    previous_position = 0
    previous_domain = None
    idx = None
    row = None
    for idx, row in df.iterrows():
        if idx == 0 and row["start"] != 1:
            # the first domain is not starting at position 1
            # record the before domain
            data["domain"].append(f"before {row['domain']}")
            data["start"].append(1)
            data["end"].append(row["start"] - 1)
        elif row["start"] - 1 != previous_position:
            # the domains are not just one after the other and the domain
            # record between the domain
            data["domain"].append(f"between_{previous_domain}_and_{row['domain']}")
            data["start"].append(previous_position + 1)
            data["end"].append(row["start"] - 1)
        # record the domain
        data["domain"].append(row["domain"])
        data["start"].append(row["start"])
        data["end"].append(row["end"])
        previous_position = row["end"]
        previous_domain = row["domain"]

    if idx == len(df) - 1 and row["end"] < seq_length:
        # after the last domain
        data["domain"].append(f"after_{row['domain']}")
        data["start"].append(previous_position + 1)
        data["end"].append(seq_length)

    df_domains = pd.DataFrame.from_dict(data)
    for idx, row in df_domains.iterrows():
        df_domains.at[idx, "start"] = data_whole_positions[smp][row["start"]]["aligned position"]
        df_domains.at[idx, "end"] = data_whole_positions[smp][row["end"]]["aligned position"]
    logging.info("\tdomains updated.")

    return df_domains


def get_contact_file_paths_by_condition(path, grouped):
    """
    By condition, get the sample names and the contact file paths.

    :param path: the CSV file path of the directory paths by condition.
    :type path: str
    :param grouped: the grouped conditions.
    :type grouped: list
    :return: the samples and paths by condition.
    :rtype: dict
    """
    logging.info("get the paths of the contacts files by positions:")
    data = {}
    conditions_paths = pd.read_csv(path, sep=",", header=0)
    pattern_fn = re.compile("^outliers_(.+)\\.csv")
    for _, row in conditions_paths.iterrows():
        condition = row["condition"]
        if grouped and condition in grouped:
            condition = "-".join(grouped)
        if not condition in data:
            data[condition] = {}
        for fn in os.listdir(row["path"]):
            match = pattern_fn.search(fn)
            if match:
                sample = match.group(1)
                data[condition][sample]= os.path.join(row["path"], fn)
    for condition in data:
        logging.info(f"\t{condition}: {len(data[condition])} files")
    return data


def get_whole_contact_positions(cond_smp_paths, data_whole_positions):
    """
    Extract the whole contacts for each condition and get the matching positions between the alignment and the original
    sequence.
    The output is a dictionary:
    {<POSITION IN THE ALIGNMENT>:
        {domain: <DOMAIN>,
         count: <NB OF ATOMS IN CONTACT>,
         original:
            {<SAMPLE>:
                position: <POSITION IN THE ORIGINAL SEQUENCE>},
             <SAMPLE>:
                position: <POSITION IN THE ORIGINAL SEQUENCE>},
              ...},
    <POSITION IN THE ALIGNMENT>:
        ...}

    :param cond_smp_paths: the samples and directories' paths of the contacts by condition.
    :type cond_smp_paths: dict
    :param data_whole_positions: the match between the original position and the aligned position.
    :type data_whole_positions: dict
    :return: the whole contact positions.
    :rtype: dict
    """
    logging.info("get the contacts positions by condition:")
    data = {}
    for condition, smp_and_paths in cond_smp_paths.items():
        data[condition] = {}
        for smp, contact_path in smp_and_paths.items():
            df = pd.read_csv(contact_path, sep=",", header=0)
            df = df.sort_values(by=["second partner position"])
            for _, row in df.iterrows():
                position_in_aln = data_whole_positions[smp][row["second partner position"]]["aligned position"]
                if position_in_aln not in data[condition]:
                    data[condition][position_in_aln] = {"domain": row["second partner domain"],
                                                        "count": row["number atoms contacts"],
                                                        "original": {smp: {"position": row["second partner position"]}}}
                else:
                    if smp not in data[condition][position_in_aln]["original"]:
                        data[condition][position_in_aln]["count"] += row["number atoms contacts"]
                        data[condition][position_in_aln]["original"][smp] = {"position": row["second partner position"]}
        logging.info(f"\t{condition}: {len(data[condition])} contacts positions.")
    return data


def get_differences(data, condition_1, condition_2, nb_smp_cond_1, output_dir, roi):
    """
    Get the contacts positions differences between two conditions.

    :param data: the whole contact positions.
    :type data: dict
    :param condition_1: the first condition.
    :type condition_1: str
    :param condition_2: the second condition.
    :type condition_2: str
    :param nb_smp_cond_1: the number of residues in contacts for the first condition.
    :type nb_smp_cond_1: int
    :param output_dir: the path to the output directory.
    :type output_dir: str
    :param roi: the region making contacts with the others.
    :type roi: str
    :return: the contacts in the first condition that are not in the second condition.
    :rtype: pandas.DataFrame
    """
    difference_keys = sorted(list(set(list(data[condition_1].keys())) - set(list(data[condition_2].keys()))))
    differences_dict = {"position alignment": [], "number of contacts": [], "domain": [],
                        "number of samples with contacts": [], "number of samples": [], "original positions": []}
    for aln_diff_pos in difference_keys:
        differences_dict["position alignment"].append(aln_diff_pos)
        differences_dict["number of contacts"].append(data[condition_1][aln_diff_pos]["count"])
        differences_dict["domain"].append(data[condition_1][aln_diff_pos]["domain"])
        differences_dict["number of samples with contacts"].append(len(data[condition_1][aln_diff_pos]["original"]))
        differences_dict["number of samples"].append(nb_smp_cond_1)
        original_contacts = None
        for smp in data[condition_1][aln_diff_pos]["original"]:
            if original_contacts:
                original_contacts = (f"{original_contacts} | "
                                     f"{data[condition_1][aln_diff_pos]['original'][smp]['position']}:{smp}")
            else:
                original_contacts = f"{data[condition_1][aln_diff_pos]['original'][smp]['position']}:{smp}"
        differences_dict["original positions"].append(original_contacts)
    df = pd.DataFrame.from_dict(differences_dict)
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"different_contacts_for_{roi}_{condition_1}_vs_{condition_2}.csv")
    df.to_csv(out_path, sep=",", index=False)
    logging.info(f"\t\t- {len(df)} contact{'s' if len(df) > 1 else ''} in {condition_1} not in {condition_2}: "
                 f"{out_path}")


def get_commons(data, condition_1, condition_2, nb_smp_cond_1, nb_smp_cond_2, output_dir, roi):
    """
    Get the common contacts positions between two conditions.

    :param data: the whole contact positions.
    :type data: dict
    :param condition_1: the first condition.
    :type condition_1: str
    :param condition_2: the second condition.
    :type condition_2: str
    :param nb_smp_cond_1: the number of residues in contacts for the first condition.
    :type nb_smp_cond_1: int
    :param nb_smp_cond_2: the number of residues in contacts for the second condition.
    :type nb_smp_cond_2: int
    :param output_dir: the path to the output directory.
    :type output_dir: str
    :param roi: the region making contacts with the others.
    :type roi: str
    :return: the contacts in the first condition that are in common with the second condition.
    :rtype: pandas.DataFrame
    """
    common_keys = sorted(list(set(list(data[condition_1].keys())).intersection(list(data[condition_2].keys()))))
    commons_dict = {"position alignment": [], "number of contacts": [], "domain": [],
                    f"number of {condition_1} samples with contacts": [], f"number of samples {condition_1}": [],
                    f"number of {condition_2} samples with contacts": [], f"number of samples {condition_2}": [],
                    f"original positions {condition_1}": [], f"original positions {condition_2}": []}
    for aln_comm_pos in common_keys:
        commons_dict["position alignment"].append(aln_comm_pos)
        commons_dict["number of contacts"].append(data[condition_1][aln_comm_pos]["count"])
        commons_dict["domain"].append(data[condition_1][aln_comm_pos]["domain"])
        commons_dict[f"number of {condition_1} samples with contacts"].append(
            len(data[condition_1][aln_comm_pos]["original"]))
        commons_dict[f"number of samples {condition_1}"].append(nb_smp_cond_1)
        original_contacts_cond1 = None
        for smp in data[condition_1][aln_comm_pos]["original"]:
            if original_contacts_cond1:
                original_contacts_cond1 = (f"{original_contacts_cond1} | "
                                           f"{data[condition_1][aln_comm_pos]['original'][smp]['position']}:{smp}")
            else:
                original_contacts_cond1 = f"{data[condition_1][aln_comm_pos]['original'][smp]['position']}:{smp}"
        commons_dict[f"original positions {condition_1}"].append(original_contacts_cond1)
        commons_dict[f"number of {condition_2} samples with contacts"].append(
            len(data[condition_2][aln_comm_pos]["original"]))
        commons_dict[f"number of samples {condition_2}"].append(nb_smp_cond_2)
        original_contacts_cond2 = None
        for smp in data[condition_2][aln_comm_pos]["original"]:
            if original_contacts_cond2:
                original_contacts_cond2 = (f"{original_contacts_cond2} | "
                                     f"{data[condition_2][aln_comm_pos]['original'][smp]['position']}:{smp}")
            else:
                original_contacts_cond2 = f"{data[condition_2][aln_comm_pos]['original'][smp]['position']}:{smp}"
        commons_dict[f"original positions {condition_2}"].append(original_contacts_cond2)
    df = pd.DataFrame.from_dict(commons_dict)
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"common_contacts_for_{roi}_{condition_1}_vs_{condition_2}.csv")
    df.to_csv(out_path, sep=",", index=False)
    logging.info(f"\t\t- {len(df)} contact{'s' if len(df) > 1 else ''} in {condition_1} and in {condition_2}: "
                 f"{out_path}")


def compare_contacts_by_condition(data_whole_contacts, out_dir, data_files_by_condition, region_of_interest):
    """
    Search for contacts different and common positions between the conditions.
    :param data_whole_contacts: the whole contact positions.
    :type data_whole_contacts: dict
    :param out_dir: the path to the output directory.
    :type out_dir: str
    :param data_files_by_condition: the samples and paths by condition.
    :type data_files_by_condition: dict
    :param region_of_interest: the region making contacts with the others.
    :type region_of_interest: str
    :return: the different and common contacts positions.
    :rtype: dict
    """
    logging.info("get the different and common positions between the conditions:")
    for idx_1 in range(len(data_whole_contacts.keys()) - 1):
        cond_1 = list(data_whole_contacts.keys())[idx_1]
        nb_smp_condition_1 = len(data_files_by_condition[cond_1])
        for idx_2 in range(idx_1 + 1, len(data_whole_contacts.keys())):
            cond_2 = list(data_whole_contacts.keys())[idx_2]
            versus = f"{cond_1} vs. {cond_2}"
            logging.info(f"\t{versus}:")
            get_differences(data_whole_contacts, cond_1, cond_2, len(data_files_by_condition[cond_1]),
                            os.path.join(out_dir, f"{cond_1}_vs_{cond_2}"), region_of_interest)
            get_commons(data_whole_contacts, cond_1, cond_2, nb_smp_condition_1, len(data_files_by_condition[cond_2]),
                        os.path.join(out_dir, f"{cond_1}_vs_{cond_2}"), region_of_interest)


def subsample_msa(whole_msa, dict_samples, conditions_ids, tmp_directory):
    """
    Sample the MSA based on the samples present in the studied conditions.

    :param whole_msa: the MSA with all the samples.
    :type whole_msa: Bio.Align.MultipleSeqAlignment
    :param dict_samples: the dictionary with samples' names by condition.
    :type dict_samples: dict
    :param conditions_ids: the two conditions that are compared.
    :param tmp_directory: str
    :return: the path to the sub-sampled alignment.
    :rtype: str
    """
    conditions = conditions_ids.split("_vs_")
    samples_to_keep = []
    for condition in conditions:
        samples_to_keep = samples_to_keep + list(dict_samples[condition].keys())
    new_msa = MultipleSeqAlignment([])
    for record in whole_msa:
        if record.id in samples_to_keep:
            new_msa.append(record)
    path_tmp = os.path.join(tmp_directory, f"aln_{'_'.join(conditions)}.fa")
    AlignIO.write(new_msa, path_tmp, "fasta")
    return path_tmp


def plot_msa(aln, region_of_interest, data_samples, domains, out_dir):
    """
    For two conditions, plot the Multiple Sequences Alignment with annotations of number of contacts.

    :param aln: the MSA.
    :type aln: Bio.Align.MultipleSeqAlignment
    :param region_of_interest: the region making contacts with the others.
    :type region_of_interest: str
    :param data_samples: the samples by condition.
    :type data_samples: dict
    :param domains: the domains' information.
    :type domains: pandas.core.frame.DataFrame
    :param out_dir: the path to the output directory.
    :type out_dir: str
    """
    logging.info("plotting the alignments with the contacts annotations:")
    tmp_dir = os.path.join(out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    for versus in os.listdir(out_dir):
        if os.path.isdir(os.path.join(out_dir, versus)) and versus != "tmp":
            logging.info(f"\t- {versus.replace('_vs_', ' vs. ')}:")
            versus_msa_file = subsample_msa(aln, data_samples, versus, tmp_dir)
            for _, row in domains.iterrows():
                msa_viz = MsaViz(versus_msa_file, start=row["start"], end=row["end"], wrap_length=60,
                                 show_consensus=False)
                for comparison in ["common", "different"]:
                    if comparison == "common":
                        annotations_color = "blue"
                    else:
                        annotations_color = "red"
                    df = pd.read_csv(os.path.join(out_dir, versus,
                                                  f"{comparison}_contacts_for_{region_of_interest}_{versus}.csv"))
                    df_by_domain = df.loc[df["domain"] == row["domain"]]
                    for _, row_in_domain in df_by_domain.iterrows():
                        msa_viz.add_markers([row_in_domain["position alignment"]], color=annotations_color,
                                            marker=f"${row_in_domain['number of contacts']}$")
                out = os.path.join(out_dir, versus, f"msa_plot_{versus}_{row['domain'].replace(' ', '_')}.svg")
                msa_viz.savefig(out)
                logging.info(f"\t\t- {row['domain']} alignment plot saved: {out}")
    shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    Retrieve the common contacts between Molecular Dynamics simulations of various systems. The input files are the ones
    from the plot_contacts.py script (https://github.com/njeanne/plot_contacts).
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-r", "--roi", required=True, type=str,
                        help="the region of interest making contacts with other domains.")
    parser.add_argument("-a", "--aln", required=True, type=str,
                        help="the path to the alignment file (fasta format).")
    parser.add_argument("-d", "--domains", required=True, type=str,
                        help="a sample CSV domains annotation file for one of the sequences of the alignment. The file "
                             "name must be <SAMPLE>_domain.csv with <SAMPLE> matching exactly the sequence name in the "
                             "alignment. The file must contain 3 comma separated columns: domain, start and end. The "
                             "start and end position for each domain must be 1 indexed.")
    parser.add_argument("-t", "--md-time", required=True, type=int,
                        help="the molecular dynamics duration in nanoseconds.")
    parser.add_argument("-g", "--group", required=False, nargs="+", type=str,
                        help="a list of conditions, separated by spaces, to group as they appear in the first column "
                             "of the input file.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("input", type=str,
                        help="the path to the CSV (comma separated without header) file which first column is the "
                             "condition, the second column the path of the directory containing the plots_contacts "
                             "script CSV output files.")
    args = parser.parse_args()

    # create output directory if necessary
    os.makedirs(args.out, exist_ok=True)
    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(args.out, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")

    msa = AlignIO.read(args.aln, "fasta")

    positions_original_alignment = match_aligned_positions_with_original(msa)
    updated_domains = get_domains(args.domains, msa, positions_original_alignment)
    conditions_samples_paths = get_contact_file_paths_by_condition(args.input, args.group)

    #todo: remove
    for k, v in conditions_samples_paths.items():
        print(f"{k}:")
        for k1, v1 in v.items():
            print(f'\t{k1}: {v1}')

    whole_contacts_by_condition = get_whole_contact_positions(conditions_samples_paths, positions_original_alignment)
    compare_contacts_by_condition(whole_contacts_by_condition, args.out, conditions_samples_paths, args.roi)
    plot_msa(msa, args.roi, conditions_samples_paths, updated_domains, args.out)
