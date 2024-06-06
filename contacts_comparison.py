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
import sys

from Bio import AlignIO
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


def match_aligned_positions_with_original(aln_path):
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

    :param aln_path: the path to the alignment (fasta format)
    :type aln_path: str
    :return: the match between the original position and the aligned position.
    :rtype: dict
    """
    logging.info("get the corresponding positions between the sequences of the alignment and the original ones.")
    data = {}
    aln = AlignIO.read(aln_path, "fasta")
    for record in aln:
        data[record.id] = {}
        aln_seq = record.seq
        for index in range(len(aln_seq)):
            gaps = aln_seq[:index].count('-')
            original_index = len(aln_seq[:index]) - gaps
            data[record.id][original_index + 1] = {"aligned position": index + 1, "residue": aln_seq[index]}
    return data


def get_contact_file_paths_by_condition(path):
    """
    By condition, get the sample names and the contact file paths.

    :param path: the CSV file path of the directory paths by condition.
    :type path: str
    :return: the samples and paths by condition.
    :rtype: dict
    """
    logging.info("get the paths of the contacts files by positions:")
    data = {}
    conditions_paths = pd.read_csv(path, sep=",", header=0)
    pattern_fn = re.compile("^outliers_(.+)\\.csv")
    for _, row in conditions_paths.iterrows():
        data[row["condition"]] = {}
        for fn in os.listdir(row["path"]):
            match = pattern_fn.search(fn)
            if match:
                sample = match.group(1)
                data[row["condition"]][sample]= os.path.join(row["path"], fn)
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
    :param nb_smp_cond_1: the number of residues in contacts for the  first condition.
    :type nb_smp_cond_1: int
    :param output_dir: the path to the output directory.
    :type output_dir: str
    :param roi: the region of interest making contacts.
    :type roi: str
    :return: the contact in the first condition that are not in the second condition.
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
    out_path = os.path.join(output_dir, f"{roi}_{condition_1}_not_in_{condition_2}.csv")
    df.to_csv(out_path, sep=",", index=False)
    logging.info(f"\t\t- {len(df)} contact{'s' if len(df) > 1 else ''} in {condition_1} not in {condition_2}: "
                 f"{out_path}")
    return df


def get_commons(data, condition_1, condition_2, nb_smp_cond_1, nb_smp_cond_2, output_dir, roi):
    common_keys = sorted(list(set(list(data[condition_1].keys())).intersection(list(data[condition_2].keys()))))
    commons_dict = {"position alignment": [], "number of contacts": [], "domain": [],
                    f"number of {condition_1} samples with contacts": [], f"number of samples {condition_1}": [],
                    f"number of {condition_2} samples with contacts": [], f"number of samples {condition_2}": [],
                    f"original positions {condition_1}": [], f"original positions {condition_2}": []}
    for aln_diff_pos in common_keys:
        commons_dict["position alignment"].append(aln_diff_pos)
        commons_dict["number of contacts"].append(data[condition_1][aln_diff_pos]["count"])
        commons_dict["domain"].append(data[condition_1][aln_diff_pos]["domain"])
        commons_dict[f"number of {condition_1} samples with contacts"].append(
            len(data[condition_1][aln_diff_pos]["original"]))
        commons_dict[f"number of samples {condition_1}"].append(nb_smp_cond_1)
        original_contacts_cond1 = None
        for smp in data[condition_1][aln_diff_pos]["original"]:
            if original_contacts_cond1:
                original_contacts_cond1 = (f"{original_contacts_cond1} | "
                                           f"{data[condition_1][aln_diff_pos]['original'][smp]['position']}:{smp}")
            else:
                original_contacts = f"{data[condition_1][aln_diff_pos]['original'][smp]['position']}:{smp}"
        commons_dict[f"original positions {condition_1}"].append(original_contacts_cond1)
        commons_dict[f"number of {condition_2} samples with contacts"].append(
            len(data[condition_2][aln_diff_pos]["original"]))
        commons_dict[f"number of samples {condition_2}"].append(nb_smp_cond_2)
        original_contacts_cond2 = None
        for smp in data[condition_2][aln_diff_pos]["original"]:
            if original_contacts_cond2:
                original_contacts = (f"{original_contacts_cond2} | "
                                     f"{data[condition_2][aln_diff_pos]['original'][smp]['position']}:{smp}")
            else:
                original_contacts = f"{data[condition_2][aln_diff_pos]['original'][smp]['position']}:{smp}"
        commons_dict[f"original positions {condition_2}"].append(original_contacts_cond2)
    df = pd.DataFrame.from_dict(commons_dict)
    out_path = os.path.join(output_dir, f"{roi}_{condition_1}_in_{condition_2}.csv")
    df.to_csv(out_path, sep=",", index=False)
    logging.info(f"\t\t- {len(df)} contact{'s' if len(df) > 1 else ''} in {condition_1} and in {condition_2}: "
                 f"{out_path}")
    return df


def compare_contacts_by_condition(data_whole_contacts, out_dir, data_files_by_condition, region_of_interest):
    """
    Search for contacts different and common positions between the conditions.
    :param data_whole_contacts: the whole contact positions.
    :type data_whole_contacts: dict
    :param out_dir: the path to the output directory.
    :type out_dir: str
    :param data_files_by_condition: the samples and paths by condition.
    :type data_files_by_condition: dict
    :return: the different and common contacts positions.
    :rtype: dict
    """
    logging.info("get the differences and common positions between the conditions:")
    out_data = {}
    for idx_1 in range(len(data_whole_contacts.keys()) - 1):
        cond_1 = list(data_whole_contacts.keys())[idx_1]
        nb_smp_condition_1 = len(data_files_by_condition[cond_1])
        for idx_2 in range(idx_1 + 1, len(data_whole_contacts.keys())):
            cond_2 = list(data_whole_contacts.keys())[idx_2]
            comparison = f"{cond_1} vs. {cond_2}"
            logging.info(f"\t{comparison}:")
            df_differences = get_differences(data_whole_contacts, cond_1, cond_2, len(data_files_by_condition[cond_1]),
                                             out_dir, region_of_interest)
            df_commons = get_commons(data_whole_contacts, cond_1, cond_2, nb_smp_condition_1,
                                     len(data_files_by_condition[cond_2]), out_dir, region_of_interest)
            out_data[comparison] = {"differences": df_differences, "commons": df_commons}
    return out_data




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
    parser.add_argument("-t", "--md-time", required=True, type=int,
                        help="the molecular dynamics duration in nanoseconds.")
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

    positions_original_alignment = match_aligned_positions_with_original(args.aln)
    conditions_samples_paths = get_contact_file_paths_by_condition(args.input)
    whole_contacts_by_condition = get_whole_contact_positions(conditions_samples_paths, positions_original_alignment)

    # for k, v in whole_contacts_by_condition.items():
    #     print(f"{k}:")
    #     sorted_v = dict(sorted(v.items()))
    #     for k1, v1 in sorted_v.items():
    #         print(f"\t{k1}: {v1}")

    differences_commons = compare_contacts_by_condition(whole_contacts_by_condition, args.out,
                                                        conditions_samples_paths, args.roi)