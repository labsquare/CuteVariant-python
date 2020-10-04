import os
import csv

import cutevariant.commons as cm

LOGGER = cm.logger()


class PedReader:
    """PED *.fam file (PLINK sample information file) parser

    Data has the same structure of a fam file object
    https://www.cog-genomics.org/plink/1.9/formats#fam

    Format description:
        A tabular text file with no header line, and one line per sample with the
        following six fields:

        - Family ID (str): ('FID'); ex: "Valse"
        - Within-family ID (str): ('IID'; cannot be '0'); ex: "Emmanuel"
        - Within-family ID of father (str): ('0' if father isn't in dataset)
        - Within-family ID of mother (str): ('0' if mother isn't in dataset)
        - Sex code (int): ('1' = male, '2' = female, '0' = unknown)
        - Phenotype value (int): ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            We expect '0' for this last one.

    Recall, DB structure:
        - id INTEGER PRIMARY KEY ASC,
        - name TEXT,
        - fam TEXT DEFAULT 'fam',
        - father_id INTEGER DEFAULT 0,
        - mother_id INTEGER DEFAULT 0,
        - sex INTEGER DEFAULT 0,
        - phenotype INTEGER DEFAULT 0
    """

    def __init__(self, filepath, samples, *args, **kwargs):
        """

        Args:
            filepath: PED filepath
            samples: Database current samples
        """
        assert os.path.isfile(filepath)

        self.filepath = filepath
        self.samples = samples

    def __iter__(self):
        """Generator on PED file

        PED file is opened as a tabulated file.
        """
        with open(self.filepath, "r") as stream:
            reader = csv.reader(
                stream,
                delimiter="\t",
            )

            yield from self.get_samples(reader)

    def get_samples(self, reader):
        """Yield samples from PED file

        Returns:
            (generator): Generator of samples
        """
        samples_mapping = {
            (sample["fam"], sample["name"]): sample["id"] for sample in self.samples
        }

        for index, line in enumerate(reader, 1):
            if len(line) < 6:
                LOGGER.error(
                    "PED file conformity line <%s>; too few fields; expected at least 6",
                    index,
                )
                continue

            # Extract and validate data
            family_id = line[0]
            individual_id = line[1]
            father_id = line[2]
            mother_id = line[3]
            sex = int(line[4]) if line[4].isdigit() else 0
            phenotype = int(line[5]) if line[5].isdigit() else 0

            if sex not in (0, 1, 2):
                LOGGER.error(
                    "PED file conformity line <%s>; sex code <%s> not expected",
                    index,
                    sex,
                )
                continue

            if phenotype not in (0, 1, 2):
                LOGGER.error(
                    "PED file conformity line <%s>; phenotype code <%s> not expected",
                    index,
                    sex,
                )
                continue

            if individual_id == "0":
                LOGGER.error(
                    "PED file conformity line <%s>; Within-family ID/Sample name can't be '0'",
                    index,
                )
                continue

            # Test presence of sample in DB
            if (family_id, individual_id) not in samples_mapping.keys():
                LOGGER.error(
                    "PED file conformity line <%s>; sample <fam:%s>, <individual_id:%s> not found in database",
                    index, family_id, individual_id,
                )
                continue

            # Test presence of parents in DB
            father_key = (family_id, father_id)
            mother_key = (family_id, mother_id)
            if {father_key, mother_key} - samples_mapping.keys():
                # If set not empty: 1 tuple is not found in DB
                # => will be replaced by 0 (unknown)
                LOGGER.warning(
                    "PED file conformity line <%s>; parents <%s> or <%s> not found in database",
                    index, father_key, mother_key,
                )

            new_sample = {
                "id": samples_mapping[(family_id, individual_id)],  # Get DB sample id
                "fam": family_id,
                "father_id": samples_mapping.get(father_key, 0),
                "mother_id": samples_mapping.get(mother_key, 0),
                "sex": sex,
                "phenotype": phenotype,
            }

            # print(new_sample)
            yield new_sample
