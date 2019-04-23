from .abstractreader import AbstractReader
import vcf


class FakeReader(AbstractReader):
    def __init__(self):
        super().__init__(None)

    def _get_variants(self):
        yield {
            "chr": "11",
            "pos": 125010,
            "ref": "T",
            "alt": "A",
            "annotations": [
                {"transcript": "NM_234234", "gene": "CFTR"},
                {"transcript": "NM_234235", "gene": "CFTR"},
            ],
            "samples": [{"name": "sacha", "gt": 1}],
        }

        yield {
            "chr": "12",
            "pos": 125010,
            "ref": "T",
            "alt": "A",
            "annotations": [
                {"transcript": "NM_234234", "gene": "CFTR"},
                {"transcript": "NM_234235", "gene": "CFTR"},
            ],
            "samples": [{"name": "sacha", "gt": 1}],
        }

        yield {
            "chr": "13",
            "pos": 125010,
            "ref": "T",
            "alt": "A",
            "annotations": [
                {"transcript": "NM_234234", "gene": "CFTR"},
                {"transcript": "NM_234235", "gene": "CFTR"},
            ],
            "samples": [{"name": "sacha", "gt": 1}],
        }

    def _get_fields(self):
        """ Extract fields informations from VCF fields """

        yield {
            "name": "chr",
            "category": "variants",
            "description": "chromosom",
            "type": "str",
        }
        yield {
            "name": "pos",
            "category": "variants",
            "description": "position",
            "type": "int",
        }

        yield {
            "name": "ref",
            "category": "variants",
            "description": "reference base",
            "type": "str",
        }
        yield {
            "name": "alt",
            "category": "variants",
            "description": "alternative base",
            "type": "str",
        }

        yield {
            "name": "gt",
            "category": "samples",
            "description": "genotype",
            "type": "int",
        }

        yield {
            "name": "af",
            "category": "samples",
            "description": "allele frequency",
            "type": "float",
        }

        yield {
            "name": "gene",
            "category": "annotations",
            "description": "gene name",
            "type": "str",
        }

        yield {
            "name": "transcript",
            "category": "annotations",
            "description": "gene transcripts",
            "type": "str",
        }

    def _get_samples(self):
        return ["sacha"]
