# Human Protein Atlas (pathology)
import os.path
import pandas as pd
import zipfile
import os.path
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class HPAParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'HPA'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        url = self.config['hpa_pathology_url']
        disease_mapping = self.get_mapping_from_ontology(ontology="Disease",
                                                         source=None)
        protein_mapping = self.get_multiple_mapping_for_entity("Protein")
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        compressed_fileName = os.path.join(directory, url.split('/')[-1])
        file_name = '.'.join(url.split('/')[-1].split('.')[0:2])
        relationships_headers = self.config['relationships_headers']

        if self.download:
            self.download_db(url, directory)

        with zipfile.ZipFile(compressed_fileName) as z:
            if file_name == "pathology.tsv":
                pathology = self.parse_pathology_file(
                    z, file_name, protein_mapping, disease_mapping)

        # self.remove_directory(directory)

        return (pathology, relationships_headers)

    def parse_pathology_file(self, fhandler, file_name, protein_mapping, disease_mapping):
        url = self.config['linkout_url']
        pathology = defaultdict(set)
        first = True
        with fhandler.open(file_name) as f:
            df = pd.read_csv(f, sep='\t', header=None,
                             low_memory=False)
            df = df.fillna(0)
            first = True
            for index, row in df.iterrows():
                if first:
                    first = False
                    continue
                identifier = row[0]
                name = row[1]
                disease_name = row[2]
                hexpr = row[3]
                mexpr = row[4]
                lexpr = row[5]
                ndetected = row[6]
                if isinstance(row[7], str):
                    row[7] = float(row[7].replace('e', 'E'))
                if isinstance(row[8], str):
                    row[8] = float(row[8].replace('e', 'E'))
                if isinstance(row[9], str):
                    row[9] = float(row[9].replace('e', 'E'))
                if isinstance(row[10], str):
                    row[10] = float(row[10].replace('e', 'E'))

                uprog_pos = row[8]
                prog_pos = row[7] if row[7] != 0 else uprog_pos
                uprog_neg = row[10]
                prog_neg = row[9] if row[9] != 0 else uprog_neg

                linkout = url.replace("GENECODE", identifier).replace(
                    "GENENAME", name).replace("DISEASE", disease_name.replace(' ', '+'))

                if identifier in protein_mapping or name in protein_mapping:
                    if identifier in protein_mapping:
                        identifiers = protein_mapping[identifier]
                    else:
                        identifiers = protein_mapping[name]
                    for identifier in identifiers:
                        if disease_name in disease_mapping:
                            disease_id = disease_mapping[disease_name]
                            pathology[("protein", "detected_in_pathology_sample")].add((identifier, disease_id, "DETECTED_IN_PATHOLOGY_SAMPLE",
                                                                                        hexpr, mexpr, lexpr, ndetected, prog_pos, prog_neg, linkout, "Human Protein Atlas pathology"))

        return pathology

    def build_stats(self):
        stats = set()
        relationships, headers = self.parse()
        for entity, relationship in relationships:
            hpa_outputfile = os.path.join(self.import_directory,
                                          self.database_name.lower()+"_"+entity.lower()+"_"+relationship.lower()+".tsv")
            self.write_relationships(relationships[(entity, relationship)],
                                     headers[relationship], hpa_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(self.database_name,
                                                                              relationship, len(relationships[(entity, relationship)])))
            stats.add(self._build_stats(len(relationships[(entity, relationship)]),
                                        "relationships", relationship, self.database_name, hpa_outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
