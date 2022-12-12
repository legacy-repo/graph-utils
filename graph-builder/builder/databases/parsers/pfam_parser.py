# Pfam database
import os.path
import pandas as pd
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class PfamParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'Pfam'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        entity_header = self.config['entity_header']
        relationship_headers = self.config['relationship_headers']

        directory = os.path.join(self.database_directory, 'Pfam')
        self.check_directory(directory)
        protein_mapping = self.get_mapping_for_entity(entity="Protein")
        valid_proteins = list(set(protein_mapping.values()))

        ftp_url = self.config['ftp_url']
        filename = self.config['full_uniprot_file']
        # url = config['test']

        if self.download:
            self.download_db(ftp_url+filename, directory)

        stats = set()
        if os.path.exists(os.path.join(directory, filename)):
            fhandler = self.read_gzipped_file(
                os.path.join(directory, filename))
            identifier = None
            description = []
            lines = []
            missed = 0
            entities = set()
            relationships = defaultdict(set)
            is_first = True
            i = 0
            read_lines = 0
            num_entities = 0
            num_relationships = {}
            try:
                for line in fhandler:
                    i += 1
                    read_lines += 1
                    if line.startswith("# STOCKHOLM"):
                        if identifier is not None:
                            entities.add(
                                (identifier, 'Functional_region', name, " ".join(description), "PFam"))
                            if len(entities) == 100:
                                outputfile = os.path.join(
                                    self.import_directory, 'Functional_region.tsv')
                                self.print_files(entities, entity_header,
                                                 outputfile=outputfile, is_first=is_first)
                                num_entities += len(entities)
                                if 'mentioned_in_publication' in relationships:
                                    outputfile = os.path.join(self.import_directory,
                                                              'Functional_region_mentioned_in_publication.tsv')
                                    self.print_files(relationships['mentioned_in_publication'],
                                                     relationship_headers['mentioned_in_publication'],
                                                     outputfile=outputfile, is_first=is_first)
                                    if 'mentioned_in_publication' not in num_relationships:
                                        num_relationships['mentioned_in_publication'] = 0
                                    num_relationships['mentioned_in_publication'] += len(
                                        relationships['mentioned_in_publication'])
                                if 'found_in_protein' in relationships:
                                    outputfile = os.path.join(self.import_directory,
                                                              'Functional_region_found_in_protein.tsv')
                                    self.print_files(relationships['found_in_protein'], relationship_headers['found_in_protein'],
                                                     outputfile=outputfile, is_first=is_first, filter_for=('END_ID', valid_proteins))
                                    if 'found_in_protein' not in num_relationships:
                                        num_relationships['found_in_protein'] = 0
                                    num_relationships['found_in_protein'] += len(
                                        relationships['found_in_protein'])
                                entities = set()
                                relationships = defaultdict(set)
                                is_first = False
                            identifier = None
                            description = []
                    elif line.startswith("#=GF"):
                        data = line.rstrip('\r\n').split()
                        if 'AC' in data:
                            identifier = data[2].split('.')[0]
                        elif 'DE' in data:
                            name = " ".join(data[2:])
                        elif 'RM' in data:
                            relationships['mentioned_in_publication'].add(
                                (identifier, data[2], "MENTIONED_IN_PUBLICATION", "PFam"))
                        elif 'CC' in data:
                            description.append(" ".join(data[2:]))
                    elif not line.startswith('//'):
                        data = line.rstrip('\r\n').split()
                        protein, positions = data[0].split('/')
                        protein = protein.replace('.', '-')
                        start, end = positions.split('-')
                        sequence = data[1]
                        relationships['found_in_protein'].add(
                            (identifier, protein, "FOUND_IN_PROTEIN", start, end, sequence, "PFam"))
                        if protein.split('-')[0] != protein:
                            relationships['found_in_protein'].add((identifier, protein.split(
                                '-')[0], "FOUND_IN_PROTEIN", start, end, sequence, "PFam"))
            except UnicodeDecodeError:
                lines.append(i)
                missed += 1

            fhandler.close()

            if len(entities) > 0:
                outputfile = os.path.join(
                    self.import_directory, 'Functional_region.tsv')
                self.print_files(entities, entity_header,
                                 outputfile=outputfile, is_first=is_first)
                num_entities += len(entities)
                outputfile = os.path.join(self.import_directory,
                                          'Functional_region_mentioned_in_publication.tsv')
                self.print_files(relationships['mentioned_in_publication'],
                                 relationship_headers['mentioned_in_publication'],
                                 outputfile=outputfile, is_first=is_first)
                num_relationships['mentioned_in_publication'] += len(
                    relationships['mentioned_in_publication'])
                outputfile = os.path.join(
                    self.import_directory, 'Functional_region_found_in_protein.tsv')
                self.print_files(relationships['found_in_protein'],
                                 relationship_headers['found_in_protein'], outputfile=outputfile, is_first=is_first)
                num_relationships['found_in_protein'] += len(
                    relationships['found_in_protein'])

            stats.add(self._build_stats(num_entities, "entity",
                      "Functional_region", "Pfam", 'Functional_region.tsv', self.updated_on))

            for rel in num_relationships:
                stats.add(self._build_stats(num_relationships[rel], "relationship", rel.upper(),
                                            "Pfam", 'Functional_region_'+rel+'.tsv', self.updated_on))

        # self.remove_directory(directory)

        return stats

    def print_files(self, data, header, outputfile, is_first, filter_for=None):
        df = pd.DataFrame(list(data), columns=header)
        if filter_for is not None:
            df = df[df[filter_for[0]].isin(filter_for[1])]
        if not df.empty:
            with open(outputfile, 'a', encoding='utf-8') as f:
                df.to_csv(path_or_buf=f, sep='\t',
                          header=is_first, index=False, quotechar='"',
                          line_terminator='\n', escapechar='\\')

    def build_stats(self):
        stats = self.parse()
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
