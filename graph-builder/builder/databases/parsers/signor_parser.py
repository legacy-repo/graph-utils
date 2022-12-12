import os.path
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class SIGNORParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'SIGNOR'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)

        url = self.config['url']
        modifications = self.config['modifications']
        amino_acids = self.config['amino_acids']
        accronyms = self.config['accronyms']
        entities_header = self.config['entities_header']
        relationships_headers = self.config['rel_headers']

        entities = set()
        relationships = defaultdict(set)

        filename = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)

        entities, relationships = self.parse_substrates(filename, modifications,
                                                        accronyms, amino_acids)

        return entities, relationships, entities_header, relationships_headers

    def parse_substrates(self, filename, modifications, accronyms, amino_acids):
        entities = set()
        relationships = defaultdict(set)
        first = True
        with open(filename, 'r', encoding="utf-8") as fhandler:
            for line in fhandler:
                if first:
                    first = False
                    continue

                data = line.rstrip("\r\n").split("\t")
                source = data[2]
                target = data[6]
                regulation = data[8]
                mechanism = data[9]
                residue_mod = data[10]
                seq_window = data[11]
                organism = data[12]
                pubmedid = data[21]
                if organism in self.organisms and mechanism in modifications and residue_mod != '':
                    if len(residue_mod) > 3:
                        residue = ''.join(residue_mod[0:3])
                        position = ''.join(residue_mod[3:])
                    if residue in amino_acids:
                        residue = amino_acids[residue]
                        modification = modifications[mechanism]
                        if mechanism in accronyms:
                            modified_protein_id = target+"_"+residue + \
                                position+"-"+accronyms[mechanism]
                            entities.add((modified_protein_id, "Modified_protein",
                                         target, seq_window, position, residue, "SIGNOR"))
                            relationships[("Protein", "has_modified_site")].add((target, modified_protein_id,
                                                                                 "HAS_MODIFIED_SITE", "SIGNOR"))
                            relationships[("Peptide", "has_modified_site")].add((seq_window.upper(), modified_protein_id,
                                                                                 "HAS_MODIFIED_SITE", "SIGNOR"))
                            relationships[("Modified_protein", "has_modification")].add((modified_protein_id, modification,
                                                                                         "HAS_MODIFICATION", "SIGNOR"))
                            relationships[('Substrate', 'is_substrate_of')].add((modified_protein_id, source,
                                                                                 "IS_SUBSTRATE_OF", regulation, "CURATED", 5, "SIGNOR"))
                            if pubmedid != '':
                                relationships['Modified_protein_Publication', 'mentioned_in_publication'].add((pubmedid, modified_protein_id,
                                                                                                               "MENTIONED_IN_PUBLICATION"))
        return entities, relationships

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_headers = self.parse()
        entity_outputfile = os.path.join(
            self.import_directory, "SIGNOR_Modified_protein.tsv")
        self.write_entities(entities, entities_header, entity_outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(
            self.database_name, "Modified_protein", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Modified_protein",
                                    self.database_name, entity_outputfile, self.updated_on))
        for entity, relationship in relationships:
            rel_header = ["START_ID", "END_ID", "TYPE", "source"]
            prefix = 'SIGNOR_' + entity.lower()
            if relationship in relationships_headers:
                rel_header = relationships_headers[relationship]
            if relationship == 'mentioned_in_publication':
                prefix = entity
            outputfile = os.path.join(
                self.import_directory, prefix+"_"+relationship.lower()+".tsv")
            self.write_relationships(relationships[(entity, relationship)],
                                     rel_header, outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(self.database_name,
                                                                              relationship, len(relationships[(entity, relationship)])))
            stats.add(self._build_stats(len(relationships[(entity, relationship)]),
                                        "relationships", relationship, self.database_name, outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
