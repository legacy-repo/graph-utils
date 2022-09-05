import os.path
import gzip
import logging
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser
from builder.databases.parsers.stitch_parser import STITCHParser

logger = logging.getLogger(__name__)


class NoSuchFiles(Exception):
    pass


class PhosphoSitePlusParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'PhosphoSitePlus'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def check_files(self, files):
        directory = os.path.join(self.database_directory, self.database_name)
        exists = True
        not_exist_files = []
        for file in files:
            file_name = os.path.join(directory, file)
            if not os.path.exists(file_name):
                exists = False
                not_exist_files.append(file)

        if not exists:
            logger.warn("Caution: registration necessary to download files")
            logger.warn("Please download the following files from PhosphoSitePlus site into %s manually." %
                        directory)
            print('- ' + '\n- '.join(not_exist_files) + '\n\n')
        return exists

    def parse(self):
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        modifications = self.config['modifications']
        annotation_files = self.config['annotation_files']
        entities_header = self.config['entities_header']
        relationships_headers = self.config['rel_headers']
        entities = set()
        relationships = defaultdict(set)
        if not self.check_files(list(self.config['annotation_files'].values()) + self.config['site_files']):
            raise NoSuchFiles("Please check the error log.")

        for site_file in self.config['site_files']:
            file_name = os.path.join(directory, site_file)
            with gzip.open(file_name, 'r') as f:
                sites, site_relationships = self.parse_sites(f, modifications)
                entities.update(sites)
                for r in site_relationships:
                    relationships[r].update(site_relationships[r])

        for er in annotation_files:
            entity, relationship_type = er.split('-')
            file_name = os.path.join(directory, annotation_files[er])
            with gzip.open(file_name, 'r') as f:
                if entity == "disease":
                    mapping = self.get_mapping_from_ontology(ontology="Disease",
                                                             source=None)
                    relationships[(entity, relationship_type)].update(
                        self.parse_disease_annotations(f, mapping))
                elif entity == "biological_process":
                    mapping = self.get_mapping_from_ontology(ontology="Gene_ontology",
                                                             source=None)
                    relationships[(entity, relationship_type)].update(
                        self.parse_regulation_annotations(f, mapping))
                elif entity == "substrate":
                    relationships[(entity, relationship_type)
                                  ] = self.parse_kinase_substrates(f)

        return entities, relationships, entities_header, relationships_headers

    def parse_sites(self, fhandler, modifications):
        entities = set()
        relationships = defaultdict(set)
        i = 0
        for line in fhandler:
            if i < 4:
                i += 1
                continue
            data = line.decode("utf-8").rstrip("\r\n").split("\t")
            protein = data[2]
            residue_mod = data[4].split('-')
            modified_protein_id = protein+'_'+data[4]
            organism = data[6]
            seq_window = data[9]
            if len(residue_mod) > 1:
                modification = modifications[residue_mod[1]]
                position = residue_mod[0][0]
                residue = ''.join(residue_mod[0][1:])
                if organism == "human":
                    #"sequence_window", "position", "Amino acid"
                    entities.add((modified_protein_id, "Modified_protein", protein,
                                  seq_window, position, residue, "PhosphositePlus"))
                    relationships[("Protein", "has_modified_site")].add((protein, modified_protein_id,
                                                                         "HAS_MODIFIED_SITE", "PhosphositePlus"))
                    relationships[("Peptide", "has_modified_site")].add((seq_window.upper(), modified_protein_id,
                                                                         "HAS_MODIFIED_SITE", "PhosphositePlus"))
                    relationships[("Modified_protein", "has_modification")].add((modified_protein_id, modification,
                                                                                 "HAS_MODIFICATION", "PhosphositePlus"))

        return entities, relationships

    def parse_kinase_substrates(self, fhandler):
        relationships = set()
        i = 0
        for line in fhandler:
            if i < 4:
                i += 1
                continue
            data = line.decode("utf-8").rstrip("\r\n").split("\t")
            kinase = data[2]
            organism = data[3]
            substrate = data[6]
            modified_protein_id = substrate+'_'+data[9]+'-p'
            if organism == "human":
                relationships.add((modified_protein_id, kinase, "IS_SUBSTRATE_OF",
                                   "NA", "CURATED", 5, "PhosphoSitePlus"))
        return relationships

    def parse_regulation_annotations(self, fhandler, mapping):
        relationships = set()
        i = 0
        for line in fhandler:
            if i < 4:
                i += 1
                continue
            data = line.decode("utf-8").rstrip("\r\n").split("\t")
            protein = data[3]
            organism = data[6]
            residue_mod = data[7].split('-')
            modified_protein_id = protein+'_'+data[7]
            functions = data[11].split('; ')
            processes = data[12].split('; ')
            pmid = data[15]
            if organism == "human":
                for process in processes:
                    if process.lower() in mapping:
                        process_code = mapping[process.lower()]
                        relationships.add((modified_protein_id, process_code, "ASSOCIATED_WITH",
                                          "CURATED", 5, "PhosphoSitePlus", pmid, "unspecified"))
                    elif process.lower().split(',')[0] in mapping:
                        process_code = mapping[process.lower().split(',')[0]]
                        relationships.add((modified_protein_id, process_code, "ASSOCIATED_WITH",
                                          "CURATED", 5, "PhosphoSitePlus", pmid, process.lower().split(',')[1]))
                    else:
                        pass
        return relationships

    def parse_disease_annotations(self, fhandler, mapping):
        relationships = set()
        i = 0
        for line in fhandler:
            if i < 4:
                i += 1
                continue
            data = line.decode("utf-8").rstrip("\r\n").split("\t")
            if len(data) > 13:
                diseases = data[0].split('; ')
                alteration = data[1]
                protein = data[4]
                organism = data[8]
                internalid = data[9]
                residue_mod = data[10].split('-')
                modified_protein_id = protein+'_'+data[10]
                pmid = data[13]
                if organism == "human":
                    for disease_name in diseases:
                        if disease_name.lower() in mapping:
                            disease_code = mapping[disease_name.lower()]
                            relationships.add((modified_protein_id, disease_code,
                                               "ASSOCIATED_WITH", "CURATED", 5, "PhosphoSitePlus", pmid))
        return relationships

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_headers = self.parse()
        entity_outputfile = os.path.join(
            self.import_directory, "psp_Modified_protein.tsv")
        self.write_entities(entities, entities_header, entity_outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(
            self.database_name, "Modified_protein", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Modified_protein",
                  self.database_name, entity_outputfile, self.updated_on))
        for entity, relationship in relationships:
            rel_header = ["START_ID", "END_ID", "TYPE", "source"]
            if entity in relationships_headers:
                rel_header = relationships_headers[entity]
            outputfile = os.path.join(self.import_directory,
                                      "psp_"+entity.lower()+"_"+relationship.lower()+".tsv")
            self.write_relationships(relationships[(entity, relationship)],
                                     rel_header, outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(self.database_name, relationship,
                                                                              len(relationships[(entity, relationship)])))
            stats.add(self._build_stats(len(relationships[(entity, relationship)]),
                                        "relationships", relationship, self.database_name, outputfile, self.updated_on))
        logger.info("Done Parsing database {}".format(self.database_name))
        return stats
