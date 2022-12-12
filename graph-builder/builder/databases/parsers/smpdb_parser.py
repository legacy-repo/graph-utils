import os.path
import zipfile
import pandas as pd
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class SMPDBParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'SMPDB'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        urls = self.config['smpdb_urls']
        entities = set()
        relationships = defaultdict(set)
        entities_header = self.config['pathway_header']
        relationships_headers = self.config['relationships_header']
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)

        for dataset in urls:
            url = urls[dataset]
            file_name = url.split('/')[-1]
            if self.download:
                self.download_db(url, directory)
            zipped_file = os.path.join(directory, file_name)
            with zipfile.ZipFile(zipped_file) as rf:
                if dataset == "pathway":
                    entities = self.parse_pathways(rf)
                elif dataset == "protein":
                    relationships.update(
                        self.parse_pathway_protein_relationships(rf))
                elif dataset == "metabolite":
                    relationships.update(
                        self.parse_pathway_metabolite_drug_relationships(rf))

        # self.remove_directory(directory)

        return entities, relationships, entities_header, relationships_headers

    def parse_pathways(self, fhandler):
        entities = set()
        url = self.config['linkout_url']
        organism = 9606
        for filename in fhandler.namelist():
            if not os.path.isdir(filename):
                with fhandler.open(filename) as f:
                    df = pd.read_csv(
                        f, sep=',', low_memory=False)
                    for index, row in df.iterrows():
                        identifier = row[0]
                        name = row[2]
                        description = row[3]
                        linkout = url.replace("PATHWAY", identifier)
                        entities.add((identifier, "Pathway", name,
                                     description, organism, linkout, "SMPDB"))

        return entities

    def parse_pathway_protein_relationships(self, fhandler):
        relationships = defaultdict(set)
        loc = "unspecified"
        evidence = "unspecified"
        organism = 9606
        for filename in fhandler.namelist():
            if not os.path.isdir(filename):
                with fhandler.open(filename) as f:
                    df = pd.read_csv(f, sep=',', low_memory=False)
                    for index, row in df.iterrows():
                        identifier = row[0]
                        protein = row[3]
                        if protein != '':
                            relationships[("protein", "annotated_to_pathway")].add(
                                (protein, identifier, "ANNOTATED_TO_PATHWAY", evidence, organism, loc, "SMPDB"))

        return relationships

    def parse_pathway_metabolite_drug_relationships(self, fhandler):
        relationships = defaultdict(set)
        loc = "unspecified"
        evidence = "unspecified"
        organism = 9606
        for filename in fhandler.namelist():
            if not os.path.isdir(filename):
                with fhandler.open(filename) as f:
                    df = pd.read_csv(
                        f, sep=',', low_memory=False)
                    for index, row in df.iterrows():
                        identifier = row[0]
                        metabolite = row[5]
                        drug = row[8]
                        if metabolite != '':
                            relationships[("metabolite", "annotated_to_pathway")].add(
                                (metabolite, identifier, "ANNOTATED_TO_PATHWAY", evidence, organism, loc, "SMPDB"))
                        if drug != "":
                            relationships[("drug", "annotated_to_pathway")].add(
                                (drug, identifier, "ANNOTATED_TO_PATHWAY", evidence, organism, loc, "SMPDB"))

        return relationships

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_header = self.parse()
        entity_outputfile = os.path.join(
            self.import_directory, self.database_name.lower()+"_Pathway.tsv")
        self.write_entities(entities, entities_header, entity_outputfile)
        stats.add(self._build_stats(len(entities), "entity", "Pathway",
                                    self.database_name, entity_outputfile, self.updated_on))
        for entity, relationship in relationships:
            smpdb_outputfile = os.path.join(self.import_directory,
                                            self.database_name.lower()+"_"+entity.lower()+"_"+relationship.lower()+".tsv")
            self.write_relationships(relationships[(entity, relationship)],
                                     relationships_header[entity], smpdb_outputfile)
            logger.info("Database {} - Number of {} {} relationships: {}".format(self.database_name,
                                                                                 entity, relationship, len(relationships[(entity, relationship)])))
            stats.add(self._build_stats(len(relationships[(entity, relationship)]),
                                        "relationships", relationship, self.database_name, smpdb_outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
