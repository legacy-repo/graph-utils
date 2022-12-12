import os.path
import zipfile
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class CORUMParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'CORUM'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        entities = set()
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)

        database_url = self.config['database_url']
        entities_header = self.config['entities_header']
        relationships_headers = self.config['relationships_headers']
        zipped_fileName = os.path.join(directory, database_url.split('/')[-1])
        fileName = '.'.join(database_url.split('/')[-1].split('.')[0:2])
        if self.download:
            self.download_db(database_url, directory)
        names = set()
        first = True
        with zipfile.ZipFile(zipped_fileName) as z:
            with z.open(fileName) as f:
                for line in f:
                    if first:
                        first = False
                        continue
                    data = line.decode("utf-8").rstrip("\r\n").split("\t")
                    identifier = data[0]
                    name = data[1]
                    organism = data[2]
                    synonyms = data[3].split(';') \
                        if data[3] != "None" else [""]
                    cell_lines = data[4].join(';')
                    subunits = data[5].split(';')
                    evidences = data[7].split(';')
                    processes = data[8].split(';')
                    pubmedid = data[14]

                    if organism == "Human":
                        # ID name organism synonyms source
                        if name not in names:
                            entities.add((identifier, name, "9606",
                                          ",".join(synonyms), "CORUM"))
                            names.add(name)
                        for subunit in subunits:
                            # START_ID END_ID type cell_lines evidences publication source
                            relationships[("Protein", "is_subunit_of")].add(
                                (subunit, identifier, "IS_SUBUNIT_OF", ",".join(cell_lines), ",".join(evidences), pubmedid, "CORUM"))
                        for process in processes:
                            # START_ID END_ID type evidence_type score source
                            relationships["Biological_process", "associated_with"].add(
                                (identifier, process, "ASSOCIATED_WITH", "CURATED", 5, "CORUM"))

        # self.remove_directory(directory)

        return entities, relationships, entities_header, relationships_headers

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_headers = self.parse()
        entity_outputfile = os.path.join(self.import_directory, "Complex.tsv")
        self.write_entities(entities, entities_header, entity_outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(self.database_name,
                                                                     "Complex", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Complex",
                                    self.database_name, entity_outputfile, self.updated_on))
        for entity, relationship in relationships:
            corum_outputfile = os.path.join(self.import_directory,
                                            self.database_name.lower()+"_"+entity.lower()+"_"+relationship.lower()+".tsv")
            self.write_relationships(relationships[(entity, relationship)],
                                     relationships_headers[entity], corum_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(self.database_name,
                                                                              relationship, len(relationships[(entity, relationship)])))
            stats.add(self._build_stats(len(relationships[(entity, relationship)]), "relationships",
                                        relationship, self.database_name, corum_outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
