import os.path
import pandas as pd
import logging
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser
from builder.databases.parsers.stitch_parser import STITCHParser
from builder.databases.parsers.string_parser import STRINGParser

logger = logging.getLogger(__name__)


class JensenLabMentionsParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'JensenLabMentions'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self):
        outputfileName = "Publications.tsv"
        url = self.config['db_url']
        ifile = self.config['organisms_file']
        organisms = str(self.config['organisms'])
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(os.path.join(directory, "textmining"))

        if self.download:
            self.download_db(url.replace("FILE", ifile), os.path.join(directory, "textmining"))

        ifile = os.path.join(directory, os.path.join("textmining", ifile))
        valid_pubs = self.read_valid_pubs(organisms, ifile)
        entities, header = self.parse_pmc_list(os.path.join(directory, "textmining"), valid_pubs=valid_pubs)
        num_entities = len(entities)
        outputfile = os.path.join(self.import_directory, outputfileName)
        self.write_entities(entities, header, outputfile)
        entities = None

        for qtype in self.config['db_mentions_types']:
            self.parse_mentions(directory, qtype)

        # self.remove_directory(os.path.join(directory, "textmining"))

        return (num_entities, outputfile)


    def read_valid_pubs(self, organisms, organisms_file):
        pubs = set()
        with open(organisms_file, 'r') as idbf:
            for line in idbf:
                data = line.rstrip('\r\n').split('\t')
                if str(data[0]) in organisms:
                    pubs.update(set(data[1].split(" ")))
        return list(pubs)


    def parse_pmc_list(self, directory, valid_pubs=None):
        url = self.config['PMC_db_url']
        plinkout = self.config['pubmed_linkout']
        entities = set()
        fileName = os.path.join(directory, url.split('/')[-1])

        if self.download:
            self.download_db(url, directory)

        entities = pd.read_csv(fileName, sep=',', dtype=str, compression='gzip', low_memory=False)
        entities = entities[config['PMC_fields']]
        entities = entities[entities.iloc[:, 0].notnull()]
        entities = entities.set_index(list(entities.columns)[0])
        if valid_pubs is not None:
            valid_pubs = set(entities.index).intersection(valid_pubs)
            entities = entities.loc[list(valid_pubs)]

        entities['linkout'] = [plinkout.replace("PUBMEDID", str(int(pubmedid))) for pubmedid in list(entities.index)]
        entities.index.names = ['ID']
        entities['TYPE'] = 'Publication'
        entities = entities.reset_index()
        header = [c.replace(' ', '_').lower() if c not in ['ID', 'TYPE'] else c for c in list(entities.columns)]
        entities = entities.replace('\\\\', '', regex=True)
        entities = list(entities.itertuples(index=False))

        return entities, header

    def parse_mentions(self, directory, qtype):
        url = self.config['db_url']
        ifile = self.config['db_mentions_files'][qtype]
        if qtype == "9606":
            string_parser = STRINGParser(self.import_directory, self.database_directory, 
                                         self.download, self.skip)
            mapping = string_parser.get_string_mapping()
        elif qtype == "-1":
            stitch_parser = STITCHParser(self.import_directory, self.database_directory, 
                                         self.download, self.skip)
            source = self.builder_config['database']['sources']["Drug"]
            mapping = stitch_parser.get_string_mapping(source=source)

        filters = []
        if qtype in config['db_mentions_filters']:
            filters = config['db_mentions_filters'][qtype]
        entity1, entity2 = config['db_mentions_types'][qtype]
        outputfile = os.path.join(self.import_directory, entity1 + "_" + entity2 + "_mentioned_in_publication.tsv")

        if self.download:
            self.download_db(url.replace("FILE", ifile), os.path.join(directory, "textmining"))
        ifile = os.path.join(directory, os.path.join("textmining", ifile))
        with open(outputfile, 'w') as f:
            f.write("START_ID\tEND_ID\tTYPE\n")
            with open(ifile, 'r') as idbf:
                for line in idbf:
                    data = line.rstrip("\r\n").split('\t')
                    id1 = data[0]
                    pubmedids = data[1].split(" ")
                    ident = []
                    if qtype == "9606":
                        id1 = "9606."+id1
                        if id1 in mapping:
                            ident = mapping[id1]
                    elif qtype == "-1":
                        if id1 in mapping:
                            ident = mapping[id1]
                    elif qtype == "-26":
                        if id1.startswith("DOID"):
                            ident = [id1]
                    else:
                        ident = [id1]

                    for i in ident:
                        if i not in filters:
                            aux = pd.DataFrame(data={"Pubmedids": list(pubmedids)})
                            aux["START_ID"] = i
                            aux["TYPE"] = "MENTIONED_IN_PUBLICATION"
                            aux.to_csv(path_or_buf=f, sep='\t', header=False, index=False, quotechar='"', line_terminator='\n', escapechar='\\')
                            aux = None

    def build_stats(self):
        stats = set()
        num_entities, outputfile = self.parse()
        logger.info("Database {} - Number of {} entities: {}".format(self.database_name, "Publication", num_entities))
        stats.add(self._build_stats(num_entities, "entity", "Publication", 
                                    self.database_name, outputfile, self.updated_on))
        logger.info("Done Parsing database {}".format(self.database_name))
        return stats
