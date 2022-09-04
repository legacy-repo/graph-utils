import os.path
import logging
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser

logger = logging.getLogger(__name__)


class HGNCParser(BaseParser):
    def __init__(self, import_directory, database_directory, download=True, skip=True) -> None:
        self.database_name = 'HGNC'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, download, skip)

    def parse(self):
        url = self.config['hgnc_url']
        entities = set()
        directory = os.path.join(self.database_directory, "HGNC")
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        taxid = 9606
        entities_header = self.config['header']

        if self.download:
            self.download_db(url, directory)

        with open(fileName, 'r', encoding="utf-8") as df:
            first = True
            for line in df:
                if first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                geneSymbol = data[1]
                geneName = data[2]
                status = data[5]
                geneFamily = data[12]
                synonyms = data[18:23]
                transcript = data[23]
                if status != "Approved":
                    continue

                entities.add((geneSymbol, "Gene", geneName,
                             geneFamily, ",".join(synonyms), taxid))
                #relationships.add((geneSymbol, transcript, "TRANSCRIBED_INTO"))

        # self.remove_directory(directory)

        return entities, entities_header

    def build_stats(self):
        stats = set()
        entities, header = self.parse()
        outputfile = os.path.join(self.import_directory, "Gene.tsv")
        self.write_entities(entities, header, outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(
            self.database_name, "Gene", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Gene",
                  self.database_name, outputfile, self.updated_on))
        logger.info("Done Parsing database {}".format(self.database_name))
        return stats
