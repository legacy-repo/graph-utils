# The Drug Gene Interaction Database
import os.path
import logging
import verboselogs
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class DGIdbParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'DGIdb'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self):
        url = self.config['DGIdb_url']
        header = self.config['header']
        output_file = "dgidb_targets.tsv"
        drugmapping = self.get_mapping_for_entity("Drug")

        relationships = set()
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)

        with open(fileName, 'r', encoding='utf-8') as associations:
            first = True
            for line in associations:
                if first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                gene = data[0]
                source = data[3]
                interactionType = data[4] if data[4] != '' else 'unknown'
                drug = data[8].lower()
                if drug == "":
                    drug = data[7]
                    if drug == "" and data[6] != "":
                        drug = data[6]
                    else:
                        continue
                if gene != "":
                    if drug in drugmapping:
                        drug = drugmapping[drug]
                        relationships.add((drug, gene, "TARGETS", "NA", "NA", "NA",
                                           interactionType, "DGIdb: "+source))

        # self.remove_directory(directory)

        return (relationships, header, output_file)

    def build_stats(self):
        stats = set()
        relationships, header, outputfileName = self.parse()
        outputfile = os.path.join(self.import_directory, outputfileName)
        self.write_relationships(relationships, header, outputfile)
        logger.info("Database {} - Number of {} relationships: {}".format(self.database_name,
                                                                          "targets", len(relationships)))
        stats.add(self._build_stats(len(relationships), "relationships",
                                    "targets", self.database_name, outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
