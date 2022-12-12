import re
import os.path
import logging
import verboselogs
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class MutationDsParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'MutationDs'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        relationships = set()
        header = self.config['header']
        output_file_name = "mutation_curated_affects_interaction_with.tsv"
        regex = r":(\w+)\("
        url = self.config['mutations_url']
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        file_name = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)

        with open(file_name, 'r') as mf:
            first = True
            for line in mf:
                if first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                if len(data) > 12:
                    internal_id = data[0]
                    pvariant = '_'.join(data[1].split(':'))
                    effect = data[5]
                    organism = data[10]
                    interaction = data[11]
                    evidence = data[12]

                    if organism.startswith("9606"):
                        matches = re.finditer(regex, interaction)
                        for matchNum, match in enumerate(matches, start=1):
                            interactor = match.group(1)
                            relationships.add((pvariant, interactor, "CURATED_AFFECTS_INTERACTION_WITH",
                                              effect, interaction, evidence, internal_id, "Intact-MutationDs"))

        # self.remove_directory(directory)

        return (relationships, header, output_file_name)

    def build_stats(self):
        stats = set()
        relationships, header, outputfileName = self.parse()
        outputfile = os.path.join(self.import_directory, outputfileName)
        self.write_relationships(relationships, header, outputfile)
        logger.info("Database {} - Number of {} relationships: {}".format(self.database_name,
                                                                          "curated_affects_interaction_with", len(relationships)))
        stats.add(self._build_stats(len(relationships), "relationships", "curated_affects_interaction_with",
                                    self.database_name, outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
