import re
import os.path
import logging
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser

logger = logging.getLogger(__name__)


class IntActParser(BaseParser):
    def __init__(self, import_directory, database_directory, download=True, skip=True) -> None:
        self.database_name = 'IntAct'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, download, skip)
        
    def parse(self):
        intact_dictionary = defaultdict()
        stored = set()
        relationships = set()
        header = self.config['header']
        outputfileName = "intact_interacts_with.tsv"
        regex = r"\((.*)\)"
        taxid_regex = r"\:(\d+)"
        url = self.config['intact_psimitab_url']
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)

        with open(fileName, 'r', encoding="utf-8") as idf:
            first = True
            for line in idf:
                if first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                intA = data[0].split(":")[1]
                intB = data[1].split(':')
                if len(intB)> 1:
                    intB = intB[1]
                else:
                    continue
                methodMatch = re.search(regex, data[6])
                method = methodMatch.group(1) if methodMatch else "unknown"
                publications = data[8]
                tAmatch = re.search(taxid_regex, data[9])
                tBmatch = re.search(taxid_regex, data[10])
                taxidA = ""
                taxidB = ""
                if tAmatch and tBmatch:
                    taxidA = tAmatch.group(1)
                    taxidB = tBmatch.group(1)
                itypeMatch = re.search(regex, data[11])
                itype = itypeMatch.group(1) if itypeMatch else "unknown"
                sourceMatch = re.search(regex, data[12])
                source = sourceMatch.group(1) if sourceMatch else "unknown"
                score = data[14].split(":")[1]
                if self.is_number(score):
                    score = float(score)
                else:
                    continue
                if taxidA == "9606" and taxidB == "9606":
                    if (intA, intB) in intact_dictionary:
                        intact_dictionary[(intA, intB)]['methods'].add(method)
                        intact_dictionary[(intA, intB)]['sources'].add(source)
                        intact_dictionary[(intA, intB)]['publications'].add(publications.replace('|', ','))
                        intact_dictionary[(intA, intB)]['itype'].add(itype)
                    else:
                        intact_dictionary[(intA, intB)] = {'methods': set([method]), 'sources': set([source]), 'publications': set([publications]), 'itype': set([itype]), 'score': score}
        for (intA, intB) in intact_dictionary:
            if (intA, intB, intact_dictionary[(intA, intB)]["score"]) not in stored:
                relationships.add((intA, intB, "CURATED_INTERACTS_WITH", intact_dictionary[(intA, intB)]['score'], ",".join(intact_dictionary[(intA, intB)]['itype']), ",".join(intact_dictionary[(intA, intB)]['methods']), ",".join(intact_dictionary[(intA, intB)]['sources']), ",".join(intact_dictionary[(intA, intB)]['publications'])))
                stored.add((intA, intB, intact_dictionary[(intA, intB)]["score"]))

        # self.remove_directory(directory)

        return (relationships, header, outputfileName)
    
    def build_stats(self):
        stats = set()
        relationships, header, outputfileName = self.parse()
        outputfile = os.path.join(self.import_directory, outputfileName)
        self.write_relationships(relationships, header, outputfile)
        logger.info("Database {} - Number of {} relationships: {}".format(self.database_name, "curated_interacts_with", len(relationships)))
        stats.add(self._build_stats(len(relationships), "relationships", "curated_interacts_with", 
                                    self.database_name, outputfile, self.updated_on))
        logger.info("Done Parsing database {}".format(self.database_name))
        return stats
    
    def is_number(self, s):
        """
        This function checks if given input is a float and returns True if so, and False if it is not.

        :param s: input
        :return: Boolean.
        """
        try:
            float(s)
            return True
        except ValueError:
            return False
