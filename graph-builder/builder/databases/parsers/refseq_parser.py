# RefSeq
import os.path
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class RefSeqParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'RefSeq'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self):
        url = self.config['refseq_url']
        ftp_dir = self.config['refseq_ftp_dir']
        entities = defaultdict(set)
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, "RefSeq")
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        headers = self.config['headerEntities']
        taxid = 9606

        if self.download:
            file_dir = self.list_ftp_directory(ftp_dir)[0]
            new_file = file_dir.split('/')[-1]+"_feature_table.txt.gz"
            url = ftp_dir + file_dir.split('/')[-1] + "/" + new_file
            self.download_db(url, directory)
            fileName = os.path.join(directory, new_file)

        if os.path.isfile(fileName):
            df = self.read_gzipped_file(fileName)
            first = True
            for line in df:
                if first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                tclass = data[1]
                assembly = data[2]
                chrom = data[5]
                geneAcc = data[6]
                start = data[7]
                end = data[8]
                strand = data[9]
                protAcc = data[10]
                name = data[13]
                symbol = data[14]

                if protAcc != "":
                    entities["Transcript"].add(
                        (protAcc, "Transcript", name, tclass, assembly, taxid))
                    if chrom != "":
                        entities["Chromosome"].add(
                            (chrom, "Chromosome", chrom, taxid))
                        relationships["LOCATED_IN"].add(
                            (protAcc, chrom, "LOCATED_IN", start, end, strand, "RefSeq"))
                    if symbol != "":
                        relationships["TRANSCRIBED_INTO"].add(
                            (symbol, protAcc, "TRANSCRIBED_INTO", "RefSeq"))
                elif geneAcc != "":
                    entities["Transcript"].add(
                        (geneAcc, "Transcript", name, tclass, assembly, taxid))
                    if chrom != "":
                        entities["Chromosome"].add(
                            (chrom, "Chromosome", chrom, taxid))
                        relationships["LOCATED_IN"].add(
                            (protAcc, chrom, "LOCATED_IN", start, end, strand, "RefSeq"))
            df.close()

        # self.remove_directory(directory)

        return (entities, relationships, headers)

    def build_stats(self):
        stats = set()
        entities, relationships, headers = self.parse()
        for entity in entities:
            header = headers[entity]
            outputfile = os.path.join(self.import_directory, entity+".tsv")
            self.write_entities(entities[entity], header, outputfile)
            logger.info("Database {} - Number of {} entities: {}".format(
                self.database_name, entity, len(entities[entity])))
            stats.add(self._build_stats(len(entities[entity]), "entity", entity,
                                        self.database_name, outputfile, self.updated_on))

        for rel in relationships:
            header = headers[rel]
            outputfile = os.path.join(
                self.import_directory, "refseq_"+rel.lower()+".tsv")
            self.write_relationships(
                relationships[rel], header, outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(
                self.database_name, rel, len(relationships[rel])))
            stats.add(self._build_stats(len(relationships[rel]), "relationships",
                                        rel, self.database_name, outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
